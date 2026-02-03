use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;
use std::thread;

use clap::Parser;
use crossbeam_channel::{bounded, Receiver};
use hashbrown::HashSet;
use nohash_hasher::BuildNoHashHasher;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use xxhash_rust::xxh3::xxh3_64;

type NoHashSet64 = HashSet<u64, BuildNoHashHasher<u64>>;
type NoHashSet32 = HashSet<u32, BuildNoHashHasher<u32>>;

const K_MIX: u64 = 0x9e3779b97f4a7c15;

#[derive(Parser)]
#[command(about = "K-mer membership using sharded hash tables")]
struct Args {
    /// Reference FASTA
    reference: PathBuf,
    /// Query FASTA
    query: PathBuf,
    /// K-mer size
    #[arg(short = 'k', long, default_value_t = 31)]
    k: usize,
    /// Number of shards (power of two)
    #[arg(long, default_value_t = 256)]
    shards: usize,
    /// Worker threads (0 = auto)
    #[arg(long, default_value_t = 0)]
    threads: usize,
    /// Use 32-bit hashes (more collisions, less memory)
    #[arg(long)]
    hash32: bool,
    /// Number of gap positions in the k-mer
    #[arg(long, default_value_t = 0)]
    gaps: usize,
    /// Optional random seed for gap placement
    #[arg(long)]
    seed: Option<u64>,
    /// Number of different gap patterns to index
    #[arg(long, default_value_t = 1)]
    gap_patterns: usize,
}

struct GappedMask {
    k: usize,
    include: Vec<usize>,
    gaps: usize,
    gap_positions: Vec<usize>,
}

fn main() -> io::Result<()> {
    let args = Args::parse();
    if args.k == 0 {
        eprintln!("Error: --k must be >= 1");
        std::process::exit(2);
    }
    if args.gaps >= args.k {
        eprintln!("Error: --gaps must be < k");
        std::process::exit(2);
    }
    if args.gaps > 0 {
        let allowed = args.k.saturating_sub(8);
        if args.gaps > allowed {
            eprintln!(
                "Error: --gaps must be <= k-8 (no gaps allowed in the first 8 positions)"
            );
            std::process::exit(2);
        }
    }
    if args.gap_patterns == 0 {
        eprintln!("Error: --gap-patterns must be >= 1");
        std::process::exit(2);
    }
    if !args.shards.is_power_of_two() {
        eprintln!("Error: --shards must be a power of two");
        std::process::exit(2);
    }

    let threads = if args.threads == 0 {
        num_cpus::get().max(1)
    } else {
        args.threads.max(1)
    };
    let shard_bits = args.shards.trailing_zeros();
    let hash_bits = if args.hash32 { 32 } else { 64 };
    if shard_bits > hash_bits {
        eprintln!("Error: --shards uses more bits than the selected hash width");
        std::process::exit(2);
    }

    let masks = Arc::new(make_gapped_masks(
        args.k,
        args.gaps,
        args.gap_patterns,
        args.seed,
    ));
    if masks[0].gaps > 0 {
        eprintln!(
            "Using gapped k-mers: k={}, gaps={}, patterns={}",
            masks[0].k,
            masks[0].gaps,
            masks.len()
        );
        for (i, mask) in masks.iter().enumerate() {
            eprintln!("gap_pattern_{}\t{:?}", i, mask.gap_positions);
        }
    }

    if args.hash32 {
        let table =
            build_table32(&args.reference, Arc::clone(&masks), args.shards, shard_bits, threads)?;
        let table = Arc::new(table);
        let (hits, seqs, total_kmers, total_seqs) =
            query_table32(&args.query, table, masks, shard_bits, threads)?;
        print_stats(hits, total_kmers, seqs, total_seqs);
    } else {
        let table =
            build_table64(&args.reference, Arc::clone(&masks), args.shards, shard_bits, threads)?;
        let table = Arc::new(table);
        let (hits, seqs, total_kmers, total_seqs) =
            query_table64(&args.query, table, masks, shard_bits, threads)?;
        print_stats(hits, total_kmers, seqs, total_seqs);
    }

    Ok(())
}

fn build_table64(
    path: &Path,
    masks: Arc<Vec<GappedMask>>,
    shards: usize,
    shard_bits: u32,
    threads: usize,
) -> io::Result<Vec<NoHashSet64>> {
    let (tx, rx) = bounded::<Vec<u8>>(threads * 4);
    let mut handles = Vec::with_capacity(threads);
    for _ in 0..threads {
        let rx = rx.clone();
        let masks = Arc::clone(&masks);
        let handle = thread::spawn(move || build_worker64(rx, &masks, shards, shard_bits));
        handles.push(handle);
    }
    drop(rx);

    read_fasta_sequences(path, |seq| {
        if !seq.is_empty() {
            let _ = tx.send(seq);
        }
    })?;
    drop(tx);

    let mut global = make_shards64(shards);
    for handle in handles {
        let local = handle.join().expect("worker thread panicked");
        merge_shards64(&mut global, local);
    }
    Ok(global)
}

fn build_table32(
    path: &Path,
    masks: Arc<Vec<GappedMask>>,
    shards: usize,
    shard_bits: u32,
    threads: usize,
) -> io::Result<Vec<NoHashSet32>> {
    let (tx, rx) = bounded::<Vec<u8>>(threads * 4);
    let mut handles = Vec::with_capacity(threads);
    for _ in 0..threads {
        let rx = rx.clone();
        let masks = Arc::clone(&masks);
        let handle = thread::spawn(move || build_worker32(rx, &masks, shards, shard_bits));
        handles.push(handle);
    }
    drop(rx);

    read_fasta_sequences(path, |seq| {
        if !seq.is_empty() {
            let _ = tx.send(seq);
        }
    })?;
    drop(tx);

    let mut global = make_shards32(shards);
    for handle in handles {
        let local = handle.join().expect("worker thread panicked");
        merge_shards32(&mut global, local);
    }
    Ok(global)
}

fn query_table64(
    path: &Path,
    table: Arc<Vec<NoHashSet64>>,
    masks: Arc<Vec<GappedMask>>,
    shard_bits: u32,
    threads: usize,
) -> io::Result<(u64, u64, u64, u64)> {
    let (tx, rx) = bounded::<Vec<u8>>(threads * 4);
    let total_hits = Arc::new(AtomicU64::new(0));
    let seq_hits = Arc::new(AtomicU64::new(0));
    let total_kmers = Arc::new(AtomicU64::new(0));
    let total_seqs = Arc::new(AtomicU64::new(0));

    let mut handles = Vec::with_capacity(threads);
    for _ in 0..threads {
        let rx = rx.clone();
        let table = Arc::clone(&table);
        let masks = Arc::clone(&masks);
        let total_hits = Arc::clone(&total_hits);
        let seq_hits = Arc::clone(&seq_hits);
        let total_kmers = Arc::clone(&total_kmers);
        let total_seqs = Arc::clone(&total_seqs);
        let handle = thread::spawn(move || {
            query_worker64(
                rx,
                &table,
                &masks,
                shard_bits,
                &total_hits,
                &seq_hits,
                &total_kmers,
                &total_seqs,
            )
        });
        handles.push(handle);
    }
    drop(rx);

    read_fasta_sequences(path, |seq| {
        if !seq.is_empty() {
            let _ = tx.send(seq);
        }
    })?;
    drop(tx);

    for handle in handles {
        handle.join().expect("worker thread panicked");
    }

    Ok((
        total_hits.load(Ordering::Relaxed),
        seq_hits.load(Ordering::Relaxed),
        total_kmers.load(Ordering::Relaxed),
        total_seqs.load(Ordering::Relaxed),
    ))
}

fn query_table32(
    path: &Path,
    table: Arc<Vec<NoHashSet32>>,
    masks: Arc<Vec<GappedMask>>,
    shard_bits: u32,
    threads: usize,
) -> io::Result<(u64, u64, u64, u64)> {
    let (tx, rx) = bounded::<Vec<u8>>(threads * 4);
    let total_hits = Arc::new(AtomicU64::new(0));
    let seq_hits = Arc::new(AtomicU64::new(0));
    let total_kmers = Arc::new(AtomicU64::new(0));
    let total_seqs = Arc::new(AtomicU64::new(0));

    let mut handles = Vec::with_capacity(threads);
    for _ in 0..threads {
        let rx = rx.clone();
        let table = Arc::clone(&table);
        let masks = Arc::clone(&masks);
        let total_hits = Arc::clone(&total_hits);
        let seq_hits = Arc::clone(&seq_hits);
        let total_kmers = Arc::clone(&total_kmers);
        let total_seqs = Arc::clone(&total_seqs);
        let handle = thread::spawn(move || {
            query_worker32(
                rx,
                &table,
                &masks,
                shard_bits,
                &total_hits,
                &seq_hits,
                &total_kmers,
                &total_seqs,
            )
        });
        handles.push(handle);
    }
    drop(rx);

    read_fasta_sequences(path, |seq| {
        if !seq.is_empty() {
            let _ = tx.send(seq);
        }
    })?;
    drop(tx);

    for handle in handles {
        handle.join().expect("worker thread panicked");
    }

    Ok((
        total_hits.load(Ordering::Relaxed),
        seq_hits.load(Ordering::Relaxed),
        total_kmers.load(Ordering::Relaxed),
        total_seqs.load(Ordering::Relaxed),
    ))
}

fn build_worker64(
    rx: Receiver<Vec<u8>>,
    masks: &[GappedMask],
    shards: usize,
    shard_bits: u32,
) -> Vec<NoHashSet64> {
    let mut local = make_shards64(shards);
    if masks.len() == 1 {
        let mask = &masks[0];
        for seq in rx {
            for_each_hash64(&seq, mask, |h| {
                let idx = shard_index64(h, shard_bits);
                local[idx].insert(h);
            });
        }
    } else {
        for seq in rx {
            for_each_hash_multi64(&seq, masks, |h| {
                let idx = shard_index64(h, shard_bits);
                local[idx].insert(h);
            });
        }
    }
    local
}

fn build_worker32(
    rx: Receiver<Vec<u8>>,
    masks: &[GappedMask],
    shards: usize,
    shard_bits: u32,
) -> Vec<NoHashSet32> {
    let mut local = make_shards32(shards);
    if masks.len() == 1 {
        let mask = &masks[0];
        for seq in rx {
            for_each_hash32(&seq, mask, |h| {
                let idx = shard_index32(h, shard_bits);
                local[idx].insert(h);
            });
        }
    } else {
        for seq in rx {
            for_each_hash_multi32(&seq, masks, |h| {
                let idx = shard_index32(h, shard_bits);
                local[idx].insert(h);
            });
        }
    }
    local
}

fn query_worker64(
    rx: Receiver<Vec<u8>>,
    table: &[NoHashSet64],
    masks: &[GappedMask],
    shard_bits: u32,
    total_hits: &AtomicU64,
    seq_hits: &AtomicU64,
    total_kmers: &AtomicU64,
    total_seqs: &AtomicU64,
) {
    let mut local_hits = 0u64;
    let mut local_seq_hits = 0u64;
    let mut local_total_kmers = 0u64;
    let mut local_total_seqs = 0u64;
    if masks.len() == 1 {
        let mask = &masks[0];
        for seq in rx {
            local_total_seqs += 1;
            let mut hits = 0u64;
            let mut total = 0u64;
            for_each_hash64(&seq, mask, |h| {
                total += 1;
                let idx = shard_index64(h, shard_bits);
                if table[idx].contains(&h) {
                    hits += 1;
                }
            });
            local_total_kmers += total;
            if hits > 0 {
                local_hits += hits;
                local_seq_hits += 1;
            }
        }
    } else {
        for seq in rx {
            local_total_seqs += 1;
            let (hits, total) = query_multi64(&seq, masks, table, shard_bits);
            local_total_kmers += total;
            if hits > 0 {
                local_hits += hits;
                local_seq_hits += 1;
            }
        }
    }
    total_hits.fetch_add(local_hits, Ordering::Relaxed);
    seq_hits.fetch_add(local_seq_hits, Ordering::Relaxed);
    total_kmers.fetch_add(local_total_kmers, Ordering::Relaxed);
    total_seqs.fetch_add(local_total_seqs, Ordering::Relaxed);
}

fn query_worker32(
    rx: Receiver<Vec<u8>>,
    table: &[NoHashSet32],
    masks: &[GappedMask],
    shard_bits: u32,
    total_hits: &AtomicU64,
    seq_hits: &AtomicU64,
    total_kmers: &AtomicU64,
    total_seqs: &AtomicU64,
) {
    let mut local_hits = 0u64;
    let mut local_seq_hits = 0u64;
    let mut local_total_kmers = 0u64;
    let mut local_total_seqs = 0u64;
    if masks.len() == 1 {
        let mask = &masks[0];
        for seq in rx {
            local_total_seqs += 1;
            let mut hits = 0u64;
            let mut total = 0u64;
            for_each_hash32(&seq, mask, |h| {
                total += 1;
                let idx = shard_index32(h, shard_bits);
                if table[idx].contains(&h) {
                    hits += 1;
                }
            });
            local_total_kmers += total;
            if hits > 0 {
                local_hits += hits;
                local_seq_hits += 1;
            }
        }
    } else {
        for seq in rx {
            local_total_seqs += 1;
            let (hits, total) = query_multi32(&seq, masks, table, shard_bits);
            local_total_kmers += total;
            if hits > 0 {
                local_hits += hits;
                local_seq_hits += 1;
            }
        }
    }
    total_hits.fetch_add(local_hits, Ordering::Relaxed);
    seq_hits.fetch_add(local_seq_hits, Ordering::Relaxed);
    total_kmers.fetch_add(local_total_kmers, Ordering::Relaxed);
    total_seqs.fetch_add(local_total_seqs, Ordering::Relaxed);
}

fn read_fasta_sequences<P: AsRef<Path>>(
    path: P,
    mut on_seq: impl FnMut(Vec<u8>),
) -> io::Result<()> {
    let file = File::open(path)?;
    let mut reader = BufReader::with_capacity(1024 * 1024, file);
    let mut line = Vec::new();
    let mut seq = Vec::new();

    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line)?;
        if n == 0 {
            break;
        }
        while matches!(line.last(), Some(b'\n' | b'\r')) {
            line.pop();
        }
        if line.is_empty() {
            continue;
        }
        if line[0] == b'>' {
            if !seq.is_empty() {
                on_seq(std::mem::take(&mut seq));
            }
            continue;
        }
        seq.extend_from_slice(&line);
    }
    if !seq.is_empty() {
        on_seq(seq);
    }
    Ok(())
}

#[inline(always)]
fn base_to_bits(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' | b'U' | b'u' => 3,
        _ => 0xFF,
    }
}

#[inline(always)]
fn normalize_base(b: u8) -> Option<u8> {
    match b {
        b'A' | b'a' => Some(b'A'),
        b'C' | b'c' => Some(b'C'),
        b'G' | b'g' => Some(b'G'),
        b'T' | b't' | b'U' | b'u' => Some(b'T'),
        _ => None,
    }
}

#[inline(always)]
fn complement_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

fn fill_canonical_window(seq: &[u8], start: usize, k: usize, canon: &mut [u8]) -> bool {
    let mut order = 0i8;
    for i in 0..k {
        let f = match normalize_base(seq[start + i]) {
            Some(b) => b,
            None => return false,
        };
        let r_base = match normalize_base(seq[start + k - 1 - i]) {
            Some(b) => b,
            None => return false,
        };
        let r = complement_base(r_base);
        canon[i] = f;
        if order == 0 {
            if f < r {
                order = -1;
            } else if f > r {
                order = 1;
            }
        }
    }
    if order == 1 {
        for i in 0..k {
            let b = match normalize_base(seq[start + k - 1 - i]) {
                Some(x) => x,
                None => return false,
            };
            canon[i] = complement_base(b);
        }
    }
    true
}

fn for_each_canonical_window(seq: &[u8], k: usize, mut f: impl FnMut(&[u8])) {
    if k == 0 || seq.len() < k {
        return;
    }
    let mut canon = vec![0u8; k];
    for start in 0..=seq.len() - k {
        if fill_canonical_window(seq, start, k, &mut canon) {
            f(&canon);
        }
    }
}

#[inline(always)]
fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^ (z >> 31)
}

fn for_each_hash64(seq: &[u8], mask: &GappedMask, mut f: impl FnMut(u64)) {
    let k = mask.k;
    if k == 0 {
        return;
    }
    if mask.gaps == 0 && k <= 32 {
        let bitmask = if k == 32 {
            u64::MAX
        } else {
            (1u64 << (2 * k)) - 1
        };
        let mut fwd = 0u64;
        let mut rev = 0u64;
        let mut valid = 0usize;
        for &b in seq {
            let bits = base_to_bits(b);
            if bits == 0xFF {
                valid = 0;
                fwd = 0;
                rev = 0;
                continue;
            }
            fwd = (fwd << 2) | bits as u64;
            if k < 32 {
                fwd &= bitmask;
            }
            let comp = (bits ^ 0b11) as u64;
            rev = (rev >> 2) | (comp << (2 * (k - 1)));
            valid += 1;
            if valid >= k {
                let canon = if fwd < rev { fwd } else { rev };
                let h = mix64(canon ^ K_MIX.wrapping_mul(k as u64));
                f(h);
            }
        }
        return;
    }
    if mask.gaps == 0 {
        for_each_canonical_window(seq, k, |canon| {
            let h = xxh3_64(canon);
            f(h);
        });
        return;
    }

    let include = &mask.include;
    let m = include.len();
    if m == 0 || seq.len() < k {
        return;
    }
    if m <= 32 {
        for_each_canonical_window(seq, k, |canon| {
            let mut val = 0u64;
            for &pos in include {
                let bits = base_to_bits(canon[pos]);
                val = (val << 2) | bits as u64;
            }
            let h = mix64(val ^ K_MIX.wrapping_mul(m as u64));
            f(h);
        });
    } else {
        let mut buf = vec![0u8; m];
        for_each_canonical_window(seq, k, |canon| {
            for (j, &pos) in include.iter().enumerate() {
                buf[j] = canon[pos];
            }
            let h = xxh3_64(&buf);
            f(h);
        });
    }
}

fn for_each_hash32(seq: &[u8], mask: &GappedMask, mut f: impl FnMut(u32)) {
    let k = mask.k;
    if k == 0 {
        return;
    }
    if mask.gaps == 0 && k <= 32 {
        let bitmask = if k == 32 {
            u64::MAX
        } else {
            (1u64 << (2 * k)) - 1
        };
        let mut fwd = 0u64;
        let mut rev = 0u64;
        let mut valid = 0usize;
        for &b in seq {
            let bits = base_to_bits(b);
            if bits == 0xFF {
                valid = 0;
                fwd = 0;
                rev = 0;
                continue;
            }
            fwd = (fwd << 2) | bits as u64;
            if k < 32 {
                fwd &= bitmask;
            }
            let comp = (bits ^ 0b11) as u64;
            rev = (rev >> 2) | (comp << (2 * (k - 1)));
            valid += 1;
            if valid >= k {
                let canon = if fwd < rev { fwd } else { rev };
                let h = mix64(canon ^ K_MIX.wrapping_mul(k as u64));
                f(h as u32);
            }
        }
        return;
    }
    if mask.gaps == 0 {
        for_each_canonical_window(seq, k, |canon| {
            let h = xxh3_64(canon);
            f(h as u32);
        });
        return;
    }

    let include = &mask.include;
    let m = include.len();
    if m == 0 || seq.len() < k {
        return;
    }
    if m <= 32 {
        for_each_canonical_window(seq, k, |canon| {
            let mut val = 0u64;
            for &pos in include {
                let bits = base_to_bits(canon[pos]);
                val = (val << 2) | bits as u64;
            }
            let h = mix64(val ^ K_MIX.wrapping_mul(m as u64));
            f(h as u32);
        });
    } else {
        let mut buf = vec![0u8; m];
        for_each_canonical_window(seq, k, |canon| {
            for (j, &pos) in include.iter().enumerate() {
                buf[j] = canon[pos];
            }
            let h = xxh3_64(&buf);
            f(h as u32);
        });
    }
}

fn for_each_hash_multi64(seq: &[u8], masks: &[GappedMask], mut f: impl FnMut(u64)) {
    if masks.is_empty() {
        return;
    }
    let k = masks[0].k;
    if k == 0 || seq.len() < k {
        return;
    }
    let mut buffers: Vec<Vec<u8>> = masks
        .iter()
        .map(|m| {
            let mlen = m.include.len();
            if mlen > 32 {
                vec![0u8; mlen]
            } else {
                Vec::new()
            }
        })
        .collect();

    for_each_canonical_window(seq, k, |canon| {
        for (mask, buf) in masks.iter().zip(buffers.iter_mut()) {
            let include = &mask.include;
            let mlen = include.len();
            if mlen <= 32 {
                let mut val = 0u64;
                for &pos in include {
                    let bits = base_to_bits(canon[pos]);
                    val = (val << 2) | bits as u64;
                }
                let h = mix64(val ^ K_MIX.wrapping_mul(mlen as u64));
                f(h);
            } else {
                for (j, &pos) in include.iter().enumerate() {
                    buf[j] = canon[pos];
                }
                let h = xxh3_64(&buf);
                f(h);
            }
        }
    });
}

fn for_each_hash_multi32(seq: &[u8], masks: &[GappedMask], mut f: impl FnMut(u32)) {
    if masks.is_empty() {
        return;
    }
    let k = masks[0].k;
    if k == 0 || seq.len() < k {
        return;
    }
    let mut buffers: Vec<Vec<u8>> = masks
        .iter()
        .map(|m| {
            let mlen = m.include.len();
            if mlen > 32 {
                vec![0u8; mlen]
            } else {
                Vec::new()
            }
        })
        .collect();

    for_each_canonical_window(seq, k, |canon| {
        for (mask, buf) in masks.iter().zip(buffers.iter_mut()) {
            let include = &mask.include;
            let mlen = include.len();
            if mlen <= 32 {
                let mut val = 0u64;
                for &pos in include {
                    let bits = base_to_bits(canon[pos]);
                    val = (val << 2) | bits as u64;
                }
                let h = mix64(val ^ K_MIX.wrapping_mul(mlen as u64));
                f(h as u32);
            } else {
                for (j, &pos) in include.iter().enumerate() {
                    buf[j] = canon[pos];
                }
                let h = xxh3_64(&buf);
                f(h as u32);
            }
        }
    });
}

fn query_multi64(
    seq: &[u8],
    masks: &[GappedMask],
    table: &[NoHashSet64],
    shard_bits: u32,
) -> (u64, u64) {
    if masks.is_empty() {
        return (0, 0);
    }
    let k = masks[0].k;
    if k == 0 || seq.len() < k {
        return (0, 0);
    }
    let mut buffers: Vec<Vec<u8>> = masks
        .iter()
        .map(|m| {
            let mlen = m.include.len();
            if mlen > 32 {
                vec![0u8; mlen]
            } else {
                Vec::new()
            }
        })
        .collect();

    let mut hits = 0u64;
    let mut total = 0u64;
    for_each_canonical_window(seq, k, |canon| {
        total += 1;
        let mut matched = false;
        for (mask, buf) in masks.iter().zip(buffers.iter_mut()) {
            let include = &mask.include;
            let mlen = include.len();
            if mlen <= 32 {
                let mut val = 0u64;
                for &pos in include {
                    let bits = base_to_bits(canon[pos]);
                    val = (val << 2) | bits as u64;
                }
                let h = mix64(val ^ K_MIX.wrapping_mul(mlen as u64));
                let idx = shard_index64(h, shard_bits);
                if table[idx].contains(&h) {
                    matched = true;
                    break;
                }
            } else {
                for (j, &pos) in include.iter().enumerate() {
                    buf[j] = canon[pos];
                }
                let h = xxh3_64(&buf);
                let idx = shard_index64(h, shard_bits);
                if table[idx].contains(&h) {
                    matched = true;
                    break;
                }
            }
        }
        if matched {
            hits += 1;
        }
    });
    (hits, total)
}

fn query_multi32(
    seq: &[u8],
    masks: &[GappedMask],
    table: &[NoHashSet32],
    shard_bits: u32,
) -> (u64, u64) {
    if masks.is_empty() {
        return (0, 0);
    }
    let k = masks[0].k;
    if k == 0 || seq.len() < k {
        return (0, 0);
    }
    let mut buffers: Vec<Vec<u8>> = masks
        .iter()
        .map(|m| {
            let mlen = m.include.len();
            if mlen > 32 {
                vec![0u8; mlen]
            } else {
                Vec::new()
            }
        })
        .collect();

    let mut hits = 0u64;
    let mut total = 0u64;
    for_each_canonical_window(seq, k, |canon| {
        total += 1;
        let mut matched = false;
        for (mask, buf) in masks.iter().zip(buffers.iter_mut()) {
            let include = &mask.include;
            let mlen = include.len();
            if mlen <= 32 {
                let mut val = 0u64;
                for &pos in include {
                    let bits = base_to_bits(canon[pos]);
                    val = (val << 2) | bits as u64;
                }
                let h = mix64(val ^ K_MIX.wrapping_mul(mlen as u64)) as u32;
                let idx = shard_index32(h, shard_bits);
                if table[idx].contains(&h) {
                    matched = true;
                    break;
                }
            } else {
                for (j, &pos) in include.iter().enumerate() {
                    buf[j] = canon[pos];
                }
                let h = xxh3_64(&buf) as u32;
                let idx = shard_index32(h, shard_bits);
                if table[idx].contains(&h) {
                    matched = true;
                    break;
                }
            }
        }
        if matched {
            hits += 1;
        }
    });
    (hits, total)
}

#[inline(always)]
fn shard_index64(hash: u64, shard_bits: u32) -> usize {
    if shard_bits == 0 {
        0
    } else {
        (hash >> (64 - shard_bits)) as usize
    }
}

#[inline(always)]
fn shard_index32(hash: u32, shard_bits: u32) -> usize {
    if shard_bits == 0 {
        0
    } else {
        (hash >> (32 - shard_bits)) as usize
    }
}

fn make_shards64(shards: usize) -> Vec<NoHashSet64> {
    let mut v = Vec::with_capacity(shards);
    for _ in 0..shards {
        v.push(HashSet::with_capacity_and_hasher(
            0,
            BuildNoHashHasher::default(),
        ));
    }
    v
}

fn make_shards32(shards: usize) -> Vec<NoHashSet32> {
    let mut v = Vec::with_capacity(shards);
    for _ in 0..shards {
        v.push(HashSet::with_capacity_and_hasher(
            0,
            BuildNoHashHasher::default(),
        ));
    }
    v
}

fn merge_shards64(global: &mut [NoHashSet64], mut local: Vec<NoHashSet64>) {
    for (g, l) in global.iter_mut().zip(local.iter_mut()) {
        g.reserve(l.len());
        g.extend(l.drain());
    }
}

fn merge_shards32(global: &mut [NoHashSet32], mut local: Vec<NoHashSet32>) {
    for (g, l) in global.iter_mut().zip(local.iter_mut()) {
        g.reserve(l.len());
        g.extend(l.drain());
    }
}

fn make_gapped_masks(k: usize, gaps: usize, patterns: usize, seed: Option<u64>) -> Vec<GappedMask> {
    if gaps == 0 {
        return vec![GappedMask {
            k,
            include: (0..k).collect(),
            gaps: 0,
            gap_positions: Vec::new(),
        }];
    }
    let mut rng = match seed {
        Some(s) => StdRng::seed_from_u64(s),
        None => StdRng::from_entropy(),
    };
    let mut masks = Vec::with_capacity(patterns);
    for _ in 0..patterns {
        let mut order: Vec<usize> = (8..k).collect();
        order.shuffle(&mut rng);
        let mut gap_positions: Vec<usize> = order.into_iter().take(gaps).collect();
        gap_positions.sort_unstable();
        let mut is_gap = vec![false; k];
        for idx in &gap_positions {
            is_gap[*idx] = true;
        }
        let include: Vec<usize> = (0..k).filter(|&i| !is_gap[i]).collect();
        masks.push(GappedMask {
            k,
            include,
            gaps,
            gap_positions,
        });
    }
    masks
}

fn print_stats(hits: u64, total_kmers: u64, seq_hits: u64, total_seqs: u64) {
    let kmer_pct = if total_kmers == 0 {
        0.0
    } else {
        (hits as f64) * 100.0 / (total_kmers as f64)
    };
    let seq_pct = if total_seqs == 0 {
        0.0
    } else {
        (seq_hits as f64) * 100.0 / (total_seqs as f64)
    };

    println!(
        "query_kmer_hits\t{} ({:.2}% of {})",
        format_u64_commas(hits),
        kmer_pct,
        format_u64_commas(total_kmers)
    );
    println!(
        "query_sequences_with_hit\t{} ({:.2}% of {})",
        format_u64_commas(seq_hits),
        seq_pct,
        format_u64_commas(total_seqs)
    );
}

fn format_u64_commas(mut n: u64) -> String {
    if n == 0 {
        return "0".to_string();
    }
    let mut out = Vec::new();
    let mut count = 0;
    while n > 0 {
        if count == 3 {
            out.push(b',');
            count = 0;
        }
        out.push(b'0' + (n % 10) as u8);
        n /= 10;
        count += 1;
    }
    out.reverse();
    String::from_utf8(out).unwrap_or_else(|_| "0".to_string())
}
