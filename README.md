# Sid

Fast k-mer membership for FASTA files with sharded hash tables, optional gapped k-mers, and canonicalization (min of forward and reverse-complement).

## Build

```bash
cargo build --release
```

The binary will be at `target/release/Sid`.

## Usage

```bash
./target/release/Sid <reference.fa> <query.fa> [options]
```

The program:
- Indexes all k-mers from the reference into sharded hash tables.
- Scans the query and counts how many of its k-mers are present.
- Reports total matching k-mer windows and the number of query sequences with at least one hit.

## Options

- `-k, --k <INT>`: K-mer length. Default: `31`.
- `--shards <INT>`: Number of hash table shards. Must be a power of two. Default: `256`.
- `--threads <INT>`: Worker threads. `0` means auto-detect (number of CPUs). Default: `0`.
- `--hash32`: Use 32-bit hashes instead of 64-bit (lower memory, higher collision risk).
- `--gaps <INT>`: Number of gap positions in each k-mer. Default: `0`. Constraints: `gaps < k`, and when `gaps > 0`, no gaps are allowed in the first 8 positions, so `gaps <= k-8`.
- `--gap-patterns <INT>`: Number of different random gap masks to index and query. Default: `1`. A window counts as a hit if any pattern matches.
- `--seed <INT>`: Optional RNG seed for reproducible gap masks.

## Canonicalization

All k-mers are canonicalized **before hashing**:
- For each window of length `k`, the forward sequence and its reverse complement are compared.
- The lexicographically smaller of the two is used for hashing.

Important consequence:
- Any non-ACGT(U) base anywhere in the k-mer causes that window to be skipped, even if the ambiguous base is in a gap position.

## Examples (runnable)

Example FASTA files are provided in `examples/`.

### Basic run

```bash
./target/release/Sid examples/ref.fa examples/query.fa -k 31
```

### With gapped k-mers

```bash
./target/release/Sid examples/ref.fa examples/query.fa -k 31 --gaps 5 --seed 123
```

### Multiple gap patterns

```bash
./target/release/Sid examples/ref.fa examples/query.fa -k 31 --gaps 5 --gap-patterns 4 --seed 123
```

### 32-bit hashes (lower memory)

```bash
./target/release/Sid examples/ref.fa examples/query.fa -k 31 --hash32
```

## Output

Two tab-separated lines:

```
query_kmer_hits\t<COUNT> (<PERCENT>% of <TOTAL>)
query_sequences_with_hit\t<COUNT> (<PERCENT>% of <TOTAL>)
```

Counts include comma separators for readability.

## Notes

- FASTA only. Headers are lines starting with `>`.
- Gaps are randomly chosen per pattern, but never in positions 0..7.
- For `k <= 32` without gaps, hashing uses a fast 2-bit rolling encoder.
- For `k > 32` or gapped k-mers with more than 32 kept positions, hashing uses `xxh3` on the compacted sequence.
