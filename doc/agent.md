# MotifScan Agent Guide

This document is for automated agents and tooling that need to run MotifScan from the command line.

## Purpose

MotifScan is a streaming, low-memory Rust CLI that counts exact motif occurrences in FASTA or FASTQ reads. It supports:

- Single motifs provided inline with `--motif`
- Multiple motifs loaded from a two-column CSV file with `--motifs`
- Optional reverse-complement scanning with `--revcomp`
- Optional read-level hit output with `--report-read-hits`
- Optional progress display with `--progress`
- Info/debug logging with `--verbose` and `--debug`

## Build

- Required dependency: Rust toolchain
- Local build command: `cargo build --release`
- Binary path: `target/release/motifscan`

## Input Formats

- FASTA
- FASTQ
- Gzipped FASTA and FASTQ files with `.gz` extension
- FASTQ is expected to use the standard four-line layout

## Matching Rules

- Exact base matching only
- Supported motif bases: `A`, `C`, `G`, `T`, `U`
- Overlapping hits are counted
- Reverse-complement hits are counted separately when enabled
- Palindromic motifs are not double-counted in reverse-complement mode

## CLI Overview

Main workflow:

```bash
motifscan count -i <reads> --motifs <motifs.csv> -o <summary.csv>
```

Single motif:

```bash
motifscan count -i reads.fastq --motif ATTATGAGAATAGTGTG --motif-name motif1 -o count.csv
```

Read-level hits:

```bash
motifscan count -i reads.fastq --motifs motifs.csv --report-read-hits read_hits.csv -o count.csv
```

Useful flags:

- `-i`, `--input <FILE>`: reads file
- `--motif <SEQUENCE>`: one motif provided inline
- `--motif-name <NAME>`: name used for `--motif`, default `motif`
- `--motifs <FILE>`: two-column CSV motif table
- `--revcomp`: scan reverse complements too
- `-t`, `--threads <INT>`: worker threads
- `--progress`: show a progress bar on stderr
- `--verbose`: enable info-level logs
- `--debug`: enable debug-level logs
- `-o`, `--output <FILE>`: summary CSV output
- `--report-read-hits <FILE>`: optional read-level hit CSV output

## Output Formats

Summary CSV columns:

```text
motif,sequence,length,reads_with_hit,total_hits,forward_hits,revcomp_hits
```

Read-hit CSV columns:

```text
read_id,motif,strand,position,matched_sequence
```

## Error Handling

- Exit code `0` indicates success.
- Non-zero exit codes indicate failure.
- Errors are printed to stderr.
- If the read-hit writer fails, the process exits non-zero instead of silently succeeding.

## Performance Notes

- The scanner streams input and does not buffer the entire dataset in memory.
- For larger motif sets, it enables an Aho-Corasick path automatically.
- Read-level hits are written through a dedicated writer thread.

## Benchmark Script

The repository includes [`doc/benchmark.sh`](./benchmark.sh), which:

- Generates random motifs and synthetic FASTA/FASTQ data under `benchmark_work/`
- Builds the release binary
- Packages and unpacks the binary archive
- Runs MotifScan against the generated data
- Verifies output against expected CSV files
- Measures wall-clock time for each phase

Environment variables supported by the script:

- `BENCH_WORKDIR`: benchmark work directory, must stay under the repository root
- `BENCH_READS`: number of synthetic reads
- `BENCH_MOTIFS`: number of synthetic motifs
- `BENCH_READ_LEN`: read length
- `BENCH_SEED`: random seed
- `BENCH_INSERT_POS`: insertion position for synthetic motif hits

Example:

```bash
BENCH_READS=20000 BENCH_MOTIFS=128 BENCH_READ_LEN=160 doc/benchmark.sh
```

## Notes for Agents

- Use `--motifs` for larger scans and `--motif` for single-motif runs.
- Use `--report-read-hits` only when read-level details are needed; it increases output size.
- Prefer `--verbose` for human-readable run summaries, `--debug` for per-chunk diagnostics.
- If you need to compare output, use exact CSV header and column ordering shown above.
- Do not rely on `test/` or any README fixture references when automating, because those files are not part of the uploaded runtime context.
