# MotifScan — Agent Usage Guide (for AI tools)

This document is intended for automated agents and tools that need to understand how to use MotifScan programmatically or from the command line.

Purpose

MotifScan is a streaming, low-memory Rust CLI that counts exact motif occurrences in FASTA/FASTQ reads. It supports multiple motifs (CSV table), optional reverse-complement scanning, and can emit per-read hit CSVs.

Key behaviors

- Input: FASTA, FASTQ, FASTA.GZ, FASTQ.GZ. FASTQ expected in standard 4-line format.
- Matching: exact base matches only (`A`, `C`, `G`, `T`, `U`). Overlapping hits are counted. Reverse complement counts are optional via `--revcomp`.
- Output summary CSV columns: `motif,sequence,length,reads_with_hit,total_hits,forward_hits,revcomp_hits`.
- Read-hit CSV columns: `read_id,motif,strand,position,matched_sequence`.
- Performance: uses SIMD window comparisons and Aho–Corasick automaton for large motif sets. Streams read hits to a writer thread to keep memory usage low.

Primary CLI usage

- Count command (main workflow):

```bash
motifscan count -i <reads> --motifs <motifs.csv> -o <summary.csv>
```

- Single motif from CLI:

```bash
motifscan count -i reads.fastq --motif ATTATGAGAATAGTGTG --motif-name motif1 -o count.csv
```

- Emit read-level hits:

```bash
motifscan count -i reads.fastq --motifs motifs.csv --report-read-hits read_hits.csv -o count.csv
```

Where

- `--motifs <FILE>`: CSV with two columns `name,sequence` (comma-separated, optional header).
- `--motif <SEQUENCE>` and `--motif-name <NAME>`: specify a single motif inline.
- `--revcomp`: include reverse complement matches.
- `-t/--threads <INT>`: number of worker threads.
- `--progress`: show a progress bar on stderr.

Interactions for agents

- Validate inputs before running: ensure motifs file uses only canonical bases; if the agent supplies IUPAC codes, the CLI will reject unless `--iupac` is explicitly supported in future versions.
- Prefer `--motifs` CSV when scanning many motifs; the scanner builds an Aho–Corasick automaton if motif count is large to speed matching.
- To collect per-read hits, set `--report-read-hits` to a writable path; the CLI streams hits so the file may be written incrementally.

Example automation steps for an agent

1. Prepare `motifs.csv` with `name,sequence` rows.
2. Ensure reads file exists (compressed allowed).
3. Run MotifScan with desired options.
4. Monitor exit code: zero indicates success; non-zero indicates error.
5. Read `summary.csv` after completion; if `--report-read-hits` was used, the hits file may be partially available during runtime but fully written on success.

Return codes and error handling

- Exit 0: success.
- Exit non-zero: failure. The CLI writes human-readable errors to stderr. Agents should capture stderr and parse messages for actionable info (e.g., invalid motif file, I/O errors).

Performance and constraints for agents

- Memory: MotifScan is designed to stream data and avoid keeping read-hit rows in memory. Agents can run it on large datasets but should size `--threads` according to CPU availability.
- Cross-platform: Binaries are available as releases; for automation across OSes, download platform-specific release artifacts.

Where to find artifacts and docs

- Releases: check the GitHub Releases page for prebuilt binaries and `.sha256` checksum files.
- Local build: `cargo build --release` produces `target/release/motifscan`.
- Docs: `README.md`, `README_CN.md`, and `doc/release.md` in the repository.

Sample agent pseudo-workflow (shell)

```bash
# prepare
# run motifscan and capture exit and stderr
motifscan count -i reads.fastq --motifs motifs.csv --report-read-hits hits.csv -o summary.csv 2>run.err
RC=$?
if [ $RC -ne 0 ]; then
  cat run.err
  exit $RC
fi
# on success, read summary.csv
cat summary.csv
```

Tips for implementers

- If running inside an orchestrated environment, stream the `--report-read-hits` file to downstream consumers as it is produced.
- When scanning many motifs, expect the Aho–Corasick path to be used; this is more CPU-friendly for large motif sets.
- For deterministic reproducible releases, verify downloaded artifacts with SHA256 before executing.

This document is intended to give automated agents enough context to run MotifScan, check for expected outputs, and handle errors programmatically.