# Changelog

All notable changes to this project will be documented in this file.

## [0.1.3] - 2026-05-06
### Added
- Streaming, low-memory motif scanner CLI with `count` subcommand.
- Support for FASTA/FASTQ and gzipped inputs.
- Aho–Corasick multi-pattern automaton for efficient many-motif scans.
- SIMD-optimized exact window comparisons with runtime feature detection.
- Configurable worker threads and progress reporting.
- Streaming read-hit writer thread to reduce peak memory usage.
- GitHub Actions CI to build artifacts for Linux/macOS/Windows, package with semantic versions, and generate SHA256 checksums.

### Changed
- Replaced unsafe `unwrap()` uses with explicit error handling (`anyhow`).
- Simplified CLI to `count`-only workflow; removed classify/group commands.

### Fixed
- Improve robustness of motif parsing and validation.
- Reduce memory use when producing per-read hits.



