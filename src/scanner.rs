//! 中文：扫描引擎模块，负责并行读取 records、执行 exact motif 匹配，并汇总输出结果。
//! English: Scan-engine module responsible for parallel record processing, exact motif matching, and result aggregation.

use aho_corasick::AhoCorasick;
use std::sync::mpsc::{channel, Sender};
use std::sync::Arc;
use std::thread;
use std::time::Duration;

use anyhow::Result;
use indicatif::{ProgressBar, ProgressStyle};
use memchr::memchr_iter;
use rayon::prelude::*;

#[cfg(target_arch = "x86")]
use std::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

use crate::cli::CountArgs;
use crate::io::{open_record_reader, ProgressSnapshot, Record, RecordReader};
use crate::motif::{compile_motifs, load_motif_file, load_single_motif, CompiledMotif, Strand};
use crate::output::{
    create_writer, write_count_summary, write_read_hit_headers, write_read_hit_rows, CountRow,
    ReadHitRow,
};

const DEFAULT_CHUNK_SIZE: usize = 512;

type HitWriterResult = (
    Option<Sender<Vec<ReadHitRow>>>,
    Option<thread::JoinHandle<()>>,
);

/// 简单封装，保存 Aho-Corasick 自动机与 pattern id -> (motif_index, strand, length) 映射
struct AhoIndex {
    ac: Arc<AhoCorasick>,
    map: Arc<Vec<(usize, Strand, usize)>>,
}

/// 中文：单条 record 扫描后的局部结果，稍后会被合并进全局计数。
/// English: Per-record scan result that is merged later into the global counters.
#[derive(Debug, Clone)]
struct RecordResult {
    motif_hits: Vec<MotifHitSummary>,
    read_hits: Vec<ReadHitRow>,
}

/// 中文：单个 motif 在单条 record 上的统计摘要。
/// English: Summary of how one motif behaved on one record.
#[derive(Debug, Clone)]
struct MotifHitSummary {
    motif_index: usize,
    total_hits: u64,
    forward_hits: u64,
    revcomp_hits: u64,
    read_has_hit: bool,
}

/// 中文：某条模式扫描的原始结果，包括命中次数以及可选的位置列表。
/// English: Raw result for scanning one pattern, including hit count and an optional list of positions.
#[derive(Debug, Clone, Default)]
struct PatternScanResult {
    hit_count: u64,
    positions: Vec<usize>,
}

/// 中文：`count` 子命令的主执行入口。
/// English: Main execution entry for the `count` subcommand.
///
/// 中文：这个函数负责验证参数、加载 motif、打开输入流、驱动扫描循环，并把最终汇总结果写成 CSV。
/// English: This function validates arguments, loads motifs, opens the input stream, drives the scan loop, and writes the final aggregated CSV output.
pub fn run_count(args: &CountArgs) -> Result<()> {
    args.validate()?;
    let raw_motifs = if let Some(sequence) = &args.motif {
        load_single_motif(&args.motif_name, sequence)?
    } else {
        load_motif_file(args.motifs.as_ref().unwrap())?
    };
    let motifs = compile_motifs(&raw_motifs, args.revcomp)?;

    let mut reader = open_record_reader(&args.input)?;
    let mut rows = initialize_rows(&motifs);
    let (hit_sender, writer_handle) = maybe_spawn_hit_writer(args.report_read_hits.as_deref())?;
    let mut progress =
        ScanProgress::new(&reader, args.progress, "count", &args.input, motifs.len())?;

    // Build Aho-Corasick automaton when many motifs exist to speed up multi-pattern scans
    let aho_index = if motifs.len() >= 8 {
        // collect patterns and mapping
        let mut patterns_owned: Vec<String> = Vec::new();
        let mut mapping: Vec<(usize, Strand, usize)> = Vec::new();
        for (i, m) in motifs.iter().enumerate() {
            let f = String::from_utf8_lossy(&m.forward.sequence).into_owned();
            mapping.push((i, Strand::Forward, f.len()));
            patterns_owned.push(f);
            if let Some(r) = &m.reverse {
                let rs = String::from_utf8_lossy(&r.sequence).into_owned();
                mapping.push((i, Strand::Reverse, rs.len()));
                patterns_owned.push(rs);
            }
        }
        let refs: Vec<&str> = patterns_owned.iter().map(|s| s.as_str()).collect();
        let ac = AhoCorasick::new(&refs);
        Some(AhoIndex {
            ac: Arc::new(ac),
            map: Arc::new(mapping),
        })
    } else {
        None
    };

    scan_records(
        &mut reader,
        &motifs,
        &mut progress,
        hit_sender.as_ref(),
        aho_index.as_ref(),
        &mut rows,
    )?;

    // close the sender so the writer thread sees EOF
    drop(hit_sender);
    if let Some(handle) = writer_handle {
        // wait for writer thread to finish flushing
        handle.join().expect("read-hit writer thread panicked");
    }
    progress.finish();
    write_count_summary(&args.output, &rows)
}

// 中文：按需创建 read-hit 明细输出器；如果用户没请求，就返回 `None`。
// English: Creates the optional read-hit writer when requested; otherwise returns `None`.
fn maybe_spawn_hit_writer(path: Option<&std::path::Path>) -> Result<HitWriterResult> {
    match path {
        Some(path) => {
            let path = path.to_path_buf();
            let (tx, rx) = channel::<Vec<ReadHitRow>>();

            let handle = thread::spawn(move || {
                // Writer thread: open writer and consume batches
                match create_writer(&path) {
                    Ok(mut writer) => {
                        if let Err(e) = write_read_hit_headers(&mut writer) {
                            eprintln!("failed to write read-hit headers: {}", e);
                            return;
                        }

                        while let Ok(batch) = rx.recv() {
                            if let Err(e) = write_read_hit_rows(&mut writer, &batch) {
                                eprintln!("failed to write read-hit rows: {}", e);
                                break;
                            }
                            if let Err(e) = writer.flush() {
                                eprintln!("failed to flush read-hit writer: {}", e);
                                break;
                            }
                        }
                    }
                    Err(e) => {
                        eprintln!("failed to create read-hit writer {}: {}", path.display(), e);
                    }
                }
            });

            Ok((Some(tx), Some(handle)))
        }
        None => Ok((None, None)),
    }
}

// 中文：根据 motif 列表预先创建汇总行，后续扫描时只需原地累加。
// English: Pre-allocates summary rows from the motif list so later scan passes can update them in place.
fn initialize_rows(motifs: &[CompiledMotif]) -> Vec<CountRow> {
    motifs
        .iter()
        .map(|motif| CountRow {
            motif: motif.name.clone(),
            sequence: motif.sequence.clone(),
            length: motif.len(),
            reads_with_hit: 0,
            total_hits: 0,
            forward_hits: 0,
            revcomp_hits: 0,
        })
        .collect()
}

// 中文：以 chunk 为单位推进整个扫描过程；每个 chunk 内部并行处理，每个 chunk 结束后统一归并结果。
// English: Advances the full scan in chunk units; each chunk is processed in parallel and merged only after the chunk finishes.
fn scan_records(
    reader: &mut RecordReader,
    motifs: &[CompiledMotif],
    progress: &mut ScanProgress,
    hit_sender: Option<&Sender<Vec<ReadHitRow>>>,
    aho_index: Option<&AhoIndex>,
    rows: &mut [CountRow],
) -> Result<()> {
    loop {
        let chunk = reader.next_chunk(DEFAULT_CHUNK_SIZE)?;
        if chunk.is_empty() {
            break;
        }
        let chunk_reads = chunk.len() as u64;
        let chunk_bases = chunk.iter().map(|record| record.seq.len() as u64).sum();

        let emit_read_hits = hit_sender.is_some();
        let record_results: Vec<RecordResult> = if let Some(aho) = aho_index {
            chunk
                .into_par_iter()
                .map(|record| scan_record_aho(&record, motifs, emit_read_hits, aho))
                .collect()
        } else {
            chunk
                .into_par_iter()
                .map(|record| scan_record(&record, motifs, emit_read_hits))
                .collect()
        };

        for record_result in record_results {
            merge_record_result(&record_result, rows);
            if emit_read_hits && !record_result.read_hits.is_empty() {
                if let Some(sender) = hit_sender {
                    // best-effort: ignore send errors if receiver is closed
                    let _ = sender.send(record_result.read_hits);
                }
            }
        }

        progress.update(chunk_reads, chunk_bases, reader.progress_snapshot());
    }

    Ok(())
}

/// 中文：进度展示状态；关闭时用 `Disabled` 避免热路径里到处判断具体组件。
/// English: Progress-report state; the `Disabled` variant keeps the hot path free from UI-specific branching details.
enum ScanProgress {
    Enabled(ProgressState),
    Disabled,
}

/// 中文：进度条运行时状态，保存累计 reads、碱基数和界面展示对象。
/// English: Runtime progress-bar state storing cumulative reads, bases, and the UI handle itself.
struct ProgressState {
    bar: ProgressBar,
    reads_processed: u64,
    bases_processed: u64,
    input_name: String,
    mode: &'static str,
    motif_count: usize,
}

impl ScanProgress {
    /// 中文：根据用户是否开启 `--progress` 来创建真实进度条或空实现。
    /// English: Creates either a real progress bar or a disabled no-op state depending on `--progress`.
    fn new(
        reader: &RecordReader,
        enabled: bool,
        mode: &'static str,
        input_path: &std::path::Path,
        motif_count: usize,
    ) -> Result<Self> {
        if !enabled {
            return Ok(Self::Disabled);
        }

        let snapshot = reader.progress_snapshot();
        let bar = ProgressBar::new(snapshot.total_bytes);
        let style = ProgressStyle::with_template(
            "{spinner:.green} {msg}\n[{bar:40.cyan/blue}] {percent:>3}% | {bytes}/{total_bytes} | eta {eta_precise}",
        )?;
        bar.set_style(style.progress_chars("=>-"));
        bar.enable_steady_tick(Duration::from_millis(120));
        bar.set_position(snapshot.bytes_read);
        let input_name = input_path
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("input")
            .to_string();
        let mut state = ProgressState {
            bar,
            reads_processed: 0,
            bases_processed: 0,
            input_name,
            mode,
            motif_count,
        };
        state.refresh_message();
        Ok(Self::Enabled(state))
    }

    /// 中文：在处理完一个 chunk 后刷新累计进度和界面文本。
    /// English: Updates cumulative progress and refreshes the progress-bar message after one chunk completes.
    fn update(&mut self, chunk_reads: u64, chunk_bases: u64, snapshot: ProgressSnapshot) {
        if let Self::Enabled(state) = self {
            state.reads_processed += chunk_reads;
            state.bases_processed += chunk_bases;
            state
                .bar
                .set_position(snapshot.bytes_read.min(snapshot.total_bytes));
            state.refresh_message();
        }
    }

    /// 中文：扫描结束时收尾进度条显示。
    /// English: Finalizes the progress display when scanning is complete.
    fn finish(&self) {
        if let Self::Enabled(state) = self {
            state.bar.finish_and_clear();
        }
    }
}

impl ProgressState {
    /// 中文：重新计算并更新进度条消息文本，例如 reads/s 和平均 read 长度。
    /// English: Recomputes and updates the progress-bar message, such as reads per second and average read length.
    fn refresh_message(&mut self) {
        let elapsed = self.bar.elapsed().as_secs_f64();
        let reads_per_sec = if elapsed > 0.0 {
            self.reads_processed as f64 / elapsed
        } else {
            0.0
        };
        let avg_read_len = if self.reads_processed > 0 {
            self.bases_processed as f64 / self.reads_processed as f64
        } else {
            0.0
        };
        let elapsed_label = format_duration(self.bar.elapsed());
        self.bar.set_message(format!(
            "{} {} | motifs {} | reads {} | avg_len {:.1} bp | {:.1} reads/s | elapsed {}",
            self.mode,
            self.input_name,
            self.motif_count,
            self.reads_processed,
            avg_read_len,
            reads_per_sec,
            elapsed_label,
        ));
    }
}

// 中文：把持续时间格式化成紧凑的 `HH:MM:SS` 或 `MM:SS`，用于进度消息展示。
// English: Formats a duration into compact `HH:MM:SS` or `MM:SS` text for the progress message.
fn format_duration(duration: Duration) -> String {
    let total_seconds = duration.as_secs();
    let hours = total_seconds / 3600;
    let minutes = (total_seconds % 3600) / 60;
    let seconds = total_seconds % 60;

    if hours > 0 {
        format!("{hours:02}:{minutes:02}:{seconds:02}")
    } else {
        format!("{minutes:02}:{seconds:02}")
    }
}

// 中文：把单条 record 的局部统计加到全局汇总表里。
// English: Merges one record's local statistics into the global summary table.
fn merge_record_result(record_result: &RecordResult, rows: &mut [CountRow]) {
    for motif_hit in &record_result.motif_hits {
        let row = &mut rows[motif_hit.motif_index];
        row.total_hits += motif_hit.total_hits;
        row.forward_hits += motif_hit.forward_hits;
        row.revcomp_hits += motif_hit.revcomp_hits;
        if motif_hit.read_has_hit {
            row.reads_with_hit += 1;
        }
    }
}

// 中文：扫描单条 record 上的所有 motif，并按需要收集 read-level hit 明细。
// English: Scans all motifs on one record and optionally collects read-level hit details.
fn scan_record(record: &Record, motifs: &[CompiledMotif], emit_read_hits: bool) -> RecordResult {
    let mut motif_hits = Vec::with_capacity(motifs.len());
    let mut read_hits = Vec::new();

    for (motif_index, motif) in motifs.iter().enumerate() {
        let forward_scan = scan_pattern(record, &motif.forward, emit_read_hits);
        if forward_scan.hit_count > 0 && emit_read_hits {
            append_read_hits(
                &mut read_hits,
                record,
                motif,
                Strand::Forward,
                &forward_scan.positions,
                motif.len(),
            );
        }

        let reverse_scan = motif
            .reverse
            .as_ref()
            .map(|pattern| scan_pattern(record, pattern, emit_read_hits))
            .unwrap_or_default();
        if reverse_scan.hit_count > 0 && emit_read_hits {
            append_read_hits(
                &mut read_hits,
                record,
                motif,
                Strand::Reverse,
                &reverse_scan.positions,
                motif.len(),
            );
        }

        let forward_hits = forward_scan.hit_count;
        let revcomp_hits = reverse_scan.hit_count;
        let total_hits = forward_hits + revcomp_hits;
        let read_has_hit = total_hits > 0;

        motif_hits.push(MotifHitSummary {
            motif_index,
            total_hits,
            forward_hits,
            revcomp_hits,
            read_has_hit,
        });
    }

    RecordResult {
        motif_hits,
        read_hits,
    }
}

// 中文：把模式命中的位置列表展开成真正的 read-hit 输出行。
// English: Expands a list of hit positions into concrete read-hit output rows.
fn append_read_hits(
    sink: &mut Vec<ReadHitRow>,
    record: &Record,
    motif: &CompiledMotif,
    strand: Strand,
    positions: &[usize],
    motif_length: usize,
) {
    for position in positions {
        let window = &record.seq[*position..(*position + motif_length)];
        sink.push(ReadHitRow {
            read_id: record.id.clone(),
            motif: motif.name.clone(),
            strand,
            position: *position,
            matched_sequence: String::from_utf8_lossy(window).into_owned(),
        });
    }
}

// 使用 Aho-Corasick 自动机扫描单条 record（多 pattern 路径）
fn scan_record_aho(
    record: &Record,
    motifs: &[CompiledMotif],
    collect_positions: bool,
    aho: &AhoIndex,
) -> RecordResult {
    // prepare per-motif accumulators
    let mut motif_acc: Vec<MotifHitSummary> = motifs
        .iter()
        .enumerate()
        .map(|(i, _)| MotifHitSummary {
            motif_index: i,
            total_hits: 0,
            forward_hits: 0,
            revcomp_hits: 0,
            read_has_hit: false,
        })
        .collect();

    let mut read_hits = Vec::new();
    let seq_str = String::from_utf8_lossy(&record.seq).into_owned();
    for mat in aho.ac.find_iter(&seq_str) {
        let pid = mat.pattern();
        if let Some((motif_idx, strand, plen)) = aho.map.get(pid).cloned() {
            let start = mat.start();
            motif_acc[motif_idx].total_hits += 1;
            match strand {
                Strand::Forward => motif_acc[motif_idx].forward_hits += 1,
                Strand::Reverse => motif_acc[motif_idx].revcomp_hits += 1,
            }
            motif_acc[motif_idx].read_has_hit = true;
            if collect_positions {
                let window = &record.seq[start..start + plen];
                read_hits.push(ReadHitRow {
                    read_id: record.id.clone(),
                    motif: motifs[motif_idx].name.clone(),
                    strand,
                    position: start,
                    matched_sequence: String::from_utf8_lossy(window).into_owned(),
                });
            }
        }
    }

    RecordResult {
        motif_hits: motif_acc,
        read_hits,
    }
}

// 中文：扫描一条具体 pattern，在 exact 模式下统计命中次数，并在需要时记录所有位置。
// English: Scans one concrete pattern in exact mode, counting hits and recording positions when requested.
fn scan_pattern(
    record: &Record,
    pattern: &crate::motif::Pattern,
    collect_positions: bool,
) -> PatternScanResult {
    if pattern.sequence.len() > record.seq.len() {
        return PatternScanResult::default();
    }

    let mut result = PatternScanResult {
        hit_count: 0,
        positions: if collect_positions {
            Vec::with_capacity(4)
        } else {
            Vec::new()
        },
    };

    for position in exact_positions_iter(&record.seq, &pattern.sequence) {
        result.hit_count += 1;
        if collect_positions {
            result.positions.push(position);
        }
    }

    result
}

#[cfg(test)]
// 中文：测试辅助函数，直接返回所有 exact 命中位置，方便断言。
// English: Test helper that materializes all exact-match positions for straightforward assertions.
fn exact_positions(sequence: &[u8], pattern: &[u8]) -> Vec<usize> {
    exact_positions_iter(sequence, pattern).collect()
}

// 中文：exact matching 的核心候选迭代器：先用 `memchr` 找首字节，再用次字节、末字节和整窗比较做快速剪枝。
// English: Core exact-match iterator: it uses `memchr` for the first byte, then prunes with second-byte, last-byte, and full-window checks.
fn exact_positions_iter<'a>(
    sequence: &'a [u8],
    pattern: &'a [u8],
) -> impl Iterator<Item = usize> + 'a {
    let pattern_len = pattern.len();
    let second_base = pattern.get(1).copied();
    let last_base = pattern.last().copied().unwrap_or(pattern[0]);

    memchr_iter(pattern[0], sequence).filter(move |&position| {
        if position + pattern_len > sequence.len() {
            return false;
        }
        if let Some(second_base) = second_base {
            if sequence[position + 1] != second_base {
                return false;
            }
        }
        if sequence[position + pattern_len - 1] != last_base {
            return false;
        }
        exact_match_window(&sequence[position..position + pattern_len], pattern)
    })
}

#[inline]
// 中文：比较一个窗口和 motif 是否完全相等；在 x86/x86_64 上会优先使用 SIMD 快路径。
// English: Compares one candidate window with the motif for exact equality, preferring SIMD fast paths on x86/x86_64.
fn exact_match_window(window: &[u8], pattern: &[u8]) -> bool {
    if window.len() != pattern.len() {
        return false;
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if window.len() >= 32 && std::is_x86_feature_detected!("avx2") {
            unsafe {
                return avx2_equal(window, pattern);
            }
        }
        if window.len() >= 16 && std::is_x86_feature_detected!("sse2") {
            unsafe {
                return sse2_equal(window, pattern);
            }
        }
    }

    window == pattern
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
// 中文：使用 AVX2 按 32 字节块比较两个切片；只在 CPU 支持且窗口足够长时调用。
// English: Uses AVX2 to compare two slices in 32-byte chunks; called only when the CPU supports it and the window is long enough.
unsafe fn avx2_equal(window: &[u8], pattern: &[u8]) -> bool {
    let mut offset = 0;
    while offset + 32 <= window.len() {
        let lhs = _mm256_loadu_si256(window.as_ptr().add(offset) as *const __m256i);
        let rhs = _mm256_loadu_si256(pattern.as_ptr().add(offset) as *const __m256i);
        let cmp = _mm256_cmpeq_epi8(lhs, rhs);
        if _mm256_movemask_epi8(cmp) != -1 {
            return false;
        }
        offset += 32;
    }
    window[offset..] == pattern[offset..]
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
// 中文：使用 SSE2 按 16 字节块比较两个切片；作为 AVX2 不可用时的次优 SIMD 路径。
// English: Uses SSE2 to compare two slices in 16-byte chunks as the fallback SIMD path when AVX2 is unavailable.
unsafe fn sse2_equal(window: &[u8], pattern: &[u8]) -> bool {
    let mut offset = 0;
    while offset + 16 <= window.len() {
        let lhs = _mm_loadu_si128(window.as_ptr().add(offset) as *const __m128i);
        let rhs = _mm_loadu_si128(pattern.as_ptr().add(offset) as *const __m128i);
        let cmp = _mm_cmpeq_epi8(lhs, rhs);
        if _mm_movemask_epi8(cmp) != 0xFFFF {
            return false;
        }
        offset += 16;
    }
    window[offset..] == pattern[offset..]
}

#[cfg(test)]
mod tests {
    use crate::io::{Record, SourceFormat};
    use crate::motif::{compile_motifs, RawMotif};

    use super::{exact_match_window, exact_positions, scan_record};

    // 中文：构造一个最小 record，方便在单元测试里直接驱动扫描逻辑。
    // English: Builds a minimal record object so unit tests can drive the scan logic directly.
    fn demo_record(id: &str, seq: &str, qual: Option<Vec<u8>>) -> Record {
        Record {
            id: id.to_string(),
            seq: seq.as_bytes().to_vec(),
            qual,
            source_format: SourceFormat::Fastq,
        }
    }

    #[test]
    // 中文：验证 exact 匹配会保留重叠命中，而不是只返回不重叠窗口。
    // English: Verifies that exact matching keeps overlapping hits instead of only returning disjoint windows.
    fn exact_matching_finds_overlapping_hits() {
        let hits = exact_positions(b"AAAAA", b"AAA");
        assert_eq!(hits, vec![0, 1, 2]);
    }

    #[test]
    // 中文：验证 `reads_with_hit` 和 `total_hits` 的统计语义不同：一条 read 多次命中时只算一次 read，但会累计多个 hit。
    // English: Verifies that `reads_with_hit` and `total_hits` have different semantics: one read can contribute once to the former and multiple times to the latter.
    fn reads_with_hit_and_total_hits_are_distinct() {
        let motifs = compile_motifs(
            &[RawMotif {
                name: "m1".to_string(),
                sequence: "AAA".to_string(),
            }],
            false,
        )
        .unwrap();
        let result = scan_record(
            &demo_record("r1", "AAAAA", Some(vec![40; 5])),
            &motifs,
            false,
        );
        assert_eq!(result.motif_hits[0].total_hits, 3);
        assert!(result.motif_hits[0].read_has_hit);
    }

    #[test]
    // 中文：验证开启反向互补后，scanner 能识别来自 reverse-complement 链的命中。
    // English: Verifies that reverse-complement hits are detected correctly when revcomp scanning is enabled.
    fn reverse_complement_hits_are_detected() {
        let motifs = compile_motifs(
            &[RawMotif {
                name: "m1".to_string(),
                sequence: "ATTATGAGAATAGTGTG".to_string(),
            }],
            true,
        )
        .unwrap();
        let reverse = "CACACTATTCTCATAAT";
        let result = scan_record(
            &demo_record("r1", reverse, Some(vec![40; reverse.len()])),
            &motifs,
            true,
        );
        assert_eq!(result.motif_hits[0].revcomp_hits, 1);
    }

    #[test]
    // 中文：验证 exact 模式不会把包含 `N` 的 read 片段误当成普通精确匹配。
    // English: Verifies that exact mode does not mistakenly treat read windows containing `N` as exact matches.
    fn exact_mode_does_not_match_ambiguous_motif_literals() {
        let motifs = compile_motifs(
            &[RawMotif {
                name: "m1".to_string(),
                sequence: "ATGAA".to_string(),
            }],
            false,
        )
        .unwrap();
        let result = scan_record(
            &demo_record("r1", "ATGNN", Some(vec![40; 5])),
            &motifs,
            false,
        );
        assert_eq!(result.motif_hits[0].total_hits, 0);
    }

    #[test]
    // 中文：验证 SIMD 快路径和标量回退路径在窗口比较上的结果一致。
    // English: Verifies that the SIMD fast path and the scalar fallback agree on window-comparison results.
    fn simd_window_match_falls_back_correctly() {
        let pattern = b"ATTATGAGAATAGTGTGATTATGAGAATAGTGTG";
        assert!(exact_match_window(pattern, pattern));
        assert!(!exact_match_window(
            pattern,
            b"ATTATGAGAATAGTGTGATTATGAGAATAGTGTA"
        ));
    }

    #[test]
    // Integration-ish test: run run_count on test/ sample data and ensure outputs are produced
    fn run_count_writes_outputs() {
        use crate::cli::CountArgs;
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write;
        use tempfile::tempdir;

        let tmp = tempdir().unwrap();
        let out_count = tmp.path().join("count.csv");
        let out_hits = tmp.path().join("read_hits.csv");

        // create a small FASTQ and gzip it
        let fq_path = tmp.path().join("reads_sample.fastq.gz");
        let fq_contents = b"@r1\nATTATGAGAATAGTGTG\n+\nFFFFFFFFFFFFFFFFF\n@r2\nATGAA\n+\nFFFFF\n";
        {
            let fq_file = std::fs::File::create(&fq_path).unwrap();
            let mut enc = GzEncoder::new(fq_file, Compression::default());
            enc.write_all(fq_contents).unwrap();
            enc.finish().unwrap();
        }

        // create a motif CSV
        let motifs_path = tmp.path().join("motifs.csv");
        std::fs::write(
            &motifs_path,
            "name,sequence\nmotif1,ATTATGAGAATAGTGTG\nmotif2,ATGAA\n",
        )
        .unwrap();

        let args = CountArgs {
            input: fq_path,
            motif: None,
            motif_name: "motif".to_string(),
            motifs: Some(motifs_path),
            revcomp: true,
            threads: 1,
            progress: false,
            output: out_count.clone(),
            report_read_hits: Some(out_hits.clone()),
        };

        // run
        super::run_count(&args).expect("run_count failed");

        // check files exist and contain expected header
        let count_txt = std::fs::read_to_string(out_count).unwrap();
        assert!(count_txt.contains("motif,sequence,length"));
        let hits_txt = std::fs::read_to_string(out_hits).unwrap();
        assert!(hits_txt.contains("read_id,motif,strand,position,matched_sequence"));
    }
}
