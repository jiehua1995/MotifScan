//! 中文：输出模块，负责把统计结果稳定地写成 CSV 文件。
//! English: Output module responsible for serializing scan results into stable CSV files.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};

use crate::motif::Strand;

/// 中文：count 汇总表的固定列顺序。
/// English: Fixed column order for the count-summary CSV.
pub const COUNT_HEADERS: &[&str] = &[
    "motif",
    "sequence",
    "length",
    "reads_with_hit",
    "total_hits",
    "forward_hits",
    "revcomp_hits",
];

/// 中文：read-level hit 明细表的固定列顺序。
/// English: Fixed column order for the read-level hit CSV.
pub const READ_HIT_HEADERS: &[&str] = &[
    "read_id",
    "motif",
    "strand",
    "position",
    "matched_sequence",
];

/// 中文：单个 motif 的汇总统计行。
/// English: One summary row aggregating hit statistics for a single motif.
#[derive(Debug, Clone, Default)]
pub struct CountRow {
    pub motif: String,
    pub sequence: String,
    pub length: usize,
    pub reads_with_hit: u64,
    pub total_hits: u64,
    pub forward_hits: u64,
    pub revcomp_hits: u64,
}

/// 中文：单次命中明细，用于可选的 read-level 报告。
/// English: One hit-detail row used in the optional read-level report.
#[derive(Debug, Clone)]
pub struct ReadHitRow {
    pub read_id: String,
    pub motif: String,
    pub strand: Strand,
    pub position: usize,
    pub matched_sequence: String,
}

/// 中文：统一的 CSV writer 类型别名，隐藏底层写入器细节。
/// English: Shared CSV writer type alias that hides the underlying buffered writer details.
pub type TableWriter = csv::Writer<Box<dyn Write>>;

const OUTPUT_DELIMITER: u8 = b',';

/// 中文：创建一个输出 writer，并确保目标目录存在。
/// English: Creates an output writer and makes sure the destination directory exists first.
pub fn create_writer(path: &Path) -> Result<TableWriter> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("failed to create output directory {}", parent.display()))?;
    }

    let file = File::create(path)
        .with_context(|| format!("failed to create output {}", path.display()))?;
    Ok(csv::WriterBuilder::new()
        .delimiter(OUTPUT_DELIMITER)
        .has_headers(false)
        .from_writer(Box::new(BufWriter::new(file))))
}

    /// 中文：把 motif 汇总结果写入最终的 count CSV 文件。
    /// English: Writes the final motif summary rows into the count CSV file.
pub fn write_count_summary(path: &Path, rows: &[CountRow]) -> Result<()> {
    let mut writer = create_writer(path)?;
    writer.write_record(COUNT_HEADERS)?;
    for row in rows {
        writer.write_record([
            row.motif.as_str(),
            row.sequence.as_str(),
            &row.length.to_string(),
            &row.reads_with_hit.to_string(),
            &row.total_hits.to_string(),
            &row.forward_hits.to_string(),
            &row.revcomp_hits.to_string(),
        ])?;
    }
    writer.flush()?;
    Ok(())
}

/// 中文：向 read-hit 输出文件写入表头。
/// English: Writes the header row for the read-hit output file.
pub fn write_read_hit_headers(writer: &mut TableWriter) -> Result<()> {
    writer.write_record(READ_HIT_HEADERS)?;
    Ok(())
}

/// 中文：批量写入 read-level hit 明细，每个命中一行。
/// English: Writes a batch of read-level hit details, one CSV row per hit.
pub fn write_read_hit_rows(writer: &mut TableWriter, rows: &[ReadHitRow]) -> Result<()> {
    for row in rows {
        writer.write_record([
            row.read_id.as_str(),
            row.motif.as_str(),
            match row.strand {
                Strand::Forward => "forward",
                Strand::Reverse => "revcomp",
            },
            &row.position.to_string(),
            row.matched_sequence.as_str(),
        ])?;
    }
    Ok(())
}
