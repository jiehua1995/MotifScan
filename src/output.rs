use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};

use crate::motif::Strand;

pub const COUNT_HEADERS: &[&str] = &[
    "motif",
    "sequence",
    "length",
    "reads_with_hit",
    "total_hits",
    "forward_hits",
    "revcomp_hits",
];

pub const READ_HIT_HEADERS: &[&str] = &[
    "read_id",
    "motif",
    "strand",
    "position",
    "matched_sequence",
];

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

#[derive(Debug, Clone)]
pub struct ReadHitRow {
    pub read_id: String,
    pub motif: String,
    pub strand: Strand,
    pub position: usize,
    pub matched_sequence: String,
}

pub type TableWriter = csv::Writer<Box<dyn Write>>;

const OUTPUT_DELIMITER: u8 = b',';

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

pub fn write_read_hit_headers(writer: &mut TableWriter) -> Result<()> {
    writer.write_record(READ_HIT_HEADERS)?;
    Ok(())
}

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
