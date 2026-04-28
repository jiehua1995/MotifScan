//! 中文：motif 解析与预编译模块，负责读取 CSV、校验序列并准备正向/反向互补模式。
//! English: Motif parsing and compilation module that loads CSV input, validates sequences, and prepares forward/reverse-complement patterns.

use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use anyhow::{anyhow, bail, Context, Result};
use csv::{ReaderBuilder, StringRecord, Trim};

/// 中文：命中链方向枚举，用于区分正向命中还是反向互补命中。
/// English: Strand-direction enum used to distinguish forward hits from reverse-complement hits.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

/// 中文：从命令行或 CSV 读取后的原始 motif，尚未做字节级预处理。
/// English: Raw motif loaded from CLI or CSV before byte-level normalization and compilation.
#[derive(Debug, Clone)]
pub struct RawMotif {
    pub name: String,
    pub sequence: String,
}

/// 中文：单条可匹配模式，目前只保存标准化后的字节序列。
/// English: One compiled match pattern; at the moment it stores only the normalized byte sequence.
#[derive(Debug, Clone)]
pub struct Pattern {
    pub sequence: Vec<u8>,
}

/// 中文：完整编译后的 motif，包含名称、原始展示序列，以及正向/反向互补模式。
/// English: Fully compiled motif containing its name, display sequence, and forward/reverse-complement patterns.
#[derive(Debug, Clone)]
pub struct CompiledMotif {
    pub name: String,
    pub sequence: String,
    pub forward: Pattern,
    pub reverse: Option<Pattern>,
    pub is_palindrome: bool,
}

impl CompiledMotif {
    /// 中文：返回 motif 长度，实际长度来自正向模式的字节数。
    /// English: Returns the motif length, derived from the forward pattern byte count.
    pub fn len(&self) -> usize {
        self.forward.sequence.len()
    }
}

/// 中文：把命令行提供的单个 motif 包装成统一的 `RawMotif` 列表接口。
/// English: Wraps a single CLI-provided motif into the common `RawMotif` list representation.
pub fn load_single_motif(name: &str, sequence: &str) -> Result<Vec<RawMotif>> {
    let sequence = sequence.trim();
    if sequence.is_empty() {
        bail!("motif sequence must not be empty")
    }
    Ok(vec![RawMotif {
        name: name.to_string(),
        sequence: sequence.to_string(),
    }])
}

/// 中文：从两列 CSV 文件加载 motif 列表，并跳过可选表头与注释行。
/// English: Loads motifs from a two-column CSV file, skipping the optional header and comment lines.
pub fn load_motif_file(path: &Path) -> Result<Vec<RawMotif>> {
    let reader = BufReader::new(
        File::open(path)
            .with_context(|| format!("failed to open motif file {}", path.display()))?,
    );
    let mut csv_reader = ReaderBuilder::new()
        .delimiter(b',')
        .has_headers(false)
        .comment(Some(b'#'))
        .trim(Trim::All)
        .from_reader(reader);
    let mut motifs = Vec::new();

    for (record_number, result) in csv_reader.records().enumerate() {
        let record = result.with_context(|| {
            format!(
                "failed to parse motif CSV record {} in {}",
                record_number + 1,
                path.display()
            )
        })?;
        if record.is_empty() {
            continue;
        }

        if is_header(&record) {
            continue;
        }

        if record.len() != 2 {
            bail!(
                "invalid motif CSV format at record {} in {}: expected exactly 2 comma-separated columns",
                record_number + 1,
                path.display()
            )
        }

        let raw = RawMotif {
            name: record.get(0).unwrap().to_string(),
            sequence: record.get(1).unwrap().to_string(),
        };

        motifs.push(raw);
    }

    if motifs.is_empty() {
        bail!("motif file {} did not contain any motifs", path.display())
    }

    Ok(motifs)
}

/// 中文：批量编译所有原始 motif，把它们转成扫描器可以直接使用的结构。
/// English: Compiles a batch of raw motifs into the structures consumed directly by the scanner.
pub fn compile_motifs(
    raw_motifs: &[RawMotif],
    include_revcomp: bool,
) -> Result<Vec<CompiledMotif>> {
    raw_motifs
        .iter()
        .map(|raw| compile_motif(raw, include_revcomp))
        .collect()
}

    /// 中文：编译单个 motif，完成标准化、合法性检查和反向互补预计算。
    /// English: Compiles a single motif by normalizing it, validating allowed bases, and precomputing the reverse complement when requested.
pub fn compile_motif(raw: &RawMotif, include_revcomp: bool) -> Result<CompiledMotif> {
    let sequence = normalize_sequence(&raw.sequence);
    if sequence.is_empty() {
        bail!("motif '{}' has an empty sequence", raw.name)
    }

    validate_motif_sequence(&raw.name, &sequence)?;

    let reverse_sequence = reverse_complement(&sequence)?;
    let is_palindrome = sequence == reverse_sequence;

    let forward = Pattern {
        sequence: sequence.clone(),
    };
    let reverse = if include_revcomp && !is_palindrome {
        Some(Pattern {
            sequence: reverse_sequence,
        })
    } else {
        None
    };

    Ok(CompiledMotif {
        name: raw.name.clone(),
        sequence: String::from_utf8(sequence.clone()).unwrap(),
        forward,
        reverse,
        is_palindrome,
    })
}

/// 中文：把序列裁掉首尾空白并统一转成大写，便于后续按字节比较。
/// English: Trims surrounding whitespace and converts a sequence to uppercase for byte-wise matching.
pub fn normalize_sequence(sequence: &str) -> Vec<u8> {
    sequence
        .trim()
        .as_bytes()
        .iter()
        .map(|byte| byte.to_ascii_uppercase())
        .collect()
}

// 中文：检查 motif 是否只包含当前 exact 模式允许的碱基字符。
// English: Verifies that a motif contains only the bases accepted by the current exact-match implementation.
fn validate_motif_sequence(name: &str, sequence: &[u8]) -> Result<()> {
    for base in sequence {
        let upper = base.to_ascii_uppercase();
        if matches!(upper, b'A' | b'C' | b'G' | b'T' | b'U') {
            continue;
        }

        bail!(
            "motif '{}' contains unsupported base '{}' ; exact matching only supports A/C/G/T/U motifs",
            name,
            upper as char
        )
    }

    Ok(())
}

/// 中文：计算给定序列的反向互补序列。
/// English: Computes the reverse-complement sequence for the provided bases.
pub fn reverse_complement(sequence: &[u8]) -> Result<Vec<u8>> {
    sequence
        .iter()
        .rev()
        .map(|base| complement(*base))
        .collect()
}

// 中文：返回单个碱基的互补碱基，用于构建反向互补 motif。
// English: Returns the complementary base for a single nucleotide when building reverse-complement motifs.
fn complement(base: u8) -> Result<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Ok(b'T'),
        b'C' => Ok(b'G'),
        b'G' => Ok(b'C'),
        b'T' | b'U' => Ok(b'A'),
        other => Err(anyhow!("unsupported motif character '{}'", other as char)),
    }
}

// 中文：判断一条 CSV 记录是否是 `name,sequence` 这样的表头。
// English: Detects whether a CSV record is the optional `name,sequence` header row.
fn is_header(record: &StringRecord) -> bool {
    matches!(
        (record.get(0), record.get(1), record.len()),
        (Some(name), Some(sequence), 2)
            if name.eq_ignore_ascii_case("name")
                && sequence.eq_ignore_ascii_case("sequence")
    )
}

#[cfg(test)]
mod tests {
    use super::{compile_motif, is_header, load_motif_file, reverse_complement, RawMotif};

    use csv::StringRecord;
    use tempfile::NamedTempFile;
    use std::io::Write;

    #[test]
    // 中文：验证标准碱基的反向互补实现是否正确。
    // English: Verifies that reverse-complement generation works correctly for canonical bases.
    fn reverse_complement_handles_canonical_bases() {
        let rc = reverse_complement(b"ATGCT").unwrap();
        assert_eq!(String::from_utf8(rc).unwrap(), "AGCAT");
    }

    #[test]
    // 中文：验证回文 motif 不会重复生成反向互补模式。
    // English: Verifies that palindromic motifs do not generate a redundant reverse-complement pattern.
    fn compiles_palindromic_reverse_once() {
        let motif = RawMotif {
            name: "pal".to_string(),
            sequence: "ATAT".to_string(),
        };
        let compiled = compile_motif(&motif, true).unwrap();
        assert!(compiled.is_palindrome);
        assert!(compiled.reverse.is_none());
    }

    #[test]
    // 中文：验证包含简并碱基的 motif 会在 exact-only 模式下被拒绝。
    // English: Verifies that motifs containing ambiguous bases are rejected in exact-only mode.
    fn rejects_non_exact_motif_characters() {
        let motif = RawMotif {
            name: "iupac".to_string(),
            sequence: "ATGRN".to_string(),
        };
        let error = compile_motif(&motif, false).unwrap_err();
        assert!(error
            .to_string()
            .contains("exact matching only supports A/C/G/T/U motifs"));
    }

    #[test]
    // 中文：验证 CSV 表头识别逻辑只匹配 `name,sequence`，不会误判普通数据行。
    // English: Verifies that header detection matches only `name,sequence` and does not misclassify ordinary data rows.
    fn detects_optional_csv_header() {
        let header = StringRecord::from(vec!["name", "sequence"]);
        let row = StringRecord::from(vec!["motif1", "ACTG"]);
        assert!(is_header(&header));
        assert!(!is_header(&row));
    }

    #[test]
    // 中文：验证 motif CSV 文件加载路径能正确跳过表头并读取数据行。
    // English: Verifies that motif CSV loading skips the header and reads the data rows correctly.
    fn loads_comma_separated_motif_file() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "name,sequence").unwrap();
        writeln!(file, "motif1,ATTATGAGAATAGTGTG").unwrap();
        writeln!(file, "motif2,TTCATTCATGGTGGCAGTAAAATGTTTATTGTG").unwrap();

        let motifs = load_motif_file(file.path()).unwrap();
        assert_eq!(motifs.len(), 2);
        assert_eq!(motifs[0].name, "motif1");
        assert_eq!(motifs[1].sequence, "TTCATTCATGGTGGCAGTAAAATGTTTATTGTG");
    }
}
