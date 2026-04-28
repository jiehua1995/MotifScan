use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{anyhow, bail, Context, Result};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

#[derive(Debug, Clone)]
pub struct RawMotif {
    pub name: String,
    pub sequence: String,
}

#[derive(Debug, Clone)]
pub struct Pattern {
    pub sequence: Vec<u8>,
    pub masks: Vec<u8>,
    pub first_base: u8,
}

#[derive(Debug, Clone)]
pub struct CompiledMotif {
    pub name: String,
    pub sequence: String,
    pub forward: Pattern,
    pub reverse: Option<Pattern>,
    pub is_palindrome: bool,
}

impl CompiledMotif {
    pub fn len(&self) -> usize {
        self.forward.sequence.len()
    }
}

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

pub fn load_motif_file(path: &Path) -> Result<Vec<RawMotif>> {
    let reader = BufReader::new(
        File::open(path)
            .with_context(|| format!("failed to open motif file {}", path.display()))?,
    );
    let mut motifs = Vec::new();

    for (line_number, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let columns = split_columns(trimmed);
        if is_header(&columns) {
            continue;
        }

        let raw = match columns.as_slice() {
            [name, sequence] => RawMotif {
                name: name.to_string(),
                sequence: sequence.to_string(),
            },
            _ => bail!(
                "invalid motif file format at line {} in {}",
                line_number + 1,
                path.display()
            ),
        };

        motifs.push(raw);
    }

    if motifs.is_empty() {
        bail!("motif file {} did not contain any motifs", path.display())
    }

    Ok(motifs)
}

pub fn compile_motifs(
    raw_motifs: &[RawMotif],
    include_revcomp: bool,
    allow_iupac: bool,
) -> Result<Vec<CompiledMotif>> {
    raw_motifs
        .iter()
        .map(|raw| compile_motif(raw, include_revcomp, allow_iupac))
        .collect()
}

pub fn compile_motif(raw: &RawMotif, include_revcomp: bool, allow_iupac: bool) -> Result<CompiledMotif> {
    let sequence = normalize_sequence(&raw.sequence);
    if sequence.is_empty() {
        bail!("motif '{}' has an empty sequence", raw.name)
    }

    validate_motif_sequence(&raw.name, &sequence, allow_iupac)?;

    let reverse_sequence = reverse_complement(&sequence)?;
    let is_palindrome = sequence == reverse_sequence;

    let forward = Pattern {
        first_base: sequence[0],
        masks: sequence.iter().copied().map(iupac_mask).collect(),
        sequence: sequence.clone(),
    };
    let reverse = if include_revcomp && !is_palindrome {
        Some(Pattern {
            first_base: reverse_sequence[0],
            masks: reverse_sequence.iter().copied().map(iupac_mask).collect(),
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

pub fn normalize_sequence(sequence: &str) -> Vec<u8> {
    sequence
        .trim()
        .as_bytes()
        .iter()
        .map(|byte| byte.to_ascii_uppercase())
        .collect()
}

pub fn iupac_mask(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => 0b0001,
        b'C' => 0b0010,
        b'G' => 0b0100,
        b'T' | b'U' => 0b1000,
        b'R' => 0b0101,
        b'Y' => 0b1010,
        b'S' => 0b0110,
        b'W' => 0b1001,
        b'K' => 0b1100,
        b'M' => 0b0011,
        b'B' => 0b1110,
        b'D' => 0b1101,
        b'H' => 0b1011,
        b'V' => 0b0111,
        b'N' => 0b1111,
        _ => 0,
    }
}

pub fn canonical_base_mask(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => 0b0001,
        b'C' => 0b0010,
        b'G' => 0b0100,
        b'T' | b'U' => 0b1000,
        _ => 0,
    }
}

fn validate_motif_sequence(name: &str, sequence: &[u8], allow_iupac: bool) -> Result<()> {
    for base in sequence {
        let upper = base.to_ascii_uppercase();
        if canonical_base_mask(upper) != 0 {
            continue;
        }

        if iupac_mask(upper) != 0 {
            if allow_iupac {
                continue;
            }
            bail!(
                "motif '{}' contains IUPAC base '{}' but --iupac was not enabled",
                name,
                upper as char
            )
        }

        bail!("motif '{}' contains unsupported base '{}'", name, upper as char)
    }

    Ok(())
}

pub fn reverse_complement(sequence: &[u8]) -> Result<Vec<u8>> {
    sequence
        .iter()
        .rev()
        .map(|base| complement(*base))
        .collect()
}

fn complement(base: u8) -> Result<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Ok(b'T'),
        b'C' => Ok(b'G'),
        b'G' => Ok(b'C'),
        b'T' | b'U' => Ok(b'A'),
        b'R' => Ok(b'Y'),
        b'Y' => Ok(b'R'),
        b'S' => Ok(b'S'),
        b'W' => Ok(b'W'),
        b'K' => Ok(b'M'),
        b'M' => Ok(b'K'),
        b'B' => Ok(b'V'),
        b'D' => Ok(b'H'),
        b'H' => Ok(b'D'),
        b'V' => Ok(b'B'),
        b'N' => Ok(b'N'),
        other => Err(anyhow!("unsupported motif character '{}'", other as char)),
    }
}

fn split_columns(line: &str) -> Vec<&str> {
    if line.contains('\t') {
        line.split('\t').map(str::trim).collect()
    } else if line.contains(',') {
        line.split(',').map(str::trim).collect()
    } else {
        line.split_whitespace().collect()
    }
}

fn is_header(columns: &[&str]) -> bool {
    match columns {
        [name, sequence] => {
            name.eq_ignore_ascii_case("name") && sequence.eq_ignore_ascii_case("sequence")
        }
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::{canonical_base_mask, compile_motif, iupac_mask, reverse_complement, RawMotif};

    #[test]
    fn reverse_complement_handles_iupac() {
        let rc = reverse_complement(b"ATGRYSWKMBDHVN").unwrap();
        assert_eq!(String::from_utf8(rc).unwrap(), "NBDHVKMWSRYCAT");
    }

    #[test]
    fn iupac_masks_match_expected_bits() {
        assert_eq!(iupac_mask(b'A'), 0b0001);
        assert_eq!(iupac_mask(b'R'), 0b0101);
        assert_eq!(iupac_mask(b'N'), 0b1111);
        assert_eq!(iupac_mask(b'Z'), 0);
    }

    #[test]
    fn compiles_palindromic_reverse_once() {
        let motif = RawMotif {
            name: "pal".to_string(),
            sequence: "ATAT".to_string(),
        };
        let compiled = compile_motif(&motif, true, false).unwrap();
        assert!(compiled.is_palindrome);
        assert!(compiled.reverse.is_none());
    }

    #[test]
    fn rejects_iupac_motif_without_flag() {
        let motif = RawMotif {
            name: "iupac".to_string(),
            sequence: "ATGRN".to_string(),
        };
        let error = compile_motif(&motif, false, false).unwrap_err();
        assert!(error
            .to_string()
            .contains("contains IUPAC base 'R' but --iupac was not enabled"));
    }

    #[test]
    fn canonical_read_mask_rejects_ambiguous_bases() {
        assert_eq!(canonical_base_mask(b'A'), 0b0001);
        assert_eq!(canonical_base_mask(b'N'), 0);
        assert_eq!(canonical_base_mask(b'R'), 0);
    }
}
