use std::fs;
use std::process::Command;

use tempfile::tempdir;

fn pseudo_dna(len: usize) -> String {
    let alphabet = [b'A', b'C', b'G', b'T'];
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut out = Vec::with_capacity(len);
    for _ in 0..len {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        out.push(alphabet[((state >> 32) & 3) as usize]);
    }
    String::from_utf8(out).unwrap()
}

fn alternate(base: u8) -> u8 {
    match base {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        b'T' => b'A',
        _ => unreachable!(),
    }
}

fn reverse_complement(seq: &str) -> String {
    seq.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => 'N',
        })
        .collect()
}

fn fastq_record(id: &str, seq: &str, q: char) -> String {
    let qual: String = std::iter::repeat_n(q, seq.len()).collect();
    format!("@{id}\n{seq}\n+\n{qual}\n")
}

#[test]
fn species_cli_classifies_synthetic_long_reads_end_to_end() {
    let dir = tempdir().unwrap();
    let motifs = dir.path().join("motifs.csv");
    let pairs = dir.path().join("pairs.csv");
    let reads = dir.path().join("reads.fastq");
    let output = dir.path().join("summary.csv");
    let snps = dir.path().join("snps.csv");
    let pair_qc = dir.path().join("pairs_qc.csv");

    let mel = pseudo_dna(140);
    let mut sim_bytes = mel.as_bytes().to_vec();
    for pos in [20usize, 50, 80, 110] {
        sim_bytes[pos] = alternate(sim_bytes[pos]);
    }
    let sim = String::from_utf8(sim_bytes).unwrap();

    fs::write(
        &motifs,
        format!("name,sequence\nmel_ref,{mel}\nsim_ref,{sim}\n"),
    )
    .unwrap();
    fs::write(&pairs, "locus,mel,sim\ntoy,mel_ref,sim_ref\n").unwrap();

    let mut mel_error = mel.as_bytes().to_vec();
    mel_error[7] = alternate(mel_error[7]);
    mel_error.insert(72, b'T');
    let mel_error = String::from_utf8(mel_error).unwrap();
    let sim_rc = reverse_complement(&sim);

    let fastq = [
        fastq_record("mel_exact", &mel, 'I'),
        fastq_record("mel_shared_errors", &mel_error, 'I'),
        fastq_record("sim_exact", &sim, 'I'),
        fastq_record("sim_revcomp", &sim_rc, 'I'),
    ]
    .concat();
    fs::write(&reads, fastq).unwrap();

    let status = Command::new(env!("CARGO_BIN_EXE_motifscan"))
        .arg("species")
        .arg("--input")
        .arg(&reads)
        .arg("--sample")
        .arg("synthetic")
        .arg("--motifs")
        .arg(&motifs)
        .arg("--pairs")
        .arg(&pairs)
        .arg("--output")
        .arg(&output)
        .arg("--snp-output")
        .arg(&snps)
        .arg("--pair-qc-output")
        .arg(&pair_qc)
        .arg("--min-shared-identity")
        .arg("0.85")
        .arg("--min-aligned-bases")
        .arg("100")
        .arg("--min-snp-baseq")
        .arg("15")
        .arg("--min-informative-snps")
        .arg("2")
        .arg("--species-fraction")
        .arg("0.75")
        .arg("--anchor-k")
        .arg("11")
        .arg("--anchors-per-locus")
        .arg("8")
        .arg("--alignment-slack")
        .arg("20")
        .status()
        .unwrap();
    assert!(status.success());

    let mut reader = csv::Reader::from_path(&output).unwrap();
    let headers = reader.headers().unwrap().clone();
    let row = reader.records().next().unwrap().unwrap();

    let col = |name: &str| headers.iter().position(|x| x == name).unwrap();
    assert_eq!(&row[col("sample")], "synthetic");
    assert_eq!(&row[col("locus")], "toy");
    assert_eq!(&row[col("reads_with_hit")], "4");
    assert_eq!(&row[col("mel_reads")], "2");
    assert_eq!(&row[col("sim_reads")], "2");
    assert_eq!(&row[col("species_assigned_reads")], "4");
    assert_eq!(&row[col("diagnostic_snp_count")], "4");

    let snp_text = fs::read_to_string(&snps).unwrap();
    assert!(snp_text.lines().count() >= 5);
    let pair_text = fs::read_to_string(&pair_qc).unwrap();
    assert!(pair_text.contains("toy,mel_ref,sim_ref,140,140,4"));
}
