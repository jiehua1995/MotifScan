use tempfile::tempdir;

#[test]
fn integration_run_count_on_sample_data() {
    let tmp = tempdir().unwrap();
    let out_count = tmp.path().join("count.csv");
    let out_hits = tmp.path().join("read_hits.csv");

    let args = motifscan::cli::CountArgs {
        input: std::path::PathBuf::from("test/reads_sample.fastq.gz"),
        motif: None,
        motif_name: "motif".to_string(),
        motifs: Some(std::path::PathBuf::from("test/motifs.csv")),
        revcomp: true,
        threads: 1,
        progress: false,
        output: out_count.clone(),
        report_read_hits: Some(out_hits.clone()),
    };

    motifscan::scanner::run_count(&args).expect("run_count failed");

    let count_txt = std::fs::read_to_string(out_count).unwrap();
    assert!(count_txt.contains("motif,sequence,length"));
    let hits_txt = std::fs::read_to_string(out_hits).unwrap();
    assert!(hits_txt.contains("read_id,motif,strand,position,matched_sequence"));
}
