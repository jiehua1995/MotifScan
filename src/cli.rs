//! Command-line interface module that parses user input into typed configuration.

use clap::{ArgAction, Args, Parser, Subcommand};

#[derive(Debug, Clone, Parser)]
#[command(
    name = "motifscan",
    version,
    disable_version_flag = true,
    disable_help_flag = true,
    disable_help_subcommand = true,
    about = "Streaming motif scanner for FASTA/FASTQ reads"
)]
pub struct Cli {
    #[arg(short = 'v', long = "version", action = ArgAction::SetTrue, global = true, help = "Print version and citation information", help_heading = "Info")]
    pub version_info: bool,
    #[arg(short = 'h', long = "help", help = "Print help", action = ArgAction::SetTrue, help_heading = "Info")]
    pub help: bool,
    #[command(subcommand)]
    pub command: Option<Command>,
}

impl Cli {
    pub fn threads(&self) -> usize {
        match &self.command {
            Some(Command::Count(args)) => args.threads,
            Some(Command::Species(args)) => args.threads,
            None => num_cpus::get(),
        }
    }

    pub fn log_level(&self) -> &'static str { "warn" }
}

pub fn version_banner() -> String {
    format!(
        "motifscan {}\n\nCitation (BibTeX):\n@software{{motifscan,\n  author = {{jiehua1995}},\n  title = {{MotifScan}},\n  url = {{https://github.com/jiehua1995/MotifScan}},\n  version = {{{}}}\n}}",
        env!("CARGO_PKG_VERSION"),
        env!("CARGO_PKG_VERSION"),
    )
}

#[derive(Debug, Clone, Subcommand)]
pub enum Command {
    #[command(about = "Count exact motif hits in reads", long_about = None)]
    Count(CountArgs),
    #[command(about = "Fuzzy long-read locus matching with automatic mel/sim SNP voting", long_about = None)]
    Species(SpeciesArgs),
}

#[derive(Debug, Clone, Args)]
#[command(
    about = "Count exact motif hits in reads",
    long_about = None,
    after_help = "\n Examples:\n  motifscan count --motif ATGCGACCGATGCGTASGGC -i reads.fq -o out.csv\n  motifscan count --motifs motifs.csv -i reads.fq -o out.csv --revcomp",
)]
pub struct CountArgs {
    #[arg(short = 'i', long, help = "Input read file in FASTA, FASTQ, FASTA.GZ, or FASTQ.GZ format", help_heading = "File")]
    pub input: Option<std::path::PathBuf>,
    #[arg(long, conflicts_with = "motifs", help = "Single motif sequence provided on the command line", help_heading = "Motif")]
    pub motif: Option<String>,
    #[arg(long, default_value = "motif", help = "Output name used with --motif", help_heading = "Motif", requires = "motif")]
    pub motif_name: String,
    #[arg(long, conflicts_with = "motif", help = "Two-column CSV file containing motif name and sequence", help_heading = "Motif")]
    pub motifs: Option<std::path::PathBuf>,
    #[arg(long, help = "Also scan the reverse complement of each motif", help_heading = "Motif")]
    pub revcomp: bool,
    #[arg(short = 't', long, default_value_t = num_cpus::get(), help = "Number of worker threads to use", help_heading = "Performance")]
    pub threads: usize,
    #[arg(long, help = "Show a live progress display on stderr", help_heading = "Behavior")]
    pub progress: bool,
    #[arg(short = 'h', long = "help", help = "Print help", action = ArgAction::SetTrue, help_heading = "Info")]
    pub help: bool,
    #[arg(short = 'o', long, help = "Output CSV file for motif summary counts", help_heading = "File")]
    pub output: Option<std::path::PathBuf>,
    #[arg(long, help = "Optional CSV file for read-level hit details", help_heading = "Behavior")]
    pub report_read_hits: Option<std::path::PathBuf>,
}

impl CountArgs {
    pub fn validate(&self) -> anyhow::Result<()> {
        if self.help { return Ok(()); }
        if self.motif.is_none() && self.motifs.is_none() { anyhow::bail!("one of --motif or --motifs is required") }
        if self.threads == 0 { anyhow::bail!("--threads must be greater than 0") }
        if self.input.is_none() { anyhow::bail!("--input is required") }
        if self.output.is_none() { anyhow::bail!("--output is required") }
        Ok(())
    }
}

#[derive(Debug, Clone, Args)]
#[command(
    about = "Fuzzy long-read locus matching with automatic species-diagnostic SNP voting",
    long_about = None,
    after_help = "\nExamples:\n  motifscan species -i reads.fastq.gz --sample HYC1 --motifs motifs.csv --pairs pairs.csv -o HYC1.csv -t 48 --progress\n  motifscan species --input-list samples.txt --motifs motifs.csv --pairs pairs.csv -o all_samples.csv -t 48 --progress",
)]
pub struct SpeciesArgs {
    #[arg(short = 'i', long, conflicts_with = "input_list", help = "Single FASTA/FASTQ(.gz) file", help_heading = "File")]
    pub input: Option<std::path::PathBuf>,
    #[arg(long, conflicts_with = "input", help = "TXT/TSV with path per line or sample<TAB>path", help_heading = "File")]
    pub input_list: Option<std::path::PathBuf>,
    #[arg(long, requires = "input", help = "Sample name for a single --input; otherwise inferred from filename", help_heading = "File")]
    pub sample: Option<String>,
    #[arg(long, help = "Two-column motif CSV: name,sequence", help_heading = "Reference")]
    pub motifs: Option<std::path::PathBuf>,
    #[arg(long, help = "Three-column pair CSV: locus,mel,sim", help_heading = "Reference")]
    pub pairs: Option<std::path::PathBuf>,
    #[arg(short = 'o', long, help = "Output sample-by-locus summary CSV", help_heading = "File")]
    pub output: Option<std::path::PathBuf>,
    #[arg(long, help = "Optional SNP-level QC CSV; default is <output stem>.snps.csv", help_heading = "File")]
    pub snp_output: Option<std::path::PathBuf>,
    #[arg(long, help = "Optional pair/SNP extraction QC CSV", help_heading = "File")]
    pub pair_qc_output: Option<std::path::PathBuf>,
    #[arg(short = 't', long, default_value_t = num_cpus::get(), help = "Number of Rayon worker threads", help_heading = "Performance")]
    pub threads: usize,
    #[arg(long, help = "Show per-file byte/read progress", help_heading = "Behavior")]
    pub progress: bool,
    #[arg(skip)]
    pub revcomp: bool,
    #[arg(long, default_value_t = 11, help = "Shared exact anchor k-mer length used only for candidate retrieval", help_heading = "Performance")]
    pub anchor_k: usize,
    #[arg(long, default_value_t = 8, help = "Number of evenly spaced shared anchors per locus", help_heading = "Performance")]
    pub anchors_per_locus: usize,
    #[arg(long, default_value_t = 20, help = "Extra read bases on each side of an anchor-estimated alignment window", help_heading = "Alignment")]
    pub alignment_slack: usize,
    #[arg(long, default_value_t = 0.85, help = "Minimum identity across non-diagnostic positions", help_heading = "Alignment")]
    pub min_shared_identity: f64,
    #[arg(long, default_value_t = 80, help = "Minimum number of aligned reference bases", help_heading = "Alignment")]
    pub min_aligned_bases: usize,
    #[arg(long, default_value_t = 15, help = "Minimum FASTQ Phred quality at a diagnostic SNP", help_heading = "Species calling")]
    pub min_snp_baseq: u8,
    #[arg(long, default_value_t = 2, help = "Minimum high-quality diagnostic SNPs required to call mel/sim", help_heading = "Species calling")]
    pub min_informative_snps: usize,
    #[arg(long, default_value_t = 0.75, help = "Fraction of informative SNPs supporting one species required for a call", help_heading = "Species calling")]
    pub species_fraction: f64,
    #[arg(long, value_parser = ["best", "all"], default_value = "best", help = "best: one best locus/read; all: keep all passing loci", help_heading = "Behavior")]
    pub locus_mode: String,
    #[arg(short = 'h', long = "help", help = "Print help", action = ArgAction::SetTrue, help_heading = "Info")]
    pub help: bool,
}

impl SpeciesArgs {
    pub fn validate(&self) -> anyhow::Result<()> {
        if self.help { return Ok(()); }
        if self.input.is_none() && self.input_list.is_none() { anyhow::bail!("one of --input or --input-list is required") }
        if self.motifs.is_none() { anyhow::bail!("--motifs is required") }
        if self.pairs.is_none() { anyhow::bail!("--pairs is required") }
        if self.output.is_none() { anyhow::bail!("--output is required") }
        if self.threads == 0 { anyhow::bail!("--threads must be greater than 0") }
        if !(0.0 < self.min_shared_identity && self.min_shared_identity <= 1.0) { anyhow::bail!("--min-shared-identity must be in (0,1]") }
        if !(0.5 < self.species_fraction && self.species_fraction <= 1.0) { anyhow::bail!("--species-fraction must be in (0.5,1]") }
        if self.min_informative_snps == 0 { anyhow::bail!("--min-informative-snps must be greater than 0") }
        Ok(())
    }
}
