//! Command-line interface module that parses user input into typed configuration.

use clap::{ArgAction, Args, Parser, Subcommand};

/// Top-level CLI object containing global flags and subcommands.
#[derive(Debug, Clone, Parser)]
#[command(
    name = "motifscan",
    version,
    disable_version_flag = true,
    about = "Streaming motif scanner for FASTA/FASTQ reads"
)]
pub struct Cli {
    #[arg(
        short = 'v',
        long = "version",
        action = ArgAction::SetTrue,
        global = true,
        help = "Print version and citation information"
    )]
    pub version_info: bool,

    #[arg(long, global = true, help = "Enable info-level logs")]
    pub verbose: bool,

    #[arg(long, global = true, help = "Enable debug-level logs")]
    pub debug: bool,

    #[command(subcommand)]
    pub command: Option<Command>,
}

impl Cli {
    /// Returns the worker-thread count for the current invocation, falling back to CPU count when no subcommand is selected.
    pub fn threads(&self) -> usize {
        match &self.command {
            Some(Command::Count(args)) => args.threads,
            None => num_cpus::get(),
        }
    }

    /// Returns the default log level derived from the CLI flags.
    pub fn log_level(&self) -> &'static str {
        if self.debug {
            "debug"
        } else if self.verbose {
            "info"
        } else {
            "warn"
        }
    }
}

/// Builds the version and citation banner printed by `-v/--version`.
pub fn version_banner() -> String {
    format!(
        "motifscan {}\n\nCitation (BibTeX):\n@software{{motifscan,\n  author = {{jiehua1995}},\n  title = {{MotifScan}},\n  url = {{https://github.com/jiehua1995/MotifScan}},\n  version = {{{}}}\n}}",
        env!("CARGO_PKG_VERSION"),
        env!("CARGO_PKG_VERSION"),
    )
}

/// Set of supported subcommands; currently only `count` is implemented.
#[derive(Debug, Clone, Subcommand)]
pub enum Command {
    #[command(about = "Count exact motif hits in reads", long_about = None)]
    Count(CountArgs),
}

/// Full argument set for the `count` subcommand, including input, motif, threading, and output targets.
#[derive(Debug, Clone, Args)]
#[command(about = "Count exact motif hits in reads", long_about = None)]
pub struct CountArgs {
    #[arg(
        short = 'i',
        long,
        help = "Input read file in FASTA, FASTQ, FASTA.GZ, or FASTQ.GZ format"
    )]
    pub input: std::path::PathBuf,
    #[arg(
        long,
        conflicts_with = "motifs",
        help = "Single motif sequence provided on the command line"
    )]
    pub motif: Option<String>,
    #[arg(long, default_value = "motif", help = "Output name used with --motif")]
    pub motif_name: String,
    #[arg(
        long,
        conflicts_with = "motif",
        help = "Two-column CSV file containing motif name and sequence"
    )]
    pub motifs: Option<std::path::PathBuf>,
    #[arg(long, help = "Also scan the reverse complement of each motif")]
    pub revcomp: bool,
    #[arg(short = 't', long, default_value_t = num_cpus::get(), help = "Number of worker threads to use")]
    pub threads: usize,
    #[arg(long, help = "Show a live progress display on stderr")]
    pub progress: bool,
    #[arg(short = 'o', long, help = "Output CSV file for motif summary counts")]
    pub output: std::path::PathBuf,
    #[arg(long, help = "Optional CSV file for read-level hit details")]
    pub report_read_hits: Option<std::path::PathBuf>,
}

impl CountArgs {
    /// Validates argument combinations, ensuring a motif source exists and thread count is non-zero.
    pub fn validate(&self) -> anyhow::Result<()> {
        if self.motif.is_none() && self.motifs.is_none() {
            anyhow::bail!("one of --motif or --motifs is required")
        }
        if self.threads == 0 {
            anyhow::bail!("--threads must be greater than 0")
        }
        Ok(())
    }
}
