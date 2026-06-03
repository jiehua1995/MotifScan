//! Command-line interface module that parses user input into typed configuration.

use clap::{ArgAction, Args, Parser, Subcommand};

/// Top-level CLI object containing global flags and subcommands.
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
    #[arg(
        short = 'v',
        long = "version",
        action = ArgAction::SetTrue,
        global = true,
        help = "Print version and citation information",
        help_heading = "Info"
    )]
    pub version_info: bool,
    #[arg(short = 'h', long = "help", help = "Print help", action = ArgAction::SetTrue, help_heading = "Info")]
    pub help: bool,

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
        "warn"
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
        help = "Input read file in FASTA, FASTQ, FASTA.GZ, or FASTQ.GZ format",
        help_heading = "File"
    )]
    pub input: Option<std::path::PathBuf>,
    #[arg(
        long,
        conflicts_with = "motifs",
        help = "Single motif sequence provided on the command line",
        help_heading = "Motif"
    )]
    pub motif: Option<String>,
    #[arg(
        long,
        default_value = "motif",
        help = "Output name used with --motif",
        help_heading = "Motif",
        requires = "motif"
    )]
    pub motif_name: String,
    #[arg(
        long,
        conflicts_with = "motif",
        help = "Two-column CSV file containing motif name and sequence",
        help_heading = "Motif"
    )]
    pub motifs: Option<std::path::PathBuf>,
    #[arg(
        long,
        help = "Also scan the reverse complement of each motif",
        help_heading = "Motif"
    )]
    pub revcomp: bool,
    #[arg(short = 't', long, default_value_t = num_cpus::get(), help = "Number of worker threads to use", help_heading = "Performance")]
    pub threads: usize,
    #[arg(
        long,
        help = "Show a live progress display on stderr",
        help_heading = "Behavior"
    )]
    pub progress: bool,
    #[arg(short = 'h', long = "help", help = "Print help", action = ArgAction::SetTrue, help_heading = "Info")]
    pub help: bool,
    #[arg(
        short = 'o',
        long,
        help = "Output CSV file for motif summary counts",
        help_heading = "File"
    )]
    pub output: Option<std::path::PathBuf>,
    #[arg(
        long,
        help = "Optional CSV file for read-level hit details",
        help_heading = "Behavior"
    )]
    pub report_read_hits: Option<std::path::PathBuf>,
}

impl CountArgs {
    /// Validates argument combinations, ensuring a motif source exists and thread count is non-zero.
    pub fn validate(&self) -> anyhow::Result<()> {
        if self.help {
            return Ok(());
        }
        if self.motif.is_none() && self.motifs.is_none() {
            anyhow::bail!("one of --motif or --motifs is required")
        }
        if self.threads == 0 {
            anyhow::bail!("--threads must be greater than 0")
        }
        if self.input.is_none() {
            anyhow::bail!("--input is required")
        }
        if self.output.is_none() {
            anyhow::bail!("--output is required")
        }

        Ok(())
    }
}
