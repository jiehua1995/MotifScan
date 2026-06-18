//! Library entry module that wires together CLI parsing, thread-pool setup, and scan execution.

pub mod cli;
pub mod io;
pub mod motif;
pub mod output;
pub mod scanner;

use anyhow::Result;
use clap::CommandFactory;
use clap::Parser;
use tracing_subscriber::EnvFilter;

fn print_top_help() {
    println!("motifscan — Streaming motif scanner for FASTA/FASTQ reads\n");
    println!("Usage: motifscan <COMMAND> [OPTIONS]\n");
    println!("Commands:");
    println!("  count    Count exact motif hits in reads\n");
    println!("Info:");
    println!("  -v, --version    Print version and citation information");
    println!("  -h, --help       Print help\n");
}

fn print_count_help() {
    println!("Count exact motif hits in reads\n");
    println!("Usage: count [OPTIONS]\n");
    println!("Required:");
    println!("  -i, --input <INPUT>    Input read file in FASTA, FASTQ, FASTA.GZ, or FASTQ.GZ format  (required)");
    println!("  -o, --output <OUTPUT>  Output CSV file for motif summary counts  (required)\n");
    println!("Motif:");
    println!("      --motif <MOTIF>            Single motif sequence provided on the command line");
    println!("      --motif-name <MOTIF_NAME>  Output name used with --motif (requires --motif)");
    println!(
        "      --motifs <MOTIFS>          Two-column CSV file containing motif name and sequence"
    );
    println!("      --revcomp                  Also scan the reverse complement of each motif\n");
    println!("Performance:");
    println!("  -t, --threads <THREADS>  Number of worker threads to use [default: auto]\n");
    println!("Behavior:");
    println!("      --progress                             Show a live progress display on stderr");
    println!("      --report-read-hits <REPORT_READ_HITS>  Optional CSV file for read-level hit details\n");
    println!("Info:");
    println!("  -h, --help  Print help");

    println!("Examples:");
    println!("  motifscan count --motif ATGCGACCGATGCGTASGGC -i reads.fq -o out.csv");
    println!("  motifscan count --motifs motifs.csv -i reads.fq -o out.csv --revcomp");
}

/// Runs the main application flow.
///
/// This function parses CLI arguments, handles version/help output, initializes the Rayon thread pool, and finally dispatches to the selected subcommand.
pub fn run() -> Result<()> {
    let cli = cli::Cli::parse();
    if cli.version_info {
        println!("{}", cli::version_banner());
        return Ok(());
    }
    if cli.help {
        print_top_help();
        return Ok(());
    }

    init_logging(&cli);

    let thread_count = cli.threads();

    let Some(command) = cli.command else {
        let mut help = cli::Cli::command();
        help.print_help()?;
        println!();
        return Ok(());
    };

    rayon::ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build_global()
        .ok();

    match command {
        cli::Command::Count(args) => {
            if args.help {
                print_count_help();
                return Ok(());
            }
            scanner::run_count(&args)
        }
    }
}

fn init_logging(_cli: &cli::Cli) {
    let filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("warn"));

    let _ = tracing_subscriber::fmt()
        .with_env_filter(filter)
        .with_target(false)
        .try_init();
}
