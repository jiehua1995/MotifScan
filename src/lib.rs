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

/// Runs the main application flow.
///
/// This function parses CLI arguments, handles version/help output, initializes the Rayon thread pool, and finally dispatches to the selected subcommand.
pub fn run() -> Result<()> {
    let cli = cli::Cli::parse();
    if cli.version_info {
        println!("{}", cli::version_banner());
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
        cli::Command::Count(args) => scanner::run_count(&args),
    }
}

fn init_logging(cli: &cli::Cli) {
    let filter = if cli.debug {
        EnvFilter::new("debug")
    } else if cli.verbose {
        EnvFilter::new("info")
    } else {
        EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("warn"))
    };

    let _ = tracing_subscriber::fmt()
        .with_env_filter(filter)
        .with_target(false)
        .try_init();
}
