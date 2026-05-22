/// Binary entrypoint that simply calls the library runner and normalizes fatal error handling.
fn main() {
    if let Err(error) = motifscan::run() {
        eprintln!("error: {error:#}");
        std::process::exit(1);
    }
}
