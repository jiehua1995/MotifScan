fn main() {
    if let Err(error) = motifscan::run() {
        eprintln!("error: {error:#}");
        std::process::exit(1);
    }
}
