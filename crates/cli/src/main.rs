fn main() -> std::process::ExitCode {
    bammetrics_cli::remove_files_on_signal();
    bammetrics_cli::main()
}
