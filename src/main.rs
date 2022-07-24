use colored::Colorize;
use env_logger::{Builder, Target};
use fibertools_rs::cli::Commands;
use fibertools_rs::*;
use log::LevelFilter;
use rayon::prelude::*;
use rust_htslib::{bam, bam::Read};
use std::collections::HashMap;
use std::time::Instant;

fn main() {
    parse_cli();
}

pub fn parse_cli() {
    let pg_start = Instant::now();
    let args = cli::make_cli_parse();
    let matches = cli::make_cli_app().get_matches();
    let subcommand = matches.subcommand_name().unwrap();

    // set the logging level
    let min_log_level = match matches.occurrences_of("verbose") {
        0 => LevelFilter::Warn,
        1 => LevelFilter::Info,
        2 => LevelFilter::Debug,
        _ => LevelFilter::Trace,
    };
    Builder::new()
        .target(Target::Stderr)
        .filter(None, min_log_level)
        .init();

    log::debug!("DEBUG logging enabled");
    log::trace!("TRACE logging enabled");

    // set up number of threads to use globally
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    match &args.command {
        //
        // Run Stats
        //
        Some(Commands::Extract { bam, reference }) => {
            let mut bam = bam::Reader::from_path(bam).unwrap();
            let _header = bam::Header::from_template(bam.header());
            // copy reverse reads to new BAM file
            for r in bam.records() {
                let record = r.unwrap();
                if record.is_reverse() {
                    println!("{:?}", record.qname());
                }
            }
        }
        //
        // no command opt
        //
        None => {}
    };
    let duration = pg_start.elapsed();
    log::info!(
        "{} done! Time elapsed: {}",
        subcommand.bright_green().bold(),
        format!("{:.2?}", duration).bright_yellow().bold()
    );
}
