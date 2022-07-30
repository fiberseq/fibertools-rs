use anyhow::{Error, Ok};
use colored::Colorize;
use env_logger::{Builder, Target};
use fibertools_rs::cli::Commands;
use fibertools_rs::*;
use log::LevelFilter;
use rust_htslib::{bam, bam::Read};
use std::time::Instant;

pub fn main() -> Result<(), Error> {
    colored::control::set_override(true);
    console::set_colors_enabled(true);
    let pg_start = Instant::now();
    let args = cli::make_cli_parse();
    let matches = cli::make_cli_app().get_matches();
    let subcommand = matches.subcommand_name().unwrap();

    // set the logging level
    let mut min_log_level = match matches.occurrences_of("verbose") {
        //0 => LevelFilter::Warn,
        0 => LevelFilter::Info,
        1 => LevelFilter::Debug,
        _ => LevelFilter::Trace,
    };
    if args.quiet {
        min_log_level = LevelFilter::Error;
    }

    Builder::new()
        .target(Target::Stderr)
        .filter(None, min_log_level)
        .init();

    log::debug!("DEBUG logging enabled");
    log::trace!("TRACE logging enabled");

    log::info!("Starting ft-{}", subcommand.bright_green().bold());

    // set up number of threads to use globally
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    match &args.command {
        Some(Commands::Extract {
            bam,
            reference,
            m6a,
            cpg,
            msp,
            nuc,
            all,
        }) => {
            // read in the bam from stdin or from a file
            let mut bam = if bam == "-" {
                bam::Reader::from_stdin().unwrap()
            } else {
                bam::Reader::from_path(bam).unwrap_or_else(|_| panic!("Failed to open {}", bam))
            };
            bam.set_threads(args.threads).unwrap();
            let out_files = FiberOutFiles::new(m6a, cpg, msp, nuc, all)?;
            extract::extract_contained(&mut bam, *reference, out_files);
        }
        Some(Commands::Center { bam, bed, wide }) => {
            // read in the bam from stdin or from a file
            let mut bam = bam::IndexedReader::from_path(bam)?;
            bam.set_threads(args.threads).unwrap();
            let center_positions = center::read_center_positions(bed)?;
            center::center_fiberdata(&mut bam, center_positions, *wide);
        }
        None => {}
    };
    let duration = pg_start.elapsed();
    log::info!(
        "{} done! Time elapsed: {}",
        subcommand.bright_green().bold(),
        format!("{:.2?}", duration).bright_yellow().bold()
    );
    Ok(())
}
