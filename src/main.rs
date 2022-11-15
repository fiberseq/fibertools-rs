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
        0 => LevelFilter::Warn,
        1 => LevelFilter::Info,
        2 => LevelFilter::Debug,
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
    log::debug!("Using {} threads.", args.threads);
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
            simplify,
            quality,
            min_ml_score,
            m6a,
            cpg,
            msp,
            nuc,
            all,
        }) => {
            // read in the bam from stdin or from a file
            let mut bam = fibertools_rs::bam_reader(bam, args.threads);
            let out_files = FiberOut::new(
                m6a,
                cpg,
                msp,
                nuc,
                all,
                *reference,
                *simplify,
                *quality,
                *min_ml_score,
            )?;
            extract::extract_contained(&mut bam, out_files);
        }
        Some(Commands::Center {
            bam,
            bed,
            min_ml_score,
            wide,
        }) => {
            // read in the bam from stdin or from a file
            let mut bam = bam::IndexedReader::from_path(bam)?;
            bam.set_threads(args.threads).unwrap();
            let center_positions = center::read_center_positions(bed)?;
            center::center_fiberdata(&mut bam, center_positions, *min_ml_score, *wide);
        }
        Some(Commands::PredictM6A {
            bam,
            out,
            keep,
            cnn,
            full_float,
        }) => {
            let mut bam = fibertools_rs::bam_reader(bam, args.threads);
            let header = bam::Header::from_template(bam.header());
            let mut out = if out == "-" {
                bam::Writer::from_stdout(&header, bam::Format::Bam).unwrap()
            } else {
                bam::Writer::from_path(out, &header, bam::Format::Bam).unwrap()
            };
            out.set_threads(args.threads).unwrap();

            let predict_options = predict_m6a::PredictOptions {
                keep: *keep,
                cnn: *cnn,
                full_float: *full_float,
                polymerase: find_pb_polymerase(&header),
            };
            predict_m6a::read_bam_into_fiberdata(&mut bam, &mut out, &predict_options);
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
