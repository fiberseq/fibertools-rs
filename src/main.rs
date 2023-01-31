use anyhow::{Error, Ok};
use colored::Colorize;
use env_logger::{Builder, Target};
use fibertools_rs::cli::Commands;
#[cfg(feature = "predict")]
use fibertools_rs::predict_m6a::PredictOptions;
use fibertools_rs::*;
use log::LevelFilter;
use rust_htslib::{bam, bam::Read};
use std::time::Instant;

/*
use std::path::Path;
// TODO FINISH THIS CHECK
pub fn yank_check() {
    let config = cargo::util::Config::default().unwrap();
    //let build = config.build_config().unwrap();
    let cur_exe_path = config.cargo_exe().unwrap();
    log::warn!("{:?}", config.cargo_exe());
    let sid = cargo::core::SourceId::for_path(cur_exe_path).unwrap();
    log::warn!("{:?}", config.registry_source_path());
    let _ws = cargo::core::Workspace::new(
        Path::new("/Users/mrvollger/Desktop/repos/fibertools-rs/Cargo.toml"),
        &config,
    )
    .unwrap();

    let pkg = cargo::core::package_id::PackageId::new("fibertools-rs", "0.0.8", sid).unwrap();
    log::warn!("{:?}", pkg);
    let src_map = cargo::core::source::SourceMap::new();
    let _pkg_set = cargo::core::package::PackageSet::new(&[pkg], src_map, &config).unwrap();
}
*/

pub fn main() -> Result<(), Error> {
    colored::control::set_override(true);
    console::set_colors_enabled(true);
    let pg_start = Instant::now();
    let args = cli::make_cli_parse();
    let matches = cli::make_cli_app().get_matches();
    let subcommand = matches.subcommand_name().unwrap();

    // set the logging level
    let mut min_log_level = match matches.get_count("verbose") {
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
            full_float,
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
                *full_float,
            )?;
            extract::extract_contained(&mut bam, out_files);
        }
        Some(Commands::Center {
            bam,
            bed,
            min_ml_score,
            dist,
            wide,
        }) => {
            // read in the bam from stdin or from a file
            let mut bam = bam::IndexedReader::from_path(bam)?;
            bam.set_threads(args.threads).unwrap();
            let center_positions = center::read_center_positions(bed)?;
            center::center_fiberdata(&mut bam, center_positions, *min_ml_score, *wide, *dist);
        }
        #[cfg(feature = "predict")]
        Some(Commands::PredictM6A {
            bam,
            out,
            keep,
            min_ml_score,
            all_calls,
            xgb,
            cnn,
            semi,
            full_float,
            batch_size,
        }) => {
            // cnn must be set to true if using semi
            let cnn = if *semi { &true } else { cnn };
            // read bam
            let mut bam = fibertools_rs::bam_reader(bam, args.threads);
            let header = bam::Header::from_template(bam.header());
            let mut out = fibertools_rs::bam_writer(out, &bam, args.threads);
            let predict_options = PredictOptions::new(
                *keep,
                *xgb,
                *cnn,
                *semi,
                *full_float,
                *min_ml_score,
                *all_calls,
                find_pb_polymerase(&header),
                *batch_size,
            );

            log::info!("{} {} {}", xgb, cnn, semi);
            log::info!("{} reads included at once in batch prediction.", batch_size);
            predict_m6a::read_bam_into_fiberdata(&mut bam, &mut out, &predict_options);
        }
        Some(Commands::ClearKinetics { bam, out }) => {
            let mut bam = fibertools_rs::bam_reader(bam, args.threads);
            let mut out = fibertools_rs::bam_writer(out, &bam, args.threads);
            fibertools_rs::clear_kinetics(&mut bam, &mut out);
        }
        Some(Commands::Completions { shell }) => {
            log::info!("Generating completion file for {:?}...", shell);
            cli::print_completions(*shell, &mut cli::make_cli_app());
        }
        Some(Commands::Man {}) => {
            let man = clap_mangen::Man::new(cli::make_cli_app());
            //let mut buffer: Vec<u8> = Default::default();
            let mut buffer = Box::new(std::io::stdout()) as Box<dyn std::io::Write>;
            man.render(&mut buffer)?;
            //man.render_subcommands_section(&mut buffer)?;
            //man.render_description_section(&mut buffer)?;
            //std::fs::write("ft.1", buffer)?;
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
