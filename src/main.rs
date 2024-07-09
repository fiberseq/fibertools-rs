use anyhow::{Error, Ok, Result};
use colored::Colorize;
use env_logger::{Builder, Target};
use fibertools_rs::cli::*;
use fibertools_rs::nucleosomes::add_nucleosomes_to_bam;
use fibertools_rs::*;
use log::LevelFilter;
use std::time::Instant;

pub fn main() -> Result<(), Error> {
    colored::control::set_override(true);
    console::set_colors_enabled(true);
    let pg_start = Instant::now();
    let mut args = cli::make_cli_parse();
    let matches = cli::make_cli_app().get_matches();
    let subcommand = matches.subcommand_name().unwrap();

    // set the logging level
    let mut min_log_level = match matches.get_count("verbose") {
        0 => LevelFilter::Warn,
        1 => LevelFilter::Info,
        2 => LevelFilter::Debug,
        _ => LevelFilter::Trace,
    };
    if args.global.quiet {
        min_log_level = LevelFilter::Error;
    }

    Builder::new()
        .target(Target::Stderr)
        .filter(None, min_log_level)
        .init();

    log::debug!("DEBUG logging enabled");
    log::trace!("TRACE logging enabled");
    log::debug!("Using {} threads.", args.global.threads);
    log::info!("Starting ft-{}", subcommand.bright_green().bold());

    // set up number of threads to use globally
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.global.threads)
        .build_global()
        .unwrap();

    #[cfg(feature = "tch")]
    {
        // setting this to 1 since I do paralyzation via processing multiple reads
        tch::set_num_threads(1);
        tch::set_num_interop_threads(1);
    }

    log::debug!("Command line options: {:?}", args.command);

    match &mut args.command {
        Some(Commands::Extract(extract_opts)) => {
            extract::extract_contained(extract_opts);
        }
        Some(Commands::Center(center_opts)) => {
            // read in the bam from stdin or from a file
            center::center_fiberdata(center_opts)?;
        }
        #[allow(unused)]
        Some(Commands::PredictM6A(predict_m6a_opts)) => {
            #[cfg(not(feature = "tch"))]
            log::warn!(
                "\n{}\n\t{}\n\t{}\n", 
                "WARNING: m6A predictions are slower without the pytorch backend.".bright_yellow().bold(),
                "Consider recompiling via cargo with: `--all-features`.",
                "For detailed instructions see: https://fiberseq.github.io/fibertools-rs/INSTALL.html."
            );

            predict_m6a::read_bam_into_fiberdata(predict_m6a_opts);
        }
        Some(Commands::ClearKinetics(clear_kinetics_opts)) => {
            fibertools_rs::clear_kinetics(clear_kinetics_opts);
        }
        Some(Commands::StripBasemods(strip_basemods_opts)) => {
            fibertools_rs::strip_basemods::strip_base_mods(strip_basemods_opts);
        }
        Some(Commands::AddNucleosomes(nuc_opts)) => {
            add_nucleosomes_to_bam(nuc_opts);
        }
        Some(Commands::Fire(fire_opts)) => {
            fibertools_rs::fire::add_fire_to_bam(fire_opts)?;
        }
        Some(Commands::Footprint(footprint_opts)) => {
            fibertools_rs::footprint::start_finding_footprints(footprint_opts)?;
        }
        Some(Commands::Pileup(pileup_opts)) => {
            fibertools_rs::pileup::pileup_track(pileup_opts)?;
        }
        Some(Commands::TrackDecorators(decorator_opts)) => {
            fibertools_rs::decorator::get_decorators_from_bam(decorator_opts)?;
        }
        Some(Commands::Completions(completion_opts)) => {
            log::info!(
                "Generating completion file for {:?}...",
                completion_opts.shell
            );
            cli::print_completions(completion_opts.shell, &mut cli::make_cli_app());
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
