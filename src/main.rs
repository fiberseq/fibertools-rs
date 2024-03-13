use anyhow::{bail, Error, Ok, Result};
use bio_io::*;
use colored::Colorize;
use env_logger::{Builder, Target};
use fibertools_rs::cli::*;
use fibertools_rs::nucleosomes::add_nucleosomes_to_bam;
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

    #[cfg(feature = "tch")]
    tch::set_num_threads(args.threads.try_into().unwrap());

    match &args.command {
        Some(Commands::Extract(extract_opts)) => {
            // read in the bam from stdin or from a file
            let mut bam = bam_reader(&extract_opts.bam, args.threads);
            extract::extract_contained(&mut bam, extract_opts);
        }
        Some(Commands::Center(center_opts)) => {
            // read in the bam from stdin or from a file
            let mut bam = bam::IndexedReader::from_path(center_opts.bam.clone())?;
            bam.set_threads(args.threads).unwrap();
            center::center_fiberdata(center_opts, &mut bam)?;
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

            let mut bam = bam_reader(&predict_m6a_opts.bam, args.threads);
            let mut out = bam_writer(&predict_m6a_opts.out, &bam, args.threads);
            predict_m6a::read_bam_into_fiberdata(&mut bam, &mut out, predict_m6a_opts);
        }
        Some(Commands::ClearKinetics { bam, out }) => {
            let mut bam = bam_reader(bam, args.threads);
            let mut out = bam_writer(out, &bam, args.threads);
            fibertools_rs::clear_kinetics(&mut bam, &mut out);
        }
        Some(Commands::StripBasemods { bam, out, basemod }) => {
            let mut bam = bam_reader(bam, args.threads);
            let mut out = bam_writer(out, &bam, args.threads);
            fibertools_rs::strip_basemods::strip_base_mods(&mut bam, &mut out, basemod);
        }
        Some(Commands::AddNucleosomes(nuc_opts)) => {
            add_nucleosomes_to_bam(nuc_opts, args.threads);
        }
        Some(Commands::Fire(fire_opts)) => {
            fibertools_rs::fire::add_fire_to_bam(fire_opts)?;
        }
        Some(Commands::Footprint(footprint_opts)) => {
            fibertools_rs::footprint::start_finding_footprints(footprint_opts)?;
        }
        Some(Commands::TrackDecorators(decorator_opts)) => {
            fibertools_rs::decorator::get_decorators_from_bam(decorator_opts)?;
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

fn _help_check_torch_env_vars() -> Result<()> {
    //std::env::set_var("LIBTORCH", "");
    //std::env::remove_var("LIBTORCH");
    //std::env::set_var("LD_LIBRARY_PATH", "aa");
    let lib = std::env::var("LIBTORCH")?;
    let ld_path = std::env::var("LD_LIBRARY_PATH")?;
    let dy_path = std::env::var("DYLD_LIBRARY_PATH")?;
    // make sure vars are not empty
    for v in [&lib, &ld_path, &dy_path] {
        if v.is_empty() {
            bail!("\t{}\n\t{}\n\t{}", lib, ld_path, dy_path)
        }
    }
    // make sure paths contain libtorch
    for path in [&ld_path, &dy_path] {
        if !path.contains(&lib) {
            bail!("\t{}\n\t{}\n\t{}", lib, ld_path, dy_path)
        }
    }
    Ok(())
}

fn _check_torch_env_vars() -> Result<()> {
    if let Err(e) = _help_check_torch_env_vars() {
        let msg = "\nEnvironment variable LIBTORCH is missing or is missing from LD_LIBRARY_PATH or DYLD_LIBRARY_PATH.".red().bold();
        let download_msg = "
Please download libtorch v1.13.0 from https://pytorch.org/get-started/ and extract the content of the zip file.

On a linux system you can download:
    wget https://download.pytorch.org/libtorch/cu116/libtorch-cxx11-abi-shared-with-deps-1.13.0%2Bcu116.zip
    
Then add the following to your .bashrc or equivalent:
    export LIBTORCH=/path/to/libtorch
    export LD_LIBRARY_PATH=${{LIBTORCH}}/lib:$LD_LIBRARY_PATH
    export DYLD_LIBRARY_PATH=${{LIBTORCH}}/lib:$LD_LIBRARY_PATH
            ".bold();
        let env_vars = format!("\nYour library variables:\n{}", e).yellow();
        bail!("{}\n{}\n{}", msg, env_vars, download_msg)
    };
    Ok(())
}
