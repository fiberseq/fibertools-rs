use anyhow::{bail, Error, Ok, Result};
use bio_io::*;
use colored::Colorize;
use env_logger::{Builder, Target};
#[cfg(feature = "predict")]
use fibertools_rs::predict_m6a::PredictOptions;
use fibertools_rs::*;
use fibertools_rs::{cli::Commands, nucleosomes::add_nucleosomes_to_bam};
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
            molecular: _molecular,
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

            // read in the bam from stdin or from a file
            let mut bam = bam_reader(bam, args.threads);
            extract::extract_contained(&mut bam, out_files);
        }
        Some(Commands::Center {
            bam,
            bed,
            min_ml_score,
            dist,
            wide,
            reference,
            simplify,
        }) => {
            // read in the bam from stdin or from a file
            let mut bam = bam::IndexedReader::from_path(bam)?;
            bam.set_threads(args.threads).unwrap();
            let center_positions = center::read_center_positions(bed)?;
            center::center_fiberdata(
                &mut bam,
                center_positions,
                *min_ml_score,
                *wide,
                *dist,
                *reference,
                *simplify,
            );
        }
        #[allow(unused)]
        Some(Commands::PredictM6A {
            bam,
            out,
            nucleosome_length,
            combined_nucleosome_length,
            min_distance_added,
            distance_from_end,
            allowed_m6a_skips,
            keep,
            min_ml_score,
            all_calls,
            xgb,
            cnn,
            semi,
            full_float,
            batch_size,
        }) => {
            //check_torch_env_vars()?;
            #[cfg(feature = "predict")]
            {
                // cnn must be set to true if using semi
                let cnn = if *semi { &true } else { cnn };
                // read bam
                let mut bam = bam_reader(bam, args.threads);
                let header = bam::Header::from_template(bam.header());
                let mut out = bam_writer(out, &bam, args.threads);
                // set up options
                let nuc_opts = fibertools_rs::cli::AddNucleosomeOptions {
                    bam: "-".to_string(),
                    out: "-".to_string(),
                    nucleosome_length: *nucleosome_length,
                    combined_nucleosome_length: *combined_nucleosome_length,
                    min_distance_added: *min_distance_added,
                    distance_from_end: *distance_from_end,
                    allowed_m6a_skips: *allowed_m6a_skips,
                    min_ml_score: 0,
                };
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
                    nuc_opts,
                );
                log::info!("{} reads included at once in batch prediction.", batch_size);
                predict_m6a::read_bam_into_fiberdata(&mut bam, &mut out, &predict_options);
            }
            #[cfg(not(feature = "predict"))]
            {
                log::error!("m6A predictions are not enabled with this install of fibertools.");
                log::error!("Please recompile with: `--all-features`.");
                log::error!("If fibertools was installed via bioconda you will not be able to do m6A predictions. Recent size limits on packages prevent us from including the pytorch libraries necessary for predictions within the bioconda install. Please install via Cargo following the directions here: https://fiberseq.github.io/fibertools-rs/INSTALL.html."
                )
            }
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
