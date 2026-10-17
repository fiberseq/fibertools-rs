use crate::cli::PgLiftOptions;
use anyhow::{Context, Result};
use duct::cmd;

pub fn run_pg_lift(opts: &PgLiftOptions) -> Result<()> {
    let input_type = if opts.bam { "BAM" } else { "BED" };

    log::info!("Starting pangenome annotation lift pipeline");
    log::info!("  Source reference: {}", opts.reference);
    log::info!("  Input file: {} ({})", opts.input, input_type);
    log::info!("  Pangenome graph: {}", opts.graph);
    log::info!("  Target sample: {}", opts.target);
    log::info!("  panSN prefix: {}", opts.prefix);

    // Pre-validate all requirements before starting pipeline
    validate_pipeline_requirements(opts)?;

    // Create the complete pipeline using pipes (like the original shell script)
    run_piped_pipeline(opts)?;

    log::info!("Pangenome annotation lift completed successfully");
    log::info!("Output written to: {}", opts.out);

    Ok(())
}

/// Validate all pipeline requirements before execution
fn validate_pipeline_requirements(opts: &PgLiftOptions) -> Result<()> {
    log::info!("Validating pipeline requirements...");

    // 1. Check vg binary exists and is functional
    which_vg(&opts.vg_binary).context("vg binary validation failed")?;

    // 2. Check input files exist and are readable
    if !std::path::Path::new(&opts.reference).exists() {
        return Err(anyhow::anyhow!(
            "Reference FASTA file not found: {}",
            opts.reference
        ));
    }

    if !std::path::Path::new(&opts.input).exists() {
        return Err(anyhow::anyhow!("Input file not found: {}", opts.input));
    }

    if !std::path::Path::new(&opts.graph).exists() {
        return Err(anyhow::anyhow!(
            "Pangenome graph file not found: {}",
            opts.graph
        ));
    }

    log::info!("✓ All pipeline requirements validated");
    Ok(())
}

/// Run the complete pipeline using pipes between commands for efficiency
fn run_piped_pipeline(opts: &PgLiftOptions) -> Result<()> {
    // Build command arguments as strings to avoid lifetime issues
    let vg_threads_str = opts.vg_threads.to_string();

    // Build verbose/quiet arguments to pass to subcommands
    let mut verbose_args = Vec::new();
    if opts.global.quiet {
        verbose_args.push("--quiet".to_string());
    }
    for _ in 0..opts.global.verbose {
        verbose_args.push("--verbose".to_string());
    }

    // Create temporary file for header transfer between first and final commands
    let temp_header_file =
        tempfile::NamedTempFile::new().context("Failed to create temporary header file")?;
    let temp_header_path = temp_header_file.path().to_string_lossy().to_string();
    log::info!("Using temporary header file: {}", temp_header_path);

    // Create the first command based on input type
    if opts.bam {
        log::info!(
            "Running BAM pipeline: pg-pansn --prefix | vg inject | vg surject | pg-pansn --strip"
        );
    } else {
        log::info!(
            "Running BED pipeline: pg-inject | vg inject | vg surject | pg-inject --extract"
        );
    }
    let first_cmd = cmd(
        std::env::current_exe()?,
        first_stage_args(opts, &temp_header_path, &verbose_args),
    );

    // vg commands are the same for both pipelines
    let vg_inject_cmd = cmd!(
        &opts.vg_binary,
        "inject",
        "-t",
        &vg_threads_str,
        "-x",
        &opts.graph,
        "-" // Read from stdin
    );

    let vg_surject_cmd = cmd!(
        &opts.vg_binary,
        "surject",
        "-C",
        "0",
        "-t",
        &vg_threads_str,
        "-x",
        &opts.graph,
        "-r",
        "-b",
        "-n",
        &opts.target,
        "-" // Read from stdin
    );

    // Final command depends on input type - always use uncompressed for final BAM output
    let final_cmd = cmd(
        std::env::current_exe()?,
        final_stage_args(opts, &temp_header_path, &verbose_args),
    );

    // Execute the complete pipeline
    let pipeline = first_cmd
        .pipe(vg_inject_cmd)
        .pipe(vg_surject_cmd)
        .pipe(final_cmd);

    pipeline.run().context("Failed to execute pipeline")?;

    Ok(())
}

/// Arguments of the first `ft` stage: `pg-pansn --prefix` for BAM input,
/// `pg-inject` (mock fibertig BAM) for BED input. Intermediate output is
/// always uncompressed and the header is written to `header_path`.
fn first_stage_args(opts: &PgLiftOptions, header_path: &str, verbose: &[String]) -> Vec<String> {
    let mut args: Vec<String> = if opts.bam {
        vec![
            "pg-pansn".into(),
            opts.input.clone(),
            "--prefix".into(),
            opts.prefix.clone(),
            "--uncompressed".into(),
            "--header-out".into(),
            header_path.into(),
        ]
    } else {
        vec![
            "pg-inject".into(),
            opts.reference.clone(),
            "--prefix".into(),
            opts.prefix.clone(),
            "--split-size".into(),
            opts.split_size.to_string(),
            "--bed".into(),
            opts.input.clone(),
            "--uncompressed".into(),
            "--header-out".into(),
            header_path.into(),
        ]
    };
    args.extend(verbose.iter().cloned());
    args
}

/// Arguments of the final `ft` stage: `pg-pansn --strip` for BAM input (its
/// output BAM is positional), `pg-inject --extract` for BED input.
fn final_stage_args(opts: &PgLiftOptions, header_path: &str, verbose: &[String]) -> Vec<String> {
    let delimiter = opts.delimiter.to_string();
    let mut args: Vec<String> = if opts.bam {
        vec![
            "pg-pansn".into(),
            "-".into(), // Read from stdin
            opts.out.clone(),
            "--strip".into(),
            "--delimiter".into(),
            delimiter,
            "--copy-header".into(),
            header_path.into(),
            "--uncompressed".into(),
        ]
    } else {
        vec![
            "pg-inject".into(),
            "-".into(), // Read from stdin
            "--extract".into(),
            "--strip".into(),
            "--delimiter".into(),
            delimiter,
            "--out".into(),
            opts.out.clone(),
            "--copy-header".into(),
            header_path.into(),
        ]
    };
    args.extend(verbose.iter().cloned());
    args
}

/// Check if vg binary exists and is executable
fn which_vg(vg_binary: &str) -> Result<()> {
    cmd!(vg_binary, "help")
        .stdout_null()
        .stderr_null()
        .run()
        .with_context(|| {
            format!(
                "Failed to find vg binary at '{}'. Please install vg or specify correct path with --vg-binary",
                vg_binary
            )
        })?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{Cli, Commands};
    use clap::Parser;

    fn lift_opts(bam: bool) -> PgLiftOptions {
        let mut argv = vec![
            "ft",
            "pg-lift",
            "ref.fa",
            "-i",
            "in",
            "-g",
            "graph.gbz",
            "-o",
            "out",
            "-p",
            "HG002#1#",
            "--target",
            "GRCh38",
        ];
        if bam {
            argv.push("--bam");
        }
        match Cli::try_parse_from(argv).unwrap().command {
            Some(Commands::PgLift(opts)) => opts,
            _ => panic!("expected pg-lift"),
        }
    }

    /// Every `ft` stage pg-lift spawns must parse with ft's own CLI.
    #[test]
    fn pipeline_stage_args_parse() {
        let verbose = vec!["--verbose".to_string()];
        for bam in [true, false] {
            let opts = lift_opts(bam);
            for args in [
                first_stage_args(&opts, "header.sam", &verbose),
                final_stage_args(&opts, "header.sam", &verbose),
            ] {
                let argv = std::iter::once("ft".to_string()).chain(args.clone());
                if let Err(e) = Cli::try_parse_from(argv) {
                    panic!("ft {} does not parse: {e}", args.join(" "));
                }
            }
        }
    }
}
