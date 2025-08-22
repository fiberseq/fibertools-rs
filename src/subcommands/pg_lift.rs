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

    // 3. Check output directory is writable (if not stdout)
    if opts.out != "-" {
        if let Some(parent) = std::path::Path::new(&opts.out).parent() {
            if !parent.exists() {
                return Err(anyhow::anyhow!(
                    "Output directory does not exist: {}",
                    parent.display()
                ));
            }
        }
    }

    log::info!("✓ All pipeline requirements validated");
    Ok(())
}

/// Run the complete pipeline using pipes between commands for efficiency
fn run_piped_pipeline(opts: &PgLiftOptions) -> Result<()> {
    // Build command arguments as strings to avoid lifetime issues
    let split_size_str = opts.split_size.to_string();
    let vg_threads_str = opts.vg_threads.to_string();
    let delimiter_str = opts.delimiter.to_string();

    // Create the first command based on input type
    let first_cmd = if opts.bam {
        log::info!("Running BAM pipeline: pg-pansn --prefix | vg inject | vg surject | pg-inject --extract");
        cmd!(
            std::env::current_exe()?,
            "pg-pansn",
            &opts.input,
            "--prefix",
            &opts.prefix,
        )
    } else {
        log::info!(
            "Running BED pipeline: pg-inject | vg inject | vg surject | pg-inject --extract"
        );
        let mut pg_inject_args = vec![
            "pg-inject",
            &opts.reference,
            "--prefix",
            &opts.prefix,
            "--split-size",
            &split_size_str,
            "--bed",
            &opts.input,
        ];

        if opts.uncompressed {
            pg_inject_args.push("--uncompressed");
        }

        cmd(std::env::current_exe()?, pg_inject_args)
    };

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

    // Last command is always pg-inject --extract (works for both BAM and BED output)
    let pg_extract_cmd = cmd!(
        std::env::current_exe()?,
        "pg-inject",
        "-", // Read from stdin
        "--extract",
        "--strip",
        "--delimiter",
        &delimiter_str,
        "--out",
        &opts.out,
    );

    // Execute the complete pipeline
    let pipeline = first_cmd
        .pipe(vg_inject_cmd)
        .pipe(vg_surject_cmd)
        .pipe(pg_extract_cmd);

    pipeline.run().context("Failed to execute pipeline")?;

    Ok(())
}

/// Check if vg binary exists and is executable
fn which_vg(vg_binary: &str) -> Result<()> {
    cmd!(vg_binary, "--version")
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
