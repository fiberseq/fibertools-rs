use std::error::Error;
use vergen_git2::*;

fn vergen() -> Result<(), Box<dyn Error>> {
    // NOTE: This will output everything, and requires all features enabled.
    // NOTE: See the specific builder documentation for configuration options.
    let build = BuildBuilder::all_build()?;
    let git2 = Git2Builder::all_git()?;
    let cargo = CargoBuilder::all_cargo()?;

    let status = Emitter::default()
        .add_instructions(&build)?
        .add_instructions(&git2)?
        .add_instructions(&cargo)?
        .fail_on_error()
        .emit();

    // set the env variable myself if the status failed
    if status.is_err() {
        eprintln!("Failed to get version information from git");
        eprintln!("Likely building a published version from cargo or bioconda.");
        std::env::set_var("VERGEN_GIT_DESCRIBE", "unknown");
    }

    Ok(())
}

fn main() {
    /*
    use std::process::Command;
    let output = Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .unwrap();
    let git_hash = String::from_utf8(output.stdout).unwrap();
    let version = env!("CARGO_PKG_VERSION");
    println!(
        "cargo:rustc-env=CARGO_GIT_HASH={}
        cargo:rustc-env=CARGO_LONG_VERSION={} commit:{}",
        git_hash, version, git_hash
    );
    */
    vergen().expect("Unable to set version and git hash");

    // Generate the model code and state file from the ONNX file.
    use burn_import::onnx::ModelGen;
    use burn_import::onnx::RecordType;
    for x in &[
        "src/m6a_burn/two_zero.onnx",
        "src/m6a_burn/two_two.onnx",
        "src/m6a_burn/three_two.onnx",
        "src/m6a_burn/revio.onnx",
    ] {
        ModelGen::new()
            .input(x) // Path to the ONNX model
            .out_dir("m6a_burn/") // Directory for the generated Rust source file (under target/)
            .record_type(RecordType::Bincode)
            .embed_states(true)
            .run_from_script();
    }
}
