use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use std::process::Command; // Run programs

#[test]
fn test_m6a() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("ft")?;
    cmd.arg("m6a")
        .arg("-v")
        .arg("tests/data/all.bam")
        .arg("/dev/null");
    cmd.assert()
        .success()
        .stderr(predicate::str::contains("done! Time elapsed:"));
    Ok(())
}

#[test]
fn test_center() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("ft")?;
    cmd.arg("center")
        .arg("-v")
        .arg("tests/data/center.bam")
        .arg("--bed")
        .arg("tests/data/center.bed");
    cmd.assert()
        .success()
        .stderr(predicate::str::contains("done! Time elapsed:"));
    Ok(())
}

#[test]
fn test_extract() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("ft")?;
    cmd.arg("extract")
        .arg("-v")
        .arg("tests/data/all.bam")
        .arg("--all")
        .arg("/dev/null");
    cmd.assert()
        .success()
        .stderr(predicate::str::contains("done! Time elapsed:"));
    Ok(())
}

#[test]
fn test_nucleosomes() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("ft")?;
    cmd.arg("add-nucleosomes")
        .arg("-v")
        .arg("tests/data/all.bam")
        .arg("/dev/null");
    cmd.assert()
        .success()
        .stderr(predicate::str::contains("done! Time elapsed:"));
    Ok(())
}

#[test]
fn test_clear() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("ft")?;
    cmd.arg("clear")
        .arg("-v")
        .arg("tests/data/all.bam")
        .arg("/dev/null");
    cmd.assert()
        .success()
        .stderr(predicate::str::contains("done! Time elapsed:"));
    Ok(())
}

#[test]
fn test_strip_basemods() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("ft")?;
    cmd.arg("strip-basemods")
        .arg("-v")
        .arg("tests/data/all.bam")
        .arg("/dev/null");
    cmd.assert()
        .success()
        .stderr(predicate::str::contains("done! Time elapsed:"));
    Ok(())
}

#[test]
fn test_footprint() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("ft")?;
    cmd.arg("footprint")
        .arg("-v")
        .arg("tests/data/ctcf.bam")
        .arg("--bed")
        .arg("tests/data/ctcf.bed.gz")
        .arg("--yaml")
        .arg("tests/data/ctcf.yaml")
        .arg("-o")
        .arg("/dev/null");
    cmd.assert()
        .success()
        .stderr(predicate::str::contains("done! Time elapsed:"));
    Ok(())
}

#[test]
fn test_pileup() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("ft")?;
    cmd.arg("pileup")
        .arg("-v")
        .arg("tests/data/all.bam")
        .arg("-o")
        .arg("/dev/null");
    cmd.assert()
        .success()
        .stderr(predicate::str::contains("done! Time elapsed:"));
    Ok(())
}

#[test]
fn test_qc() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("ft")?;
    cmd.arg("qc")
        .arg("-v")
        .arg("--acf")
        .arg("--m6a-per-msp")
        .arg("tests/data/all.bam")
        .arg("/dev/null");
    cmd.assert()
        .success()
        .stderr(predicate::str::contains("done! Time elapsed:"));
    Ok(())
}
