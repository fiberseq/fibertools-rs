use crate::cli::ValidateOptions;
use anyhow::{Error, Result};
use std::io::IsTerminal;

/// print a message to stderr in bright red
/// but only if the output is a terminal
pub fn eprintln_red(message: &str) {
    if std::io::stderr().is_terminal() {
        eprintln!("\x1b[1;31m{}\x1b[0m", message);
    } else {
        eprintln!("{}", message);
    }
}

/// same but with bright green
pub fn eprintln_green(message: &str) {
    if std::io::stderr().is_terminal() {
        eprintln!("\x1b[1;32m{}\x1b[0m", message);
    } else {
        eprintln!("{}", message);
    }
}

pub fn validate_fiberseq_bam(opts: &mut ValidateOptions) -> Result<()> {
    let mut n_reads = 0;
    let mut n_valid = 0;
    let mut n_m6a = 0;
    let mut n_fire_calls = 0;
    let mut n_nucleosomes = 0;

    let mut bam = opts.bam.bam_reader();

    for fiber in opts.bam.fibers(&mut bam) {
        let m6a_okay = !fiber.m6a.starts.is_empty();
        let nuc_okay = !fiber.nuc.starts.is_empty();
        n_fire_calls += fiber.msp.qual.iter().filter(|q| **q > 0).count();

        if m6a_okay {
            n_m6a += 1;
        }
        if nuc_okay {
            n_nucleosomes += 1;
        }
        if m6a_okay && nuc_okay {
            n_valid += 1;
        }
        n_reads += 1;
    }
    // frac with m6a, frac with nucleosomes, and frac with fire
    eprintln!(
        "Fraction with m6A: {:.2}%\nFraction with nucleosomes: {:.2}%\nNumer of FIRE calls: {}\n",
        n_m6a as f64 / n_reads as f64 * 100.0,
        n_nucleosomes as f64 / n_reads as f64 * 100.0,
        n_fire_calls,
    );

    // total reads, total valid reads, and percent valid reads
    let frac_valid = n_valid as f64 / n_reads as f64 * 100.0;
    eprintln!(
        "Total reads tested: {}\nValid reads: {}\nFraction valid: {:.2}%\nMinimum validation rate to pass {:.2}%\n",
        n_reads, n_valid, frac_valid, 100.0*opts.min_valid_fraction,
    );
    let passes_fire = n_fire_calls > 0 || !opts.check_fire;

    if frac_valid >= opts.min_valid_fraction * 100.0 && passes_fire {
        eprintln_green("Fiber-seq BAM file is valid");
    } else {
        eprintln_red("Fiber-seq BAM file is invalid");

        if n_m6a == 0 {
            eprintln!(
                "\t- No m6A calls found, please check the m6A calling was performed correctly",
            );
        }
        if n_nucleosomes == 0 {
            eprintln!("\t- No nucleosome calls found, please check the nucleosome calling was performed correctly. Nucleosome calls can be added with `ft add-nucleosomes`");
        }
        if n_fire_calls == 0 && opts.check_fire {
            eprintln!("\t- No FIRE calls found, please check FIRE calling was performed. FIRE calls can be added with `ft fire`");
        }
        eprintln!();
        return Err(Error::msg("Fiber-seq BAM file is invalid"));
    }
    Ok(())
}
