use crate::cli::ValidateOptions;
use anyhow::{Error, Result};
use std::io::IsTerminal;

/// print a message to stderr in bright red
/// but only if the output is a terminal
pub fn eprintln_red(message: &str) {
    if std::io::stderr().is_terminal() {
        eprintln!("\x1b[1;31m{message}\x1b[0m");
    } else {
        eprintln!("{message}");
    }
}

/// same but with bright green
pub fn eprintln_green(message: &str) {
    if std::io::stderr().is_terminal() {
        eprintln!("\x1b[1;32m{message}\x1b[0m");
    } else {
        eprintln!("{message}");
    }
}

pub fn validate_fiberseq_bam(opts: &mut ValidateOptions) -> Result<()> {
    let mut n_reads = 0;
    let mut n_m6a = 0;
    let mut n_fire_calls = 0;
    let mut n_nucleosomes = 0;
    let mut n_aligned = 0;
    let mut n_phased = 0;
    let mut n_kinetics = 0;

    let mut bam = opts.bam.bam_reader();

    for fiber in opts.bam.fibers(&mut bam) {
        let m6a_okay = !fiber.m6a.annotations.is_empty();
        let nuc_okay = !fiber.nuc.annotations.is_empty();
        n_fire_calls += fiber.msp.qual().iter().filter(|q| **q > 0).count();

        n_aligned += !fiber.record.is_unmapped() as usize;
        n_phased += fiber.record.aux(b"HP").is_ok() as usize;

        let mut seen_kinetics = 0;
        seen_kinetics += fiber.record.aux(b"fp").is_ok() as usize;
        seen_kinetics += fiber.record.aux(b"fi").is_ok() as usize;
        seen_kinetics += fiber.record.aux(b"rp").is_ok() as usize;
        seen_kinetics += fiber.record.aux(b"ri").is_ok() as usize;
        n_kinetics += (seen_kinetics == 4) as usize;

        if m6a_okay {
            n_m6a += 1;
        }
        if nuc_okay {
            n_nucleosomes += 1;
        }
        n_reads += 1;

        // break if we have reached the number of reads to validate
        if n_reads == opts.reads {
            break;
        }
    }

    // total reads
    eprintln!("Total reads tested: {n_reads}");

    // frac with m6a, frac with nucleosomes, and frac with fire
    eprintln!(
        "Fraction with m6A: {:.2}%\nFraction with nucleosomes: {:.2}%\nNumer of FIRE calls: {}",
        n_m6a as f64 / n_reads as f64 * 100.0,
        n_nucleosomes as f64 / n_reads as f64 * 100.0,
        n_fire_calls,
    );

    // alignment and phasing
    eprintln!(
        "Fraction aligned: {:.2}%\nFraction phased: {:.2}%\nFraction with kinetics: {:.2}%\n",
        n_aligned as f64 / n_reads as f64 * 100.0,
        n_phased as f64 / n_reads as f64 * 100.0,
        n_kinetics as f64 / n_reads as f64 * 100.0,
    );

    let passes_m6a = n_m6a as f64 / n_reads as f64 >= opts.m6a;
    let passes_nucleosomes = n_nucleosomes as f64 / n_reads as f64 >= opts.nuc;
    let passes_fire = n_fire_calls > 0 || !opts.fire;
    let passes_phased = n_phased as f64 / n_reads as f64 >= opts.phased;
    let passes_aligned = n_aligned as f64 / n_reads as f64 >= opts.aligned;
    let passes_kinetics = n_kinetics as f64 / n_reads as f64 >= opts.kinetics;

    if passes_m6a
        && passes_nucleosomes
        && passes_fire
        && passes_phased
        && passes_aligned
        && passes_kinetics
    {
        eprintln_green("Fiber-seq BAM file is valid");
    } else {
        eprintln_red("Fiber-seq BAM file is invalid");

        if !passes_m6a {
            eprintln!(
                "\t- No m6A calls found, please check the m6A calling was performed correctly",
            );
        }
        if !passes_nucleosomes {
            eprintln!("\t- No nucleosome calls found, please check the nucleosome calling was performed correctly. Nucleosome calls can be added with `ft add-nucleosomes`");
        }
        if !passes_fire {
            eprintln!("\t- No FIRE calls found, please check FIRE calling was performed. FIRE calls can be added with `ft fire`");
        }
        if !passes_aligned {
            eprintln!(
                "\t- Less than {}% of reads are aligned, please check the alignment was performed correctly",
                opts.aligned * 100.0
            );
        }
        if !passes_phased {
            eprintln!(
                "\t- Less than {}% of reads are phased, please check the phasing was performed correctly",
                opts.phased * 100.0
            );
        }
        if !passes_kinetics {
            eprintln!(
                "\t- Less than {}% of reads have kinetics information, please check the kinetics calling was performed correctly",
                opts.kinetics * 100.0
            );
        }
        eprintln!();
        return Err(Error::msg("Fiber-seq BAM file is invalid"));
    }
    Ok(())
}
