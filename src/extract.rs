use super::fiber::FiberseqData;
use super::*;
use crate::cli::ExtractOptions;
use rayon::prelude::*;

pub struct FiberOut {
    pub m6a: Option<Box<dyn Write>>,
    pub cpg: Option<Box<dyn Write>>,
    pub msp: Option<Box<dyn Write>>,
    pub nuc: Option<Box<dyn Write>>,
    pub all: Option<Box<dyn Write>>,
    pub reference: bool,
    pub simplify: bool,
    pub quality: bool,
    pub min_ml_score: u8,
}

impl FiberOut {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        m6a: &Option<String>,
        cpg: &Option<String>,
        msp: &Option<String>,
        nuc: &Option<String>,
        all: &Option<String>,
        reference: bool,
        simplify: bool,
        quality: bool,
        min_ml_score: u8,
    ) -> Result<Self> {
        let m6a = match m6a {
            Some(m6a) => Some(writer(m6a)?),
            None => None,
        };
        let cpg = match cpg {
            Some(cpg) => Some(writer(cpg)?),
            None => None,
        };
        let msp = match msp {
            Some(msp) => Some(writer(msp)?),
            None => None,
        };
        let nuc = match nuc {
            Some(nuc) => Some(writer(nuc)?),
            None => None,
        };
        let all = match all {
            Some(all) => Some(writer(all)?),
            None => None,
        };

        Ok(FiberOut {
            m6a,
            cpg,
            msp,
            nuc,
            all,
            reference,
            simplify,
            quality,
            min_ml_score,
        })
    }
}

pub fn process_bam_chunk(fiber_data: Vec<FiberseqData>, out_files: &mut FiberOut) {
    match &mut out_files.m6a {
        Some(m6a) => {
            let out: Vec<String> = fiber_data
                .par_iter()
                .map(|r| r.write_m6a(out_files.reference))
                .collect();
            for line in out {
                write_to_file(&line, m6a);
            }
        }
        None => {}
    }
    match &mut out_files.cpg {
        Some(cpg) => {
            let out: Vec<String> = fiber_data
                .par_iter()
                .map(|r| r.write_cpg(out_files.reference))
                .collect();
            for line in out {
                write_to_file(&line, cpg);
            }
        }
        None => {}
    }
    match &mut out_files.msp {
        Some(msp) => {
            let out: Vec<String> = fiber_data
                .par_iter()
                .map(|r| r.write_msp(out_files.reference))
                .collect();
            for line in out {
                write_to_file(&line, msp);
            }
        }
        None => {}
    }
    match &mut out_files.nuc {
        Some(nuc) => {
            let out: Vec<String> = fiber_data
                .par_iter()
                .map(|r| r.write_nuc(out_files.reference))
                .collect();
            for line in out {
                write_to_file(&line, nuc);
            }
        }
        None => {}
    }
    match &mut out_files.all {
        Some(all) => {
            let out: Vec<String> = fiber_data
                .par_iter()
                .map(|r| r.write_all(out_files.simplify, out_files.quality))
                .collect();
            for line in out {
                write_to_file(&line, all);
            }
        }
        None => {}
    }
}

pub fn extract_contained(extract_opts: &mut ExtractOptions) {
    let mut bam = extract_opts.input.bam_reader();
    let mut out_files = FiberOut::new(
        &extract_opts.m6a,
        &extract_opts.cpg,
        &extract_opts.msp,
        &extract_opts.nuc,
        &extract_opts.all,
        extract_opts.reference,
        extract_opts.simplify,
        extract_opts.quality,
        extract_opts.input.min_ml_score,
    )
    .unwrap();

    // print the header if in all mode
    match &mut out_files.all {
        Some(all) => {
            write!(
                all,
                "{}",
                FiberseqData::all_header(out_files.simplify, out_files.quality)
            )
            .unwrap();
        }
        None => {}
    }

    // process bam in chunks
    // keeps mem pretty low, about 1GB per thread
    let fibers = extract_opts.input.fibers(&mut bam);
    for chunk in &fibers.chunks(1_000) {
        process_bam_chunk(chunk.collect(), &mut out_files);
    }
}
