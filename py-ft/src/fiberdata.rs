use fibertools_rs::extract::FiberseqData;
use pyo3::iter::IterNextOutput;
use pyo3::prelude::*;
use rust_htslib::{bam, bam::ext::BamRecordExtensions, bam::record::Aux, bam::Read};
use std::vec::IntoIter;

#[pyclass]
/// Record class for fiberseq data, corresponding to a single bam record
pub struct Fiberdata {
    /// Number of ccs passes
    #[pyo3(get, set)]
    pub ec: i64,
    /// Name of the read
    #[pyo3(get, set)]
    pub qname: String,
    /// SAM flag
    #[pyo3(get, set)]
    pub sam_flag: u16,
    /// Read quality
    #[pyo3(get, set)]
    pub rq: String,
    /// Haplotype tag
    #[pyo3(get, set)]
    pub hp: String,
    /// Read sequence
    #[pyo3(get, set)]
    pub seq: String,
    /// Chromosome
    #[pyo3(get, set)]
    pub chrom: String,
    /// Chromosome start
    #[pyo3(get, set)]
    pub start: i64,
    /// Chromosome end
    #[pyo3(get, set)]
    pub end: i64,
    /// Strand
    #[pyo3(get, set)]
    pub strand: char,
    #[pyo3(get, set)]
    /// Read group
    pub rg: String,
    /// m6a :class:`~pyft.Basemods` object
    #[pyo3(get, set)]
    pub m6a: Basemods,
    /// cpg :class:`~pyft.Basemods` object
    #[pyo3(get, set)]
    pub cpg: Basemods,
    /// :class:`~pyft.Ranges` object for msp features
    #[pyo3(get, set)]
    pub msp: Ranges,
    /// :class:`~pyft.Ranges` object for nuc features
    #[pyo3(get, set)]
    pub nuc: Ranges,
    /// Aligned block pairs from the bam record
    #[pyo3(get)]
    pub aligned_block_pairs: Vec<([i64; 2], [i64; 2])>,
}
#[pymethods]
impl Fiberdata {
    pub fn __str__(&self) -> String {
        format!(
            "fiber: {}\t\
            chrom: {}\tstart: {}\tend {}\t\
            num m6a: {}\t num cpg: {}\t\
            num nuc: {}\t num msp: {}",
            self.qname,
            self.chrom,
            self.start,
            self.end,
            self.m6a.starts.len(),
            self.cpg.starts.len(),
            self.nuc.starts.len(),
            self.msp.starts.len()
        )
    }

    /// return the length of the sequence
    pub fn get_seq_length(&self) -> usize {
        self.seq.len()
    }

    /// liftover query (fiber) positions to the reference
    pub fn lift_query_positions(&self, positions: Vec<i64>) -> Vec<i64> {
        bamlift::lift_query_positions(&self.aligned_block_pairs, &positions)
    }

    /// liftover reference positions to the query (fiber)
    pub fn lift_reference_positions(&self, positions: Vec<i64>) -> Vec<i64> {
        bamlift::lift_reference_positions(&self.aligned_block_pairs, &positions)
    }

    /// Return a :class:`~pyft.Ranges` object for the continuous aligned regions of the fiber
    pub fn get_aligned_blocks_as_ranges(&self) -> Ranges {
        let (fiber, genome): (Vec<[i64; 2]>, Vec<[i64; 2]>) =
            self.aligned_block_pairs.iter().cloned().unzip();
        let starts = fiber.iter().map(|x| x[0]).collect();
        let lengths = fiber.iter().map(|x| x[1] - x[0]).collect();
        let reference_starts = genome.iter().map(|x| x[0]).collect();
        let reference_lengths = genome.iter().map(|x| x[1] - x[0]).collect();
        Ranges::new(starts, lengths, reference_starts, reference_lengths)
    }
}

fn new_py_fiberdata(fiber: FiberseqData) -> Fiberdata {
    // PB features
    let ec = fiber.ec.round() as i64;

    // record features
    let qname = std::str::from_utf8(fiber.record.qname()).unwrap();
    let sam_flag = fiber.record.flags();
    let rq = match fiber.get_rq() {
        Some(x) => format!("{}", x),
        None => ".".to_string(),
    };
    let hp = fiber.get_hp();
    let seq = String::from_utf8_lossy(&fiber.record.seq().as_bytes()).to_string();
    let (chrom, start, end, strand): (&str, i64, i64, char) = if fiber.record.is_unmapped() {
        (".", 0, 0, '.')
    } else {
        let strand = if fiber.record.is_reverse() { '-' } else { '+' };
        (
            &fiber.target_name,
            fiber.record.reference_start(),
            fiber.record.reference_end(),
            strand,
        )
    };
    let rg = if let Ok(Aux::String(f)) = fiber.record.aux(b"RG") {
        log::trace!("{f}");
        f
    } else {
        "."
    };

    // fiberseq features
    let m6a = fiber.base_mods.m6a();
    let m6a = Basemods::new(m6a.0, m6a.1, m6a.2);
    let cpg = fiber.base_mods.cpg();
    let cpg = Basemods::new(cpg.0, cpg.1, cpg.2);
    let msp_starts = fiber.get_msp(false, true);
    let msp_lengths = fiber.get_msp(false, false);
    let ref_msp_starts = fiber.get_msp(true, true);
    let ref_msp_lengths = fiber.get_msp(true, false);
    let msp = Ranges::new(msp_starts, msp_lengths, ref_msp_starts, ref_msp_lengths);

    let nuc_starts = fiber.get_nuc(false, true);
    let nuc_lengths = fiber.get_nuc(false, false);
    let ref_nuc_starts = fiber.get_nuc(true, true);
    let ref_nuc_lengths = fiber.get_nuc(true, false);
    let nuc = Ranges::new(nuc_starts, nuc_lengths, ref_nuc_starts, ref_nuc_lengths);

    /*
    let nuc = Ranges {
        starts: fiber.nuc.starts,
        ends: fiber.nuc.ends,
        lengths: fiber.nuc.lengths,
        reference_starts: fiber.nuc.reference_starts,
        reference_ends: fiber.nuc.reference_ends,
        reference_lengths: fiber.nuc.reference_lengths,
    };
    let msp = Ranges {
        starts: fiber.msp.starts,
        ends: fiber.msp.ends,
        lengths: fiber.msp.lengths,
        reference_starts: fiber.msp.reference_starts,
        reference_ends: fiber.msp.reference_ends,
        reference_lengths: fiber.msp.reference_lengths,
    };
    */
    let aligned_block_pairs = fiber.record.aligned_block_pairs().collect();
    Fiberdata {
        ec,
        qname: qname.to_string(),
        sam_flag,
        rq,
        hp,
        seq,
        chrom: chrom.to_string(),
        start,
        end,
        strand,
        rg: rg.to_string(),
        m6a,
        cpg,
        msp,
        nuc,
        aligned_block_pairs,
    }
}

#[pyclass]
/// Create a fibertools iterator from an indexed bam file.
/// Must provide a valid chrom, start, and end.
/// Returns an iterator over :class:`~pyft.Fiberdata` objects.
pub struct FiberdataFetch {
    fiberdata: IntoIter<FiberseqData>,
    //_inner: Rc<RefCell<bam::IndexedReader>>,
    //fiberdata: Box<dyn Iterator<Item = FiberseqData> + Send + 'static>,
    //fiberdata: Box<dyn Iterator<Item = FiberseqData> + Send + 'static>,
}

#[pymethods]
impl FiberdataFetch {
    #[new]
    pub fn new(f: &str, chrom: &str, start: i64, end: i64) -> Self {
        let mut bam = bam::IndexedReader::from_path(f).expect("unable to open indexed bam file");
        bam.set_threads(8).unwrap();
        let header = bam::Header::from_template(bam.header());
        let head_view = bam::HeaderView::from_header(&header);
        bam.fetch((chrom, start, end))
            .expect("unable to fetch region");
        let records: Vec<bam::Record> = bam.records().map(|r| r.unwrap()).collect();
        log::info!("{} records fetched", records.len());
        let fiberdata = FiberseqData::from_records(records, &head_view, 0).into_iter();
        log::info!("fiberdata created from records");
        /*let fiberdata = bam
        .records()
        .map(|rec| FiberseqData::new(&rec.unwrap(), None, 0))
        .into_iter();*/
        Self { fiberdata }
    }

    fn __next__(&mut self) -> IterNextOutput<Fiberdata, &'static str> {
        let data = self.fiberdata.next();
        match data {
            Some(fiber) => IterNextOutput::Yield(new_py_fiberdata(fiber)),
            None => {
                log::info!("\n\nDone iterating over fibers");
                IterNextOutput::Return("Ended")
            }
        }
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
}

/// Class for describing base modifications within fiberseq data.
/// For example, m6a and 5mC (cpg) features.
#[pyclass]
#[derive(Clone)]
pub struct Basemods {
    /// Basemod starts
    #[pyo3(get, set)]
    pub starts: Vec<i64>,
    /// Basemod reference starts
    #[pyo3(get, set)]
    pub reference_starts: Vec<i64>,
    /// Basemod ML
    #[pyo3(get, set)]
    pub ml: Vec<u8>,
}

#[pymethods]
impl Basemods {
    #[new]
    pub fn new(starts: Vec<i64>, reference_starts: Vec<i64>, ml: Vec<u8>) -> Self {
        Self {
            starts,
            reference_starts,
            ml,
        }
    }
}

/// Class for describing ranges within fiberseq data.
/// For example, nucleosomes and msp features.
#[pyclass]
#[derive(Clone)]
pub struct Ranges {
    /// Range starts
    #[pyo3(get, set)]
    pub starts: Vec<i64>,
    /// Range ends
    #[pyo3(get, set)]
    pub ends: Vec<i64>,
    /// Range lengths
    #[pyo3(get, set)]
    pub lengths: Vec<i64>,
    /// Reference starts
    #[pyo3(get, set)]
    pub reference_starts: Vec<i64>,
    /// Reference ends
    #[pyo3(get, set)]
    pub reference_ends: Vec<i64>,
    /// Reference lengths
    #[pyo3(get, set)]
    pub reference_lengths: Vec<i64>,
}

#[pymethods]
impl Ranges {
    #[new]
    pub fn new(
        starts: Vec<i64>,
        lengths: Vec<i64>,
        reference_starts: Vec<i64>,
        reference_lengths: Vec<i64>,
    ) -> Self {
        let ends = starts
            .iter()
            .zip(lengths.iter())
            .map(|(s, l)| s + l)
            .collect();
        let reference_ends = reference_starts
            .iter()
            .zip(reference_lengths.iter())
            .map(|(s, l)| s + l)
            .collect();
        Self {
            starts,
            ends,
            lengths,
            reference_starts,
            reference_ends,
            reference_lengths,
        }
    }
}
