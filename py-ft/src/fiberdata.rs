use fibertools_rs::fiber::FiberseqData;
use pyo3::iter::IterNextOutput;
use pyo3::prelude::*;
use rust_htslib::{bam, bam::ext::BamRecordExtensions, bam::record::Aux, bam::Read};
use std::time;
use std::vec::IntoIter;

#[pyclass]
/// Class for fiberseq data. This class corresponds to a single record in the bam file.
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
    /// offset for centering the fiber
    #[pyo3(get, set)]
    offset: Option<i64>,
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

impl Fiberdata {
    fn new(fiber: FiberseqData, offset: Option<i64>) -> Self {
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

        let aligned_block_pairs = fiber.record.aligned_block_pairs().collect();
        Self {
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
            offset,
        }
    }
}

#[pyclass]
/// Open a fiberseq bam file. Must have an index.
pub struct Fiberbam {
    bam: bam::IndexedReader,
    header: bam::Header,
    start: time::Instant,
}

/// pure rust functions for Fiberbam
impl Fiberbam {
    fn _fetch_helper(
        &mut self,
        chrom: &str,
        start: Option<i64>,
        end: Option<i64>,
    ) -> Vec<FiberseqData> {
        let fetch_args = match (start, end) {
            (Some(start), Some(end)) => bam::FetchDefinition::from((chrom, start, end)),
            _ => bam::FetchDefinition::from(chrom),
        };

        let head_view = bam::HeaderView::from_header(&self.header);
        self.bam.fetch(fetch_args).expect("unable to fetch region");
        let records: Vec<bam::Record> = self.bam.records().map(|r| r.unwrap()).collect();

        log::info!(
            "{} records fetched in {:.2}s",
            records.len(),
            self._time_from_last()
        );
        let fiberdata = FiberseqData::from_records(records, &head_view, 0);
        log::info!(
            "Fiberdata made for {} records in {:.2}s",
            fiberdata.len(),
            self._time_from_last()
        );
        fiberdata
    }
}

#[pymethods]
impl Fiberbam {
    #[new]
    #[pyo3(signature = (bam_file, threads = 8))]
    fn new(bam_file: &str, threads: usize) -> Self {
        let mut bam =
            bam::IndexedReader::from_path(bam_file).expect("unable to open indexed bam file");
        bam.set_threads(threads).unwrap();
        let header = bam::Header::from_template(bam.header());
        let start = time::Instant::now();
        Self { bam, header, start }
    }

    /// Returns an iterator over :class:`~pyft.Fiberdata` objects for the selected region.
    /// Arguments for the region to fetch can be provided in multiple ways, e.g.:
    ///
    ///     fiberdata.fetch("chr1", 100, 200)
    ///
    ///     fiberdata.fetch("chr1", start=100, end=200)
    ///
    ///     fiberdata.fetch("chr1:100-200")
    ///
    ///     fiberdata.fetch("chr1")
    #[pyo3(signature = (chrom, start = None, end=None))]
    pub fn fetch(&mut self, chrom: &str, start: Option<i64>, end: Option<i64>) -> Fiberiter {
        Fiberiter::build_fiberdata_iter(self._fetch_helper(chrom, start, end))
    }

    /// Returns an iterator over :class:`~pyft.Fiberdata` objects; however, the data is centered around the region that has been fetched (`should` work the same as **ft center**).
    #[pyo3(signature = (chrom, start, end, strand = '+'))]
    pub fn center(&mut self, chrom: &str, start: i64, end: i64, strand: char) -> Fiberiter {
        let position = if strand == '-' { end - 1 } else { start };
        let center_position = fibertools_rs::center::CenterPosition {
            chrom: chrom.to_string(),
            position,
            strand,
        };
        let fiberdata = self._fetch_helper(chrom, Some(position), Some(position + 1));
        let fibderdata: Vec<FiberseqData> = fiberdata
            .into_iter()
            .filter_map(|fiber| fiber.center(&center_position))
            .collect();
        Fiberiter::build_fiberdata_iter(fibderdata)
    }

    fn _time_from_last(&mut self) -> f64 {
        let elapsed = self.start.elapsed().as_secs_f64();
        self.start = time::Instant::now();
        elapsed
    }
}

#[pyclass]
/// An iterator over :class:`~pyft.Fiberdata` objects.
pub struct Fiberiter {
    fiberdata: IntoIter<FiberseqData>,
    length: usize,
}

/// rust only functions for Fiberiter
impl Fiberiter {
    fn build_fiberdata_iter(fiberdata: Vec<FiberseqData>) -> Self {
        let length = fiberdata.len();
        Self {
            fiberdata: fiberdata.into_iter(),
            length,
        }
    }
}

#[pymethods]
impl Fiberiter {
    fn __next__(&mut self) -> IterNextOutput<Fiberdata, &'static str> {
        let data = self.fiberdata.next();
        match data {
            Some(fiber) => IterNextOutput::Yield(Fiberdata::new(fiber, None)),
            None => IterNextOutput::Return("Ended"),
        }
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __len__(&self) -> usize {
        self.length
    }

    fn len(&self) -> usize {
        self.length
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
