use fibertools_rs::utils::bamlift;
//use fibertools_rs::bio_io;
use fibertools_rs::fiber::FiberseqData;
use fibertools_rs::utils::input_bam::FiberFilters;
use pyo3::iter::IterNextOutput;
use pyo3::prelude::*;
use rust_htslib::{bam, bam::ext::BamRecordExtensions, bam::record::Aux, bam::Read};
use std::collections::HashMap;
use std::time;
use std::vec::IntoIter;
#[pyclass]
/// Class for fiberseq data. This class corresponds to a single record in the bam file.
pub struct Fiberdata {
    /// Number of c_s
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
    /// record
    fiber: FiberseqData,
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
    pub fn lift_query_positions(&self, positions: Vec<i64>) -> Vec<Option<i64>> {
        bamlift::lift_query_positions(&self.aligned_block_pairs, &positions)
    }

    /// liftover reference positions to the query (fiber)
    pub fn lift_reference_positions(&self, positions: Vec<i64>) -> Vec<Option<i64>> {
        bamlift::lift_reference_positions(&self.aligned_block_pairs, &positions)
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
        let rg = if let std::result::Result::Ok(Aux::String(f)) = fiber.record.aux(b"RG") {
            log::trace!("{f}");
            f
        } else {
            "."
        };

        // fiberseq features
        let m6a = Basemods::new(
            fiber.m6a.starts.clone(),
            fiber.m6a.reference_starts.clone(),
            fiber.m6a.qual.clone(),
        );
        let cpg = Basemods::new(
            fiber.cpg.starts.clone(),
            fiber.cpg.reference_starts.clone(),
            fiber.cpg.qual.clone(),
        );

        let nuc = Ranges {
            starts: fiber.nuc.starts.clone(),
            ends: fiber.nuc.ends.clone(),
            lengths: fiber.nuc.lengths.clone(),
            qual: fiber.nuc.qual.clone(),
            reference_starts: fiber.nuc.reference_starts.clone(),
            reference_ends: fiber.nuc.reference_ends.clone(),
            reference_lengths: fiber.nuc.reference_lengths.clone(),
        };
        let msp = Ranges {
            starts: fiber.msp.starts.clone(),
            ends: fiber.msp.ends.clone(),
            lengths: fiber.msp.lengths.clone(),
            qual: fiber.msp.qual.clone(),
            reference_starts: fiber.msp.reference_starts.clone(),
            reference_ends: fiber.msp.reference_ends.clone(),
            reference_lengths: fiber.msp.reference_lengths.clone(),
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
            fiber,
        }
    }
}

#[pyclass]
pub struct Fiberwriter {
    writer: bam::Writer,
}

#[pymethods]
impl Fiberwriter {
    #[new]
    pub fn new(output_bam_path: &str, template_bam_path: &str) -> Self {
        let header = bam::Header::from_template(
            bam::Reader::from_path(template_bam_path)
                .expect("unable to open template bam file")
                .header(),
        );
        let writer = bam::Writer::from_path(output_bam_path, &header, bam::Format::Bam)
            .expect("unable to open output bam file");
        Self { writer }
    }
    /// Write a Fiberdata object to a bam file.
    #[pyo3(signature = (fiberdata))]
    pub fn write(&mut self, fiberdata: &Fiberdata) -> PyResult<()> {
        self.writer
            .write(&fiberdata.fiber.record)
            .expect("unable to write record");
        pyo3::prelude::PyResult::Ok(())
    }
}

#[pyclass]
/// Open a fiberseq bam file. Must have an index.
pub struct Fiberbam {
    bam: bam::IndexedReader,
    reader: bam::Reader,
    header: bam::Header,
    target_dict: HashMap<i32, String>,
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
        let fiberdata = FiberseqData::from_records(records, &head_view, &FiberFilters::default());
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

        let mut reader = bam::Reader::from_path(bam_file).expect("unable to open bam file");
        reader.set_threads(threads).unwrap();

        let header = bam::Header::from_template(bam.header());
        let head_view = bam::HeaderView::from_header(&header);
        let target_dict = FiberseqData::dict_from_head_view(&head_view);
        let start = time::Instant::now();
        Self {
            bam,
            reader,
            header,
            target_dict,
            start,
        }
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
        let center_position = fibertools_rs::subcommands::center::CenterPosition {
            chrom: chrom.to_string(),
            position,
            strand,
        };
        let fiberdata = self._fetch_helper(chrom, Some(position), Some(position + 1));
        let fiberdata: Vec<FiberseqData> = fiberdata
            .into_iter()
            .filter_map(|fiber| fiber.center(&center_position))
            .collect();
        log::info!(
            "Fiberdata centered for {} records in {:.2}s",
            fiberdata.len(),
            self._time_from_last()
        );
        Fiberiter::build_fiberdata_iter(fiberdata)
    }

    fn _time_from_last(&mut self) -> f64 {
        let elapsed = self.start.elapsed().as_secs_f64();
        self.start = time::Instant::now();
        elapsed
    }

    fn __next__(&mut self) -> IterNextOutput<Fiberdata, &'static str> {
        let data = self.reader.records().next();
        match data {
            Some(record) => {
                let record = record.unwrap();
                let target_name = self.target_dict.get(&record.tid());
                let fiber = FiberseqData::new(record, target_name, &FiberFilters::default());
                IterNextOutput::Yield(Fiberdata::new(fiber, None))
            }
            None => IterNextOutput::Return("Ended"),
        }
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
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
    #[pyo3(get)]
    starts: Vec<Option<i64>>,
    /// Basemod reference starts
    #[pyo3(get)]
    reference_starts: Vec<Option<i64>>,
    /// Basemod ML
    #[pyo3(get)]
    ml: Vec<u8>,
}

#[pymethods]
impl Basemods {
    #[new]
    pub fn new(starts: Vec<Option<i64>>, reference_starts: Vec<Option<i64>>, ml: Vec<u8>) -> Self {
        Self {
            starts,
            reference_starts,
            ml,
        }
    }

    /// return reference ends, start + 1
    pub fn get_reference_ends(&self) -> Vec<Option<i64>> {
        self.reference_starts
            .iter()
            .map(|x| x.as_ref().map(|x| x + 1))
            .collect()
    }

    /// return ends, start + 1
    pub fn get_ends(&self) -> Vec<Option<i64>> {
        self.starts
            .iter()
            .map(|x| x.as_ref().map(|x| x + 1))
            .collect()
    }
}

/// Class for describing ranges within fiberseq data.
/// For example, nucleosomes and msp features.
#[pyclass]
#[derive(Clone)]
pub struct Ranges {
    /// Range starts
    #[pyo3(get, set)]
    pub starts: Vec<Option<i64>>,
    /// Range ends
    #[pyo3(get, set)]
    pub ends: Vec<Option<i64>>,
    /// Range lengths
    #[pyo3(get, set)]
    pub lengths: Vec<Option<i64>>,
    /// quals
    #[pyo3(get, set)]
    pub qual: Vec<u8>,
    /// Reference starts
    #[pyo3(get, set)]
    pub reference_starts: Vec<Option<i64>>,
    /// Reference ends
    #[pyo3(get, set)]
    pub reference_ends: Vec<Option<i64>>,
    /// Reference lengths
    #[pyo3(get, set)]
    pub reference_lengths: Vec<Option<i64>>,
}

#[pymethods]
impl Ranges {
    #[new]
    pub fn new(
        starts: Vec<Option<i64>>,
        lengths: Vec<Option<i64>>,
        qual: Vec<u8>,
        reference_starts: Vec<Option<i64>>,
        reference_lengths: Vec<Option<i64>>,
    ) -> Self {
        let ends = starts
            .iter()
            .zip(lengths.iter())
            .map(|(s, l)| match (s, l) {
                (Some(s), Some(l)) => Some(s + l),
                _ => None,
            })
            .collect();
        let reference_ends = reference_starts
            .iter()
            .zip(reference_lengths.iter())
            .map(|(s, l)| match (s, l) {
                (Some(s), Some(l)) => Some(s + l),
                _ => None,
            })
            .collect();
        Self {
            starts,
            ends,
            lengths,
            qual,
            reference_starts,
            reference_ends,
            reference_lengths,
        }
    }
}
