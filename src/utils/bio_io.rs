use anyhow::Result;
use colored::Colorize;
use gzp::deflate::Bgzf; //, Gzip, Mgzip, RawDeflate};
use gzp::{Compression, ZBuilder};
use indicatif::{ProgressBar, ProgressStyle};
use itertools::Itertools;
use lazy_static::lazy_static;
use linear_map::LinearMap;
use niffler::get_reader;
use rayon::current_num_threads;
use regex::Regex;
use rust_htslib::bam;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::Header;
use rust_htslib::bam::Read;
use std::collections::HashMap;
use std::env;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{self, stdout, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::exit;
use std::time::Instant;

const BUFFER_SIZE: usize = 32 * 1024;
const COMPRESSION_THREADS: usize = 8;
const COMPRESSION_LEVEL: u32 = 6;

/*
PROGRESS BARS
*/
pub fn no_length_progress_bar() -> ProgressBar {
    let progress_style_no_length = format!(
        "{}    {} {} {} {} {}",
        "[{elapsed_precise:.yellow.bold}]",
        "Records".cyan(),
        "{per_sec:<15.cyan}",
        "Read".blue(),
        "{human_pos:.blue}",
        "records".blue(),
    );
    let style = ProgressStyle::default_bar()
        .template(&progress_style_no_length)
        .unwrap();
    let bar = ProgressBar::new(0);
    bar.set_style(style);
    let finish = indicatif::ProgressFinish::AndLeave;
    bar.with_finish(finish)
}

/*
STANDARD FILE IO
*/

/// Get a buffered output writer from stdout or a file
fn get_output(path: Option<PathBuf>) -> Result<Box<dyn Write + Send + 'static>> {
    let writer: Box<dyn Write + Send + 'static> = match path {
        Some(path) => {
            if path.as_os_str() == "-" {
                Box::new(BufWriter::with_capacity(BUFFER_SIZE, stdout()))
            } else {
                Box::new(BufWriter::with_capacity(BUFFER_SIZE, File::create(path)?))
            }
        }
        None => Box::new(BufWriter::with_capacity(BUFFER_SIZE, stdout())),
    };
    Ok(writer)
}

/// Write normal or compressed files seamlessly
/// Uses the presence of a `.gz` extension to decide
// Attempting to have a file writer too
pub fn writer(filename: &str) -> Result<Box<dyn Write>> {
    let ext = Path::new(filename).extension();
    let path = PathBuf::from(filename);
    let buffer = get_output(Some(path))?;
    if ext == Some(OsStr::new("gz")) || ext == Some(OsStr::new("bgz")) {
        let writer = ZBuilder::<Bgzf, _>::new()
            .num_threads(COMPRESSION_THREADS)
            .compression_level(Compression::new(COMPRESSION_LEVEL))
            .from_writer(buffer);
        Ok(Box::new(writer))
    } else {
        Ok(buffer)
    }
}

/// write to a file, but don't error on broken pipes
pub fn write_to_file(out: &str, buffer: &mut Box<dyn Write>) {
    let out = write!(buffer, "{out}");
    if let Err(err) = out {
        if err.kind() == io::ErrorKind::BrokenPipe {
            exit(0);
        } else {
            panic!("Error: {err}");
        }
    }
}

/// write a BAM record, but don't error on broken pipes
pub fn write_record(writer: &mut bam::Writer, record: &bam::Record) -> anyhow::Result<()> {
    match writer.write(record) {
        Ok(_) => Ok(()),
        Err(err) => {
            // Check for the specific WriteRecord error that indicates broken pipe
            if matches!(err, rust_htslib::errors::Error::WriteRecord) {
                let error_msg = err.to_string();
                if error_msg.contains("failed to write BAM/BCF record (out of disk space?)") {
                    exit(0);
                }
            }

            // Convert to anyhow error first to access the error chain
            let anyhow_err: anyhow::Error = err.into();

            // Also check for actual BrokenPipe IO errors in the error chain
            let mut current_err = anyhow_err.source();
            while let Some(err) = current_err {
                if let Some(io_err) = err.downcast_ref::<io::Error>() {
                    if io_err.kind() == io::ErrorKind::BrokenPipe {
                        exit(0);
                    }
                }
                current_err = err.source();
            }

            // For all other errors, return them normally
            Err(anyhow_err)
        }
    }
}

/// a reader that can read compressed files but also stdin (indicated by -)
/// ```
/// use fibertools_rs::utils::bio_io::buffer_from;
/// use std::io;
/// let reader = buffer_from("tests/data/test.txt.gz").expect("Error: cannot open file");
/// let msg = io::read_to_string(reader).unwrap();
/// assert_eq!(msg, "Hello World!\n");
/// let reader = buffer_from("tests/data/test.txt").expect("Error: cannot open file");
/// let msg = io::read_to_string(reader).unwrap();
/// assert_eq!(msg, "Hello World!\n");
/// ```
pub fn buffer_from<P: AsRef<Path>>(
    path: P,
) -> Result<BufReader<Box<dyn std::io::Read>>, anyhow::Error> {
    let path = path.as_ref();
    let readable: Box<dyn std::io::Read> = if path == Path::new("-") {
        Box::new(io::BufReader::new(io::stdin()))
    } else {
        Box::new(io::BufReader::new(std::fs::File::open(path)?))
    };
    let (reader, compression) = get_reader(readable)?;
    log::debug!("Compression: {compression:?}");
    let buffer = BufReader::new(reader);
    Ok(buffer)
}

/*
BAM IO
*/

pub fn program_bam_writer_from_header(
    out: &str,
    mut header: bam::Header,
    program_name: &str,
    program_id: &str,
    program_version: &str,
) -> bam::Writer {
    // add to the header
    let header_string = String::from_utf8_lossy(&header.to_bytes()).to_string();
    let mut header_rec = bam::header::HeaderRecord::new(b"PG");
    // ID
    let tag = format!("PN:{program_name}");
    let program_count = header_string.matches(&tag).count();
    let id_tag = format!("{}.{}", program_id, program_count + 1);
    header_rec.push_tag(b"ID", &id_tag);
    // PN
    header_rec.push_tag(b"PN", program_name);
    // PP
    let re_pp = Regex::new(r"@PG\tID:([^\t]+)").unwrap();
    let last_program = re_pp.captures_iter(&header_string).last();
    if let Some(last_program) = last_program {
        let last_program = last_program[1].to_string();
        log::trace!("last program {last_program}");
        header_rec.push_tag(b"PP", &last_program);
    };
    // VN
    header_rec.push_tag(b"VN", program_version);
    let cli = env::args().join(" ");
    // CL
    header_rec.push_tag(b"CL", &cli);
    header.push_record(&header_rec);
    log::trace!("{:?}", String::from_utf8_lossy(&header.to_bytes()));

    // make the writer
    if out == "-" {
        bam::Writer::from_stdout(&header, bam::Format::Bam).unwrap()
    } else {
        bam::Writer::from_path(out, &header, bam::Format::Bam).unwrap()
    }
}

/// Write to a bam file.
pub fn program_bam_writer(
    out: &str,
    template_bam: &bam::Reader,
    program_name: &str,
    program_id: &str,
    program_version: &str,
) -> bam::Writer {
    let header = bam::Header::from_template(template_bam.header());
    program_bam_writer_from_header(out, header, program_name, program_id, program_version)
}

/// Open bam file
pub fn bam_reader(bam: &str) -> bam::Reader {
    if bam == "-" {
        bam::Reader::from_stdin().unwrap_or_else(|_| panic!("Failed to open bam from stdin"))
    } else {
        bam::Reader::from_path(bam).unwrap_or_else(|_| panic!("Failed to open {bam}"))
    }
}
// This is a bam chunk reader
pub struct BamChunk<'a> {
    pub bam: bam::Records<'a, bam::Reader>,
    pub chunk_size: usize,
    pub pre_chunk_done: u64,
    pub bar: ProgressBar,
    pub bit_flag_filter: u16,
}

impl<'a> BamChunk<'a> {
    pub fn new(bam: bam::Records<'a, bam::Reader>, chunk_size: Option<usize>) -> Self {
        let chunk_size = std::cmp::min(
            chunk_size.unwrap_or_else(|| current_num_threads() * 100),
            2_500,
        );
        let bar = no_length_progress_bar();
        Self {
            bam,
            chunk_size,
            pre_chunk_done: 0,
            bar,
            bit_flag_filter: 0,
        }
    }

    pub fn set_bit_flag_filter(&mut self, bit_flag: u16) {
        self.bit_flag_filter = bit_flag;
    }
}

// The `Iterator` trait only requires a method to be defined for the `next` element.
impl Iterator for BamChunk<'_> {
    // We can refer to this type using Self::Item
    type Item = Vec<bam::Record>;

    // The return type is `Option<T>`:
    //     * When the `Iterator` is finished, `None` is returned.
    //     * Otherwise, the next value is wrapped in `Some` and returned.
    // We use Self::Item in the return type, so we can change
    // the type without having to update the function signatures.
    fn next(&mut self) -> Option<Self::Item> {
        // update progress bar with results from previous iteration
        self.bar.inc(self.pre_chunk_done);

        let start = Instant::now();
        let mut cur_vec = vec![];
        for r in self.bam.by_ref().take(self.chunk_size) {
            let r = r.unwrap();
            let has_mm_and_ml = r.aux(b"MM").is_ok() && r.aux(b"ML").is_ok();
            if has_mm_and_ml
                && (r.cigar().leading_hardclips() > 0 || r.cigar().trailing_hardclips() > 0)
            {
                log::warn!(
                    "Skipping read ({}) because it has been hard clipped and has ML and MM tags. This read will be excluded from calculations and any output.",
                    String::from_utf8_lossy(r.qname())
                );
                continue;
            }
            // filter by bit flag
            if r.flags() & self.bit_flag_filter != 0 {
                continue;
            }
            cur_vec.push(r);
        }

        // extend progress bar
        self.pre_chunk_done = cur_vec.len() as u64;
        self.bar.inc_length(self.pre_chunk_done);

        // return
        if cur_vec.is_empty() {
            None
        } else {
            let duration = start.elapsed().as_secs_f64();
            log::debug!(
                "Read {} bam records at {}.",
                format!("{:}", cur_vec.len()).bright_magenta().bold(),
                format!("{:.2?} reads/s", cur_vec.len() as f64 / duration)
                    .bright_cyan()
                    .bold(),
            );
            Some(cur_vec)
        }
    }
}

#[derive(Clone, Debug)]
pub enum PbChem {
    Two,
    TwoPointTwo,
    ThreePointTwo,
    Revio,
}

pub fn find_pb_polymerase(header: &bam::Header) -> PbChem {
    lazy_static! {
        static ref CHEMISTRY_MAP: HashMap<String, PbChem> = HashMap::from([
            // regular 2.0
            ("101-789-500".to_string(), PbChem::Two),
            // this is actually 2.1 we didn't have training data for it
            ("101-820-500".to_string(), PbChem::Two),
            // regular 2.2
            ("101-894-200".to_string(), PbChem::TwoPointTwo),
            // is really 3.1 but has polymerase of 2.1, and we need to make that 2.0
            ("102-194-200".to_string(), PbChem::Two),
            // regular 3.2
            ("102-194-100".to_string(), PbChem::ThreePointTwo),
            // Revio has kinetics most similar to 2.2
            ("102-739-100".to_string(), PbChem::Revio),
            // Vega currently using Revio chemistry
            ("103-426-500".to_string(), PbChem::Revio)
        ]);
    }
    lazy_static! {
        static ref MM_DS: regex::Regex =
            regex::Regex::new(r".*READTYPE=([^;]+);.*BINDINGKIT=([^;]+);").unwrap();
    }
    let z = header.to_hashmap();
    let rg = z.get("RG").expect("RG tag missing from bam file");
    let mut read_type = "";
    let mut binding_kit = "";
    for tag in rg {
        for (tag, val) in tag {
            if tag == "DS" {
                for cap in MM_DS.captures_iter(val) {
                    read_type = cap.get(1).map_or("", |m| m.as_str());
                    binding_kit = cap.get(2).map_or("", |m| m.as_str());
                }
            }
        }
    }
    // force revio model
    if env::var("FT_REVIO").is_ok() {
        binding_kit = "102-739-100";
    }
    // make sure read-type is CCS
    assert_eq!(read_type, "CCS");
    // grab chemistry
    let chemistry = CHEMISTRY_MAP.get(binding_kit).unwrap_or_else(|| {
        log::error!("Model for BINDINGKIT={binding_kit} not available. Unable to run predictions.");
        std::process::exit(1);
    });

    // log the chem being used
    let chem = match chemistry {
        PbChem::Two => "2.0",
        PbChem::TwoPointTwo => "2.2",
        PbChem::ThreePointTwo => "3.2",
        PbChem::Revio => "Revio",
    };
    log::info!("Bam header implies PacBio chemistry {chem} binding kit {binding_kit}.");
    chemistry.clone()
}

pub fn get_u32_tag(record: &bam::Record, tag: &[u8; 2]) -> Vec<i64> {
    if let Ok(Aux::ArrayU32(array)) = record.aux(tag) {
        let read_array = array.iter().map(|x| x as i64).collect::<Vec<_>>();
        read_array
    } else {
        vec![]
    }
}

pub fn get_u8_tag(record: &bam::Record, tag: &[u8; 2]) -> Vec<u8> {
    if let Ok(Aux::ArrayU8(array)) = record.aux(tag) {
        let read_array = array.iter().collect::<Vec<_>>();
        read_array
    } else {
        vec![]
    }
}

/// Convert the PacBio u16 tag into the u8 tag we normally expect.
pub fn get_pb_u16_tag_as_u8(record: &bam::Record, tag: &[u8; 2]) -> Vec<u8> {
    if let Ok(Aux::ArrayU16(array)) = record.aux(tag) {
        let read_array = array.iter().collect::<Vec<_>>();
        read_array
            .iter()
            .map(|&x| {
                if x < 64 {
                    x
                } else if x < 191 {
                    (x - 64) / 2 + 64
                } else if x < 445 {
                    (x - 192) / 4 + 128
                } else if x < 953 {
                    (x - 448) / 8 + 192
                } else {
                    255
                }
            })
            .map(|x| x as u8)
            .collect()
    } else {
        vec![]
    }
}

pub fn get_f32_tag(record: &bam::Record, tag: &[u8; 2]) -> Vec<f32> {
    if let Ok(Aux::ArrayFloat(array)) = record.aux(tag) {
        let read_array = array.iter().collect::<Vec<_>>();
        read_array
    } else {
        vec![]
    }
}

pub fn header_from_hashmap(hash_header: HashMap<String, Vec<LinearMap<String, String>>>) -> Header {
    let mut header = Header::new();
    for (key, values) in hash_header.iter() {
        for value in values {
            let mut record = HeaderRecord::new(key.as_bytes());
            for (tag, val) in value.iter() {
                record.push_tag(tag.as_bytes(), val);
            }
            header.push_record(&record);
        }
    }
    header
}

/// Convert a BAM header to a string representation including comments
///
/// # Arguments
///
/// * `header_view` - The BAM HeaderView to convert
///
/// # Returns
///
/// * String representation of the header including any comments
pub fn bam_header_to_string(header_view: &bam::HeaderView) -> String {
    // Create a Header from the HeaderView to access full functionality
    let header = bam::Header::from_template(header_view);

    // Use to_hashmap to get structured header data
    let header_hashmap = header.to_hashmap();

    let mut header_lines = Vec::new();

    // Process header records in a specific order: HD first, then SQ, then others
    let record_order = ["HD", "SQ", "RG", "PG", "CO"];

    for record_type in &record_order {
        if let Some(records) = header_hashmap.get(*record_type) {
            for record in records {
                let mut line = format!("@{}", record_type);
                for (key, value) in record {
                    line.push_str(&format!("\t{}:{}", key, value));
                }
                header_lines.push(line);
            }
        }
    }

    // Add any remaining record types not in the standard order
    for (record_type, records) in &header_hashmap {
        if !record_order.contains(&record_type.as_str()) {
            for record in records {
                let mut line = format!("@{}", record_type);
                for (key, value) in record {
                    line.push_str(&format!("\t{}:{}", key, value));
                }
                header_lines.push(line);
            }
        }
    }

    // Add comments
    for comment in header.comments() {
        header_lines.push(format!("@CO\t{}", comment));
    }

    // Join all lines with newlines and add final newline
    let mut result = header_lines.join("\n");
    result.push('\n');
    result
}

/// Converts seq to uppercase leaving characters other than a,c,g,t,n unchanged.
///
/// # Arguments
///
/// * `seq` - Input vector of bases in the sequence, all u8.
///
/// # Returns
///
/// * Output vector same as `seq` with a,c,g,t,n converted to uppercase.
///
/// # Example
///
/// ```
/// use fibertools_rs::utils::bio_io::convert_seq_uppercase;
/// let x = vec![b'A', b'C', b'G', b'T', b'N', b'a', b'c', b'g', b't', b'n', b'='];
/// let out = convert_seq_uppercase(x);
/// let expected_out = vec![b'A', b'C', b'G', b'T', b'N', b'A', b'C', b'G', b'T', b'N', b'='];
/// let matching = out.iter().zip(expected_out.iter()).filter(|&(a, b)| a == b).count();
/// assert!(matching == out.len() && matching == expected_out.len());
/// ```
pub fn convert_seq_uppercase(mut seq: Vec<u8>) -> Vec<u8> {
    for base in &mut seq {
        match *base {
            b'a' => *base = b'A',
            b'c' => *base = b'C',
            b'g' => *base = b'G',
            b't' => *base = b'T',
            b'n' => *base = b'N',
            _ => {}
        }
    }
    seq
}
