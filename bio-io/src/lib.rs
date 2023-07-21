use anyhow::Result;
use colored::Colorize;
use gzp::deflate::Bgzf; //, Gzip, Mgzip, RawDeflate};
use gzp::{Compression, ZBuilder};
use itertools::Itertools;
use lazy_static::lazy_static;
use niffler::get_reader;
use regex::Regex;
use rust_htslib::bam;
use rust_htslib::bam::record::Aux;
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
    let out = write!(buffer, "{}", out);
    if let Err(err) = out {
        if err.kind() == io::ErrorKind::BrokenPipe {
            exit(0);
        } else {
            panic!("Error: {}", err);
        }
    }
}

/// a reader that can read compressed files but also stdin (indicated by -)
/// ```
/// use bio_io::buffer_from;
/// use std::io;
/// let reader = buffer_from("../tests/data/test.txt.gz").expect("Error: cannot open file");
/// let msg = io::read_to_string(reader).unwrap();
/// assert_eq!(msg, "Hello World!\n");
/// let reader = buffer_from("../tests/data/test.txt").expect("Error: cannot open file");
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
    log::debug!("Compression: {:?}", compression);
    let buffer = BufReader::new(reader);
    Ok(buffer)
}

/*
BAM IO
*/

/// Write to a bam file.
pub fn program_bam_writer(
    out: &str,
    template_bam: &bam::Reader,
    threads: usize,
    program_name: &str,
    program_id: &str,
    program_version: &str,
) -> bam::Writer {
    let mut header = bam::Header::from_template(template_bam.header());

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
        log::trace!("last program {}", last_program);
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
    let mut out = if out == "-" {
        bam::Writer::from_stdout(&header, bam::Format::Bam).unwrap()
    } else {
        bam::Writer::from_path(out, &header, bam::Format::Bam).unwrap()
    };
    out.set_threads(threads).unwrap();
    out
}

/// Open bam file
pub fn bam_reader(bam: &str, threads: usize) -> bam::Reader {
    let mut bam = if bam == "-" {
        bam::Reader::from_stdin().unwrap_or_else(|_| panic!("Failed to open bam from stdin"))
    } else {
        bam::Reader::from_path(bam).unwrap_or_else(|_| panic!("Failed to open {}", bam))
    };
    bam.set_threads(threads).unwrap();
    bam
}
// This is a bam chunk reader
pub struct BamChunk<'a> {
    pub bam: bam::Records<'a, bam::Reader>,
    pub chunk_size: usize,
}

// The `Iterator` trait only requires a method to be defined for the `next` element.
impl<'a> Iterator for BamChunk<'a> {
    // We can refer to this type using Self::Item
    type Item = Vec<bam::Record>;

    // The return type is `Option<T>`:
    //     * When the `Iterator` is finished, `None` is returned.
    //     * Otherwise, the next value is wrapped in `Some` and returned.
    // We use Self::Item in the return type, so we can change
    // the type without having to update the function signatures.
    fn next(&mut self) -> Option<Self::Item> {
        let start = Instant::now();
        let mut cur_vec = vec![];
        for r in self.bam.by_ref().take(self.chunk_size) {
            cur_vec.push(r.unwrap())
        }
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
            ("102-739-100".to_string(), PbChem::Revio)
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
        log::warn!(
            "Polymerase for BINDINGKIT={} not found. Defaulting to ML model made for REVIO.",
            binding_kit
        );
        &PbChem::Revio
    });

    // log the chem being used
    let chem = match chemistry {
        PbChem::Two => "2.0",
        PbChem::TwoPointTwo => "2.2",
        PbChem::ThreePointTwo => "3.2",
        PbChem::Revio => "Revio",
    };
    log::info!(
        "Bam header implies PacBio chemistry {} binding kit {}.",
        chem,
        binding_kit
    );
    chemistry.clone()
}

///```
/// use rust_htslib::{bam, bam::Read};
/// use log;
/// let mut bam = bam::Reader::from_path(&"../tests/data/all.bam").unwrap();
/// for record in bam.records() {
///     let record = record.unwrap();
///     let n_s = bio_io::get_u32_tag(&record, b"ns");
///     let n_l = bio_io::get_u32_tag(&record, b"nl");
///     let a_s = bio_io::get_u32_tag(&record, b"as");
///     let a_l = bio_io::get_u32_tag(&record, b"al");
///     log::debug!("{:?}", a_s);
/// }
///```
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
