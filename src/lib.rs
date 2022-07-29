/// Command line interface for fibertools-rs.
pub mod bamlift;
pub mod center;
pub mod cli;
pub mod extract;

use anyhow::Result;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;
// MY IO TOOLS
const BUFFER_SIZE: usize = 32 * 1024;

/// Get a buffered output writer from stdout or a file
fn get_output(path: Option<PathBuf>) -> Result<Box<dyn Write + Send + 'static>> {
    let writer: Box<dyn Write + Send + 'static> = match path {
        Some(path) => {
            if path.as_os_str() == "-" {
                Box::new(BufWriter::with_capacity(BUFFER_SIZE, io::stdout()))
            } else {
                Box::new(BufWriter::with_capacity(BUFFER_SIZE, File::create(path)?))
            }
        }
        None => Box::new(BufWriter::with_capacity(BUFFER_SIZE, io::stdout())),
    };
    Ok(writer)
}

/// unzip a vector of tupples
pub fn unzip_to_vectors<T>(vec: Vec<(T, T)>) -> (Vec<T>, Vec<T>) {
    vec.into_iter().unzip()
}

/// join a vector with commas
pub fn join_by_str<'a, I, Z>(vals: I, sep: &str) -> String
where
    I: IntoIterator<Item = Z>,
    Z: ToString + 'a,
{
    vals.into_iter().map(|v| v.to_string() + sep).collect()
}

/// Write to stdout if - or the file specified by a path
pub fn writer(filename: &str) -> Result<Box<dyn Write>> {
    //let ext = Path::new(filename).extension();
    let path = PathBuf::from(filename);
    let buffer = get_output(Some(path))?; //.expect("Error: cannot create output file");
    Ok(buffer)
}
pub struct FiberOutFiles {
    pub m6a: Option<Box<dyn Write>>,
    pub cpg: Option<Box<dyn Write>>,
    pub msp: Option<Box<dyn Write>>,
    pub nuc: Option<Box<dyn Write>>,
    pub all: Option<Box<dyn Write>>,
}

impl FiberOutFiles {
    pub fn new(
        m6a: &Option<String>,
        cpg: &Option<String>,
        msp: &Option<String>,
        nuc: &Option<String>,
        all: &Option<String>,
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
        Ok(FiberOutFiles {
            m6a,
            cpg,
            msp,
            nuc,
            all,
        })
    }
}
