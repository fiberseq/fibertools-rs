/// A module for using fibertools-rs
mod fiberdata;

use env_logger::{Builder, Target};
use log::LevelFilter;
use pyo3::prelude::*;

/// A python module for using rust to access data from fiberseq BAM files.
#[pymodule]
fn pyft(_py: Python, m: &PyModule) -> PyResult<()> {
    Builder::new()
        .format_timestamp_secs()
        .target(Target::Stderr)
        .filter(None, LevelFilter::Info)
        .init();

    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<fiberdata::FiberdataFetch>()?;
    m.add_class::<fiberdata::Fiberdata>()?;
    m.add_class::<fiberdata::Basemods>()?;
    m.add_class::<fiberdata::Ranges>()?;
    Ok(())
}
