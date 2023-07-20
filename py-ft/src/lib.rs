/// A module for using fibertools-rs
mod fiberdata;

use pyo3::prelude::*;

/// A python module for using rust to access data from fiberseq BAM files.
#[pymodule]
fn pyft(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<fiberdata::FiberdataIter>()?;
    m.add_class::<fiberdata::PyFiberdata>()?;
    m.add_class::<fiberdata::Ranges>()?;
    Ok(())
}
