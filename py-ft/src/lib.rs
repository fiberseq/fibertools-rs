/// A module for using fibertools-rs
pub mod fiberdata;

use pyo3::prelude::*;

/// A Python module implemented in Rust.
#[pymodule]
fn py_ft(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<fiberdata::FiberdataIter>()?;
    m.add_class::<fiberdata::PyFiberdata>()?;
    Ok(())
}
