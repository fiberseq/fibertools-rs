use bio_io;
use pyo3::iter::IterNextOutput;
use pyo3::prelude::*;
use rust_htslib::bam;
//use rust_htslib::bam::Read;
//use std::io::Result;

#[pyclass]
struct PyClassIter {
    count: usize,
    bam: bam::Reader,
}

#[pymethods]
impl PyClassIter {
    #[new]
    pub fn new(f: &str) -> Self {
        let bam = bio_io::bam_reader(f, 8);
        PyClassIter { count: 0, bam }
    }

    fn __next__(&mut self) -> IterNextOutput<usize, &'static str> {
        if self.count < 5 {
            self.count += 1;
            // Given an instance `counter`, First five `next(counter)` calls yield 1, 2, 3, 4, 5.
            IterNextOutput::Yield(self.count)
        } else {
            // At the sixth time, we get a `StopIteration` with `'Ended'`.
            //     try:
            //         next(counter)
            //     except StopIteration as e:
            //         assert e.value == 'Ended'
            IterNextOutput::Return("Ended")
        }
    }
}
