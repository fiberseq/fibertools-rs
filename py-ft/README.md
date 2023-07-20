# pyft: python bindings for fibertools-rs
[![Documentation Status](https://readthedocs.org/projects/py-ft/badge/?version=latest)](https://py-ft.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://badge.fury.io/py/pyft.svg)](https://badge.fury.io/py/pyft)

`pyft` provides a python API for the rust library [fibertools-rs](https://github.com/fiberseq/fibertools-rs). The inspiration for this API is to make analysis in python easier and faster; therefore, only extraction of data from a fiberseq bam is supported and not writing. 

# Install
```bash
pip install pyft
```

# Example
```python
import pyft
bam_f = "../tests/data/center.bam"
fiberdata = pyft.FiberdataFetch(bam_f, "chr1", 0, 1000)
for fiber in fiberdata:
    # print some info about the fiber
    print(fiber)
    # print the fiber number of ccs passes
    print(fiber.ec)    
    # print the mps start positions
    print(fiber.msp.starts)    
```