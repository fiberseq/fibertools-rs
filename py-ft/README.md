# py-ft: Python bindings for fibertools-rs
[![Documentation Status](https://readthedocs.org/projects/py-ft/badge/?version=latest)](https://py-ft.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://badge.fury.io/py/py-ft.svg)](https://badge.fury.io/py/py-ft)

## Example
```python
import py_ft
bam_f = "../tests/data/center.bam"
fiberdata = py_ft.FiberdataIter(bam_f, "chr1", 0, 1000)
for fiber in fiberdata:
    # print some info about the fiber
    print(fiber)
    # print the fiber number of ccs passes
    print(fiber.ec)    
```

## Docs
See the [docs](https://py-ft.readthedocs.io/en/latest/).