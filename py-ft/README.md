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
# bam file with ~3X coverage of chr20
bam_f = "../tmp.bam"
fiberdata = pyft.FiberdataFetch(bam_f, "chr20", 0, 10_000_000)
for idx, fiber in enumerate(fiberdata):
    if idx < 10:
        # print some info about the fiber
        print(fiber)
    # the number of ccs passes
    fiber.ec
    # the mps start positions
    fiber.msp.starts
    # print the nuc reference starts
    fiber.nuc.reference_starts
    # lift query (fiber) positions to reference positions
    fiber.lift_query_positions([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    # lift reference positions to query (fiber) positions
    fiber.lift_reference_positions([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    # aligned blocks between query (fiber) and reference 
    fiber.get_aligned_blocks_as_ranges()
```

# Documentation
The documentation for `pyft` can be found on [readthedocs](https://py-ft.readthedocs.io/en/latest/).