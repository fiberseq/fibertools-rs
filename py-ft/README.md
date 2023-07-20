# py_ft: Python bindings for fibertools-rs

```python
import py_ft
bam_f = "../tests/data/center.bam"
fiberdata = py_ft.py_ft.FiberdataIter(bam_f, "chr1", 0, 1000)
for fiber in fiberdata:
    # print some info about the fiber
    print(fiber)
    # print the fiber number of ccs passes
    print(fiber.ec)    
```