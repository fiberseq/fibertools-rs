# `ft-footprint`

## description of column in the footprinting output table

The footprinting output table is a tab-separated file with the following columns:


| Column               | Description                                                        |
| -------------------- | ------------------------------------------------------------------ |
| chrom                | Chromosome                                                         |
| start               | The start position of the footprint                                 |
| end               | The end position of the footprint                                 |
| strand               | The strand of the footprint.                              |
| n_spanning_fibers      | The number of fibers that span the footprint.            |
| n_spanning_msps | The number of msp that span the footprint.                     |
| module_X | The number of fibers that are footprinted in module X. |
| footprint_codes | Comma separated list of footprint codes for each fiber. |
| fiber_names | Comma separated list of fiber names that span the footprint. |

## footprint codes
The footprint codes are a bit flag similar to how filtering is done with samtools. 

If the first bit is set (1) then there is an MSP that spans the footprint.

For each following bit the bit is set of that module is footprinted by that fiber.

Here are some examples in python for how you could test a footprint code in a few ways:
```python
fp_code = 0b1001
# test if the first bit is set
fp_code & 1

# test if the first module is footprinted, is false in this example
(fp_code & (1 << 1)) > 0 

# test if the third module is footprinted, is true in this example
(fp_code & (1 << 3)) > 0
```