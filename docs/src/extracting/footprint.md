# `ft footprint`

## Usage
```bash
ft footprint [OPTIONS] <BAM> --bed <BED> --yaml <YAML>
```
The `BAM` file is an indexed fiber-seq bam file.

The `BED` file is a bed file with the motifs you'd like to test for footprints. This should include the strand the motif is on.

The `YAML` file is a file that describes the modules within the motif that can be footprinted. e.g. a CTCF `yaml` with its multiple binding sites might look like:
```yaml
modules:
  - [0, 8]
  - [8, 16]
  - [16, 23]
  - [23, 29]
  - [29, 35]
```
Modules must start at zero, end at the length of the motif, be sorted, and be contiguous with one another. At most 15 modules are allowed, and the intervals are 0-based, half-open (like `BED`).

## Description of output columns

The footprinting output table is a tab-separated file with the same number of entries as the input BED file and the following columns:


| Column               | Description                                                        |
| -------------------- | ------------------------------------------------------------------ |
| chrom                | Chromosome                                                         |
| start               | The start position of the motif                                 |
| end               | The end position of the motif                                 |
| strand               | The strand of the motif.                              |
| n_spanning_fibers      | The number of fibers that span the motif.            |
| n_spanning_msps | The number of msp that span the motif.                     |
| n_overlapping_nucs | The number of fibers that have an intersecting nucleosome. |
| module_X | The number of fibers that are footprinted in module X. The number of module columns is determined by the footprinting yaml. |
| footprint_codes | Comma separated list of footprint codes for each fiber. See details below. |
| fire_quals | Comma separated list of fire qualities for each fiber. -1 if the MSP is not spanning or present. Note all fire_quals will be 0 or -1 if FIRE has not been applied to the bam. |
| fiber_names | Comma separated list of fiber names that span the motif. Names share the same index as the previous column, so they can be matched with footprint codes. |

## Footprint codes
The footprint codes are an encoded bit flag similar to how filtering is done with `samtools`. If the first bit is set (1) then there is an MSP that spans the footprint. For each following bit, the bit is set if that module is footprinted by that fiber.

Here are some examples in python for how you could test a footprint code in a few ways:
```python
fp_code = 0b1001 # this is a value of 9, but in binary it is 1001

# test if the first bit is set, there is a spanning MSP, true in this example
(fp_code & 1) > 0

# test if the first module is footprinted, false in this example
(fp_code & (1 << 1)) > 0 

# test if the third module is footprinted, true in this example
(fp_code & (1 << 3)) > 0
```