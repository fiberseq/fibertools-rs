# `ft fire`

This command identifies **<ins>F</ins>iber-seq <ins>I</ins>nferred <ins>R</ins>egulatory <ins>E</ins>lements** (FIREs) from a Fiber-seq BAM. The input is a Fiber-seq BAM file with m6A and nucleosome calls and the output is a Fiber-seq bam file with the FIREs encoded in the `aq` tags.

This command can be run in isolation; however, it is usually preferable to run the [FIRE pipeline](https://github.com/fiberseq/FIRE), which runs `ft fire` and performs many additional analyses and visualizations.


[**The help page**](../help.md#ft-fire).

## Extracting from a FIRE BAM
`ft fire` can also be used as an extraction tool to extract Fiber-seq data from an already processed FIRE BAM file. 
```bash
ft fire --extract fire.bam > all.bed
```
This produces a file in BED format that contains all the MSPs, FIREs, and nucleosomes in the FIRE BAM file. This command produces output analogous to the [now removed](https://github.com/fiberseq/FIRE/issues/24) `model.results.bed.gz` result from older versions of FIRE pipeline. 
