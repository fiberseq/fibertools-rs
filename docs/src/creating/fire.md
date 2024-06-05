# `ft fire`

This command identifies **<ins>F</ins>iber-seq <ins>I</ins>nferred <ins>R</ins>egulatory <ins>E</ins>lements** (FIREs) from a Fiber-seq BAM. The input is a Fiber-seq BAM file with m6A and nucleosome calls and the output is a Fiber-seq bam file with the FIREs encoded in the `aq` tags.

This command can be run in isolation; however, it is usually preferable to run the [FIRE pipeline](https://github.com/fiberseq/FIRE), which runs `ft fire` and performs many additional analyses and visualizations.


[**The help page**](../help.md#ft-fire).