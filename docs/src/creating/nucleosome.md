# `ft add-nucleosomes`

`ft add-nucleosomes` adds nucleosome positions and MSP positions to the Fiber-seq bam. The input is a Fiber-seq bam file with m6A calls and the output is a Fiber-seq bam file with the nucleosome positions encoded in the `ns` and `nl` tags and the MSP positions encoded in the `as` and `al` tags.

It is usually unnecessary to run this command manually, as it is called by `ft predict-m6a` automatically. However, you can run it manually if you want to adjust the parameters of the nucleosome calling algorithm.

[**The help page**](../help.md#ft-add-nucleosomes)