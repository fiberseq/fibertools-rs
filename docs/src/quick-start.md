# Quick start for `fibertools`

<!-- toc -->

# Fiber-seq starting with PacBio

HiFi kinetics are required for predicting [m6A](glossary.md#m6a) with `fibertools`. **Check with your sequencing provider prior to sequencing** to ensure that the output file will have kinetics. Additionally, many of `fibertools` commands are compatible with CpG methylation which can be completed on instrument (if requested) or later with [Jasmine](https://github.com/PacificBiosciences/jasmine), e.g. `jasmine --keep-kinetics input.ccs.bam output.ccs.bam`. This command should be run prior to `fibertools` if CpG methylation information is desired as Jasmine will overwrite the m6A predictions in the MM and ML tags.


### Predict m6A and infer nucleosomes
To create useable Fiber-seq data you must first call m6A base-mods on the PacBio CCS bam using `fibertools`. First [install fibertools](install/install.md) and then process your bam file using the prediction command. 

```bash
ft predict-m6a -t 16 input.ccs.bam output.fiberseq.bam 
```
This will both make m6A calls and identify [nucleosomes](glossary.md#inferred-nucleosome) on each [fiber](glossary.md#fiber-seq-read-or-fiber).

**Note**, the input **CCS bam must have average kinetics** to be able to call m6A. 


### Alignment and phasing
We recommend aligning with [pbmm2](https://github.com/PacificBiosciences/pbmm2) and phasing with [HiPhase](https://github.com/PacificBiosciences/hiphase). Please see their respective documentation for more information.

Alternatively, we have written a [snakemake pipeline](https://github.com/mrvollger/k-mer-variant-phasing) to align and phase Fiber-seq data; however, this pipeline is not officially supported.

After this point you will have a Fiber-seq BAM file that is compatible with all the [extraction](extracting/extracting.md) commands in `fibertools`.

### Fiber-seq peaks and UCSC browser tracks
Once you have a phased bam file, you can identify [Fiber-seq inferred regulatory elements (FIREs)](glossary.md#fires) to call Fiber-seq peaks and make a UCSC trackHub. Please see the [FIRE repository](https://github.com/fiberseq/FIRE) for more details.

# Fiber-seq starting with Oxford Nanopore Technologies (ONT)

### Predict m6A

**ft predict-m6a** does not include a model for ONT data; however, you can use software, such as [Dorado](https://github.com/nanoporetech/dorado), to add CpG and m6A to your ONT BAM file.

### Infer nucleosomes and MSPs

Once you have CpG and m6A information in your ONT BAM file, you can use `ft add-nucleosomes` to infer nucleosomes and MSPs. With Dorado, we find the best results when restricting to m6A modifications with an ML score of 250 or higher.
```bash
ft add-nucleosomes --ml 250 input.bam output.bam
```

### Alignment and phasing 
You can either use [Dorado](https://github.com/nanoporetech/dorado) to align your ONT data or use a tool like [minimap2](https://github.com/lh3/minimap2) to align your data. If you do use `minimap2` be sure to include the flag `-y` to preserve the CpG and m6A information in the output BAM file.

We recommend using [WhatsHap](https://whatshap.readthedocs.io/en/latest/) for phasing ONT data. Please see their documentation for more information.

After this point you will have a Fiber-seq BAM file that is compatible with all the [extraction](extracting/extracting.md) commands in `fibertools`.

### Fiber-seq peaks and UCSC browser tracks
Some users report reasonable success in applying the [FIRE pipeline](https://github.com/fiberseq/FIRE) to ONT data. However, **please note that ..FIRE models were not trained or validated for ONT data.** With that said all the instructions for applying the ..FIRE pipeline to a PacBio BAM should work for an ONT BAM as well.