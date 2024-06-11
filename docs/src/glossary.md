# Glossary of Fiber-seq terms

### Fiber-seq protocol
A treatment of DNA chromatin fibers with the Hia5 enzyme which preferentially methylates accessible adenosines. The fibers are then sequenced using a long-read platform that can detect the [m6A](#m6a) modifications. 

### Fiber-seq read or fiber
A DNA fiber that has been sequenced using the Fiber-seq protocol and has [m6A](#m6a) modifications called. 

### m6A
**N6-methyladenosine**: A DNA modification made during the Fiber-seq protocol on accessible adenosines. Native [m6A](#m6a)  sites are exceedingly rare or non-existent in eukaryotic DNA.

### Inferred nucleosome 
A nucleosome inferred from the Fiber-seq data. The nucleosome position along the DNA is inferred by the `ft add-nucleosomes` algorithm using the [m6A](#m6a) modifications (not a direct observation of a nucleosome).

### MSP
**Methyltransferase sensitive patch**: A stretch of a Fiber-seq read that has a high density of [m6A](#m6a) sites. Specifically it is defined as a region that is inferred to not be occluded by a nucleosome. 

### FIRE
**Fiber-seq Inferred Regulatory Element**: A [MSP](#msp) that is inferred to be a regulatory element based on features of the [MSP](#msp) including [m6A](#m6a) density. 

### Internucleosomal linker region
Any [MSPs](#msp) that is **not** a [FIRE](#fires) element.