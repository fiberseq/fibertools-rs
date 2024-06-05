# Glossary of Fiber-seq terms

### m6A
**N6-methyladenosine**: A DNA modification made during the Fiber-seq protocol on accessible adenosines. Native m6A sites are exceedingly rare or non-existent in eukaryotic DNA.

### Inferred nucleosome 
A nucleosome inferred from the Fiber-seq data. The nucleosome position along the DNA is inferred by the `ft add-nucleosomes` algorithm using the m6A modifications (not a direct observation of a nucleosome).

### MSP
**Methyltransferase sensitive patch**: A stretch of a Fiber-seq read that has a high density of m6A sites. Specifically it is defined as a region that is inferred to be not occluded by a nucleosome. 

### FIREs
**Fiber-seq Inferred Regulatory Elements**: MSPs that are inferred to be regulatory elements based on features of the MSP including m6A density. 

### Internucleosomal linker region
Any MSP that is **not** a FIRE element.