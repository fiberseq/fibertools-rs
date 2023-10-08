# `ft extract`

Extract Fiber-seq data into plain text files

## Inputs and options

See the [help message](./ft-extract-help.md) for details.

## Output description

| Column          | Description                                                     |
| --------------- | --------------------------------------------------------------- |
| ct              | Chromosome or contig                                            |
| st              | Start position of the read on the chromosome                    |
| en              | End position of the read on the chromosome                      |
| fiber           | The fiber/read name                                             |
| score           | The number of ccs passes for the read (rounded)                 |
| strand          | The strand of the read alignment                                |
| sam_flag        | The sam flag of the read alignment                              |
| HP              | The haplotype tag for the read                                  |
| RG              | The read group tag for the read                                 |
| fiber_length    | The length of the read in bp                                    |
| fiber_sequence  | The sequence of the read                                        |
| ec              | The number of ccs passes for the read (no rounding)             |
| rq              | The estimated accuracy of the read                              |
| total_AT_bp     | The total number of AT bp in the read                           |
| total_m6a_bp    | The total number of m6a bp in the read                          |
| total_nuc_bp    | The total number of nucleosome bp in the read                   |
| total_msp_bp    | The total number of MSP bp in the read                          |
| total_5mC_bp    | The total number of 5mC bp in the read                          |
| nuc_starts      | The start positions of the nucleosomes in molecular coordinates |
| nuc_lengths     | The lengths of the nucleosomes in molecular coordinates         |
| ref_nuc_starts  | The start positions of the nucleosomes in reference coordinates |
| ref_nuc_lengths | The lengths of the nucleosomes in reference coordinates         |
| msp_starts      | The start positions of the MSPs in molecular coordinates        |
| msp_lengths     | The lengths of the MSPs in molecular coordinates                |
| ref_msp_starts  | The start positions of the MSPs in reference coordinates        |
| ref_msp_lengths | The lengths of the MSPs in reference coordinates                |
| m6a             | The start positions of the m6a in molecular coordinates         |
| ref_m6a         | The start positions of the m6a in reference coordinates         |
| m6a_qual        | The quality of the m6a positions (ML value)                     |
| 5mC             | The start positions of the 5mC in molecular coordinates         |
| ref_5mC         | The start positions of the 5mC in reference coordinates         |
| 5mC_qual        | The quality of the 5mC positions (ML value)                     |
