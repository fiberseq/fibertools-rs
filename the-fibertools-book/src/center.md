# `ft-center`

This command centers Fiber-seq data around given reference positions. This is useful for making aggregate m6A and CpG observations, as well as visualization of SVs

## Inputs and options

See the [help message](./ft-center-help.md) for details.

## Output description

This command writes Fiber-seq data in a tab-delimited format to `stdout` that has been centered relative to positions specified in the input bed file.

| Column               | Description                                                        |
| -------------------- | ------------------------------------------------------------------ |
| chrom                | Chromosome                                                         |
| centering_position   | The position on the chromosome about which data is being centered. |
| strand               | The strand of the centering position.                              |
| subset_sequence      | The sequence of the read around the centering position.            |
| reference_start      | The start position of the read in the reference.                   |
| reference_end        | The end position of the read in the reference.                     |
| query_name           | The name of the sequencing read                                    |
| RG                   | The read group the read belongs to                                 |
| centered_query_start | The start position of the read relative to the centering position  |
| centered_query_end   | The end position of the read relative to the centering position    |
| query_length         | The length of the read                                             |

### Additional columns specific to long (default) format

| Column                 | Description                                                                |
| ---------------------- | -------------------------------------------------------------------------- |
| centered_position_type | The type of position being centered. One of: `m6A`, `5mC`, `nuc`, `msp`.   |
| centered_start         | The start position of the "feature" relative to the centering position     |
| centered_end           | The end position of the "feature" relative to the centering position       |
| centered_quality       | The quality of the "feature" relative to the centering position (ML value) |

### Additional columns specific to the wide format

| Column                 | Description                                 |
| ---------------------- | ------------------------------------------- |
| centered_m6a_positions | A comma separated list of m6a positions     |
| m6a_qual               | The quality of the m6a positions (ML value) |
| centered_5mC_positions | A comma separated list of 5mC positions     |
| 5mC_qual               | The quality of the 5mC positions (ML value) |
| centered_nuc_starts    | A comma separated list of nuc starts        |
| centered_nuc_ends      | A comma separated list of nuc ends          |
| centered_msp_starts    | A comma separated list of msp starts        |
| centered_msp_ends      | A comma separated list of msp ends          |
| query_sequence         | The sequence of the read                    |

Note if the `--wide` flag is used with the `--reference` flag some positions in the comma separated lists can be `NA` when the reference sequence has an insertion or deletion relative to the read sequence at that position.
