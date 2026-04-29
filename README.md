# Molecular annotation tags

With the advent of functional single-molecule sequencing, there is an emerging need to annotate non-genetic elements on individual DNA molecules. This is handled for base modifications using the MM and ML tags. However, we have a need for a more extensible format beyond base modifications that allows for generic annotation of any segment of a molecule. To accomplish this we are proposing the `MA` tag, and we describe the format of this tag below.

Integration of this spec or a similar one into the SAM/BAM/CRAM is needed to standardize the representation of molecular annotations across different sequencing platforms and analysis tools. Importantly, a standardized format will allow for integration in a wide range of genomic tools, e.g. visualization in genome browsers and allow developers to be justified in spending time integrating this format into their tools.

## Format Specification

### Structure

The molecular annotation format uses three related tags:

- **MA:Z:** - Read length followed by annotation positions with lengths (required)
- **AQ:B:C** - Annotation Quality scores (optional, u8 array; only present if any annotation type specifies P or Q)
- **AN:Z:** - Annotation Names (optional labels for individual annotations)

```
MA:Z:read_length;annotation_type1+P:start1-len1,start2-len2;annotation_type2-:start1-len1
AQ:B:C,qual1,qual2
AN:Z:name1,name2,name3
```

Regex for MA tag:

```
^\d+;(([a-zA-Z0-9_]+)[+-.][PQ]*:((\d+-\d+)(,\d+-\d+)*);?)+$
```

### Read Length

The first value in the MA tag is the length of the read at the time the annotation was made. This is important because:

- The read sequence may be modified after annotation (e.g., adapter trimming, quality trimming)
- Storing the original read length allows tools to detect and handle such modifications
- It enables validation that annotation coordinates are within bounds

### Molecular Coordinates

All coordinates in the MA tag are "molecular coordinates" meaning:

- **1-based positions**: Start positions use 1-based indexing (first base of the read is position 1), consistent with SAM format
- **Closed intervals**: The annotated region spans [start, start+length-1], where both start and end are inclusive
- **Read orientation**: Coordinates are always in the orientation of the sequenced molecule, counting from the left
- **Alignment independent**: For reverse-strand alignments, coordinates do not change; they reflect the original molecule orientation

### Delimiters

**MA tag delimiters:**

| Delimiter | Purpose                                      | Example                                 |
| --------- | -------------------------------------------- | --------------------------------------- |
| `;`       | Separates different annotation types         | `1000;msp+P:...;nuc-:...;fire.:...`     |
| `:`       | Separates annotation type name from data     | `msp+P:100-50,200-60`                   |
| `,`       | Separates annotations within a type          | `100-50,200-60,300-70`                  |
| `-`       | Separates start position from length         | `100-50` (start 100, length 50)         |

**AQ, AN tag delimiters:**

| Delimiter | Purpose          | Example    |
| --------- | ---------------- | ---------- |
| `,`       | Separates values | `40,35,200` |

### Annotation Type within the MA Tag

The annotation type name is an alphanumeric string (including underscores) that describes the type of annotation followed by a strand indicator and an optional quality type indicator:

**Strand indicator:**

- **`+`**: Annotation is on the forward strand of the sequenced molecule
- **`-`**: Annotation is on the reverse strand of the sequenced molecule
- **`.`**: Strand information is not applicable or unknown

**Quality type indicator:**

Zero or more characters specifying the number and type of quality values per annotation:

- **`P`**: One phred-scaled quality value
- **`Q`**: One linearly-scaled quality value (like ML tag)
- **Multiple characters** (e.g., `PQ`, `PQQP`): Each character adds one quality value per annotation. The character defines the scaling for that position.
- **(omitted)**: No quality scores for this annotation type

e.g., `msp+P` (one phred quality per annotation), `nuc-` (no quality), `fire.PQ` (two quality values per annotation, first phred-scaled, second linearly-scaled), `ctcf+PQQP` (four quality values per annotation)

When the quality type indicator is omitted, no quality values are stored for that annotation type in the AQ tag. The number of quality values per annotation equals the length of the quality type indicator string.

#### Strand Convention

Strand information describes the biology of an annotation feature, not the alignment orientation, and follows the same convention as the MM and ML tags for base modifications in the SAM specification:

- **Per-annotation property**: Strand is a property of each individual annotation, not of the annotation type. One annotation type may contain annotations on different strands.
- **All coordinates are in read orientation**: Coordinates always refer to positions on the forward strand of the original read sequence (i.e. starting from the leftmost base of the unaligned read), regardless of the strand a given annotation describes.
- **Alignment independent**: If a read aligns to the reverse strand of a reference, the MA tag strand indicators remain unchanged.

**Example of strand-specific annotations** on a read of length 10 (where `#` marks the annotated region, `-` marks unannotated positions, and coordinates are 1-based):

This shows CTCF annotations on both the forward strand (start 1 length 4) and reverse strand (start 6 length 3). They belong to the same logical `ctcf` annotation type — the on-disk format groups them into separate sections by strand.

```
MA:Z:10;ctcf+Q:1-4;ctcf-Q:6-3
AQ:B:C,200,180

Position: 1234567890
Forward:  ####------
Reverse:  -----###--
```

#### Annotation type uniqueness

Within a single MA tag, an annotation type is uniquely identified by its `name` alone. Strand is a property of individual annotations within the type, and `quality_spec` is a property of the type but is not part of its identity for uniqueness purposes (see below).

The on-disk MA tag format groups annotations by `(name, strand, quality_spec)` into sections for compactness. Producers MAY emit multiple sections sharing the same `name` when the annotations span different strands, e.g. `ctcf+Q:1-4;ctcf-Q:6-3` — these represent the same logical type with strand-split annotations.

Producers MUST NOT emit two sections with the same `name` and conflicting `quality_spec` within one MA tag. Parsers MUST treat repeated `name` sections as the same type and merge their annotations; parsers MUST error on conflicting `quality_spec` for the same `name`.

### Annotation Data

Each annotation is represented with the following components:

1. **MA tag** - Read length (first value) followed by annotations in `start-length` format (1-based coordinates)
2. **AQ tag** - Quality scores (0-255, u8; only present if any annotation type has a quality indicator)
3. **AN tag** - Optional name/label for the annotation (string)

The values in AN tags correspond positionally to all annotations defined in the MA tag. The AQ tag contains values only for annotation types that specify quality indicators. For each such type, each annotation contributes as many quality values as there are characters in the quality indicator. Values are grouped per-annotation: all quality values for one annotation appear consecutively, then all values for the next annotation, and so on, following MA tag order.

For example, if the MA tag contains `1000;msp+PQ:100-50,200-60;nuc+:150-103`, the AQ tag would contain 4 values: 2 per MSP annotation (phred then linear), and none for nuc (no quality indicator). The AQ array would be `[msp1_P, msp1_Q, msp2_P, msp2_Q]`.

**Note:** When using the AN tag, if some annotations have names and others don't, use an empty string to represent missing names. This maintains positional correspondence across all tags.

### Converting to Reference Coordinates

Reference coordinates are computed on-the-fly using the BAM alignment (CIGAR string) and are not stored in the `MA:Z` tag. i.e. alignment does not change the `MA` tag.

## Quality Scores

Quality scores (0-255) represent confidence in the annotation using the same convention as base modification quality scores.

## Common Annotation Types

### Standard Fiber-seq Annotations

| Type   | Description                                                         |
| ------ | ------------------------------------------------------------------- |
| `msp`  | Methylation-accessible patches (linker regions between nucleosomes) |
| `nuc`  | Nucleosome positions                                                |
| `fire` | Fiber-seq inferred regulatory elements (e.g., promoters, enhancers) |

### Custom Annotations

The format supports arbitrary annotation type names, allowing for:

- Custom chromatin features
- User-defined regions of interest
- Experimental annotations

## Examples

### Single Annotation Type (No Quality)

```
MA:Z:1000;msp+:100-50,200-60
```

- Read length: 1000 bp
- Two MSP annotations on the forward strand (no quality scores)
- First MSP: start position 100, length 50
- Second MSP: start position 200, length 60
- No AQ tag since no quality indicator was specified

### Single Annotation Type (Linear Quality)

```
MA:Z:1000;msp+Q:100-50,200-60
AQ:B:C,255,200
```

- Read length: 1000 bp
- Two MSP annotations on the forward strand with linear quality scores
- First MSP: start position 100, length 50, quality 255
- Second MSP: start position 200, length 60, quality 200

### Single Annotation Type with Phred Quality

```
MA:Z:1000;msp+P:100-50,200-60
AQ:B:C,40,30
```

- Read length: 1000 bp
- Two MSP annotations on the forward strand with phred-scaled quality scores
- First MSP: start position 100, length 50, phred quality 40
- Second MSP: start position 200, length 60, phred quality 30

### Multiple Annotation Types with Mixed Quality

```
MA:Z:1000;msp+P:100-50,200-60;nuc+:150-103,300-100;fire.Q:500-75
AQ:B:C,40,35,200
```

- Read length: 1000 bp
- 2 MSP annotations (forward strand, phred quality): positions 100 and 200, lengths 50 and 60, phred qualities 40 and 35
- 2 nucleosome annotations (forward strand, no quality): positions 150 and 300, lengths 103 and 100
- 1 FIRE annotation (no strand, linear quality): position 500, length 75, quality 200
- AQ tag contains 3 values: 2 for MSP (phred) + 1 for FIRE (linear); nucleosome annotations have no quality values

### Multiple Quality Values Per Annotation

```
MA:Z:1000;msp+PQ:100-50,200-60;nuc+:150-103
AQ:B:C,40,255,30,200
```

- Read length: 1000 bp
- 2 MSP annotations with quality spec `PQ` (2 quality values each: phred then linear)
- First MSP: start 100, length 50, phred quality 40, linear quality 255
- Second MSP: start 200, length 60, phred quality 30, linear quality 200
- 2 nucleosome annotations (no quality)
- AQ tag contains 4 values: `[msp1_P, msp1_Q, msp2_P, msp2_Q]`

### With Partial Names

```
MA:Z:1000;msp+P:100-50,200-60;nuc+:150-103,300-100
AQ:B:C,40,35
AN:Z:msp1,,,nuc2
```

- Read length: 1000 bp
- 2 MSP annotations with phred quality, 2 nucleosome annotations without quality
- AQ tag only contains 2 values (for the MSP annotations)
- Only the first MSP and second nucleosome have names
- Unnamed annotations use empty strings to maintain positional correspondence in AN

