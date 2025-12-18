# Molecular annotation tags

With the advent of functional single-molecule sequencing, there is an emerging need to annotate non-genetic elements on individual DNA molecules. This is handled for base modifications using the MM and ML tags. However, we have a need for a more extensible format beyond base modifications that allows for generic annotation of any segment of a molecule. To accomplish this we are proposing the `MA` tag, and we describe the format of this tag below.

Integration of this spec or a similar one into the SAM/BAM/CRAM is needed to standardize the representation of molecular annotations across different sequencing platforms and analysis tools. Importantly, a standardized format will allow for integration in a wide range of genomic tools, e.g. visualization in genome browsers and allow developers to be justified in spending time integrating this format into their tools.

## Format Specification

### Structure

The molecular annotation format uses four related tags:

- **MA:Z:** - Read length followed by annotation start positions (required; first value is read length, then u32 arrays separated by annotation type prefixes)
- **AL:B:I** - Annotation Lengths (required, u32 array)
- **AQ:B:C** - Annotation Quality scores (optional, u8 array; only present if any annotation type specifies P or Q)
- **AN:Z:** - Annotation Names (optional labels for individual annotations)

```
MA:Z:read_length;annotation_type1+P:start1,start2;annotation_type2-:start1
AL:B:I,len1,len2,len3
AQ:B:C,qual1,qual2,qual3
AN:Z:name1,name2,name3
```

Regex for MA tag:

```
^\d+;(([a-zA-Z0-9_]+)[+-.][PQ]?:((\d+)(,\d+)*);?)+$
```

### Read Length

The first value in the MA tag is the length of the read at the time the annotation was made. This is important because:

- The read sequence may be modified after annotation (e.g., adapter trimming, quality trimming)
- Storing the original read length allows tools to detect and handle such modifications
- It enables validation that annotation coordinates are within bounds

### Molecular Coordinates

All coordinates in the MA tag are "molecular coordinates" meaning:

- **0-based positions**: Start positions use 0-based indexing (first base of the read is position 0)
- **Half-open intervals**: The annotated region spans [start, start+length), where start is inclusive and end is exclusive
- **Read orientation**: Coordinates are always in the orientation of the sequenced molecule, counting from the left
- **Alignment independent**: For reverse-strand alignments, coordinates do not change; they reflect the original molecule orientation

### Delimiters

**MA tag delimiters:**

| Delimiter | Purpose                                  | Example                             |
| --------- | ---------------------------------------- | ----------------------------------- |
| `;`       | Separates different annotation types     | `1000;msp+P:...;nuc-:...;fire.:...` |
| `:`       | Separates annotation type name from data | `msp+P:100,200`                     |
| `,`       | Separates start positions                | `100,200,300`                       |

**AL, AQ, AN tag delimiters:**

| Delimiter | Purpose          | Example    |
| --------- | ---------------- | ---------- |
| `,`       | Separates values | `50,60,70` |

### Annotation Type within the MA Tag

The annotation type name is an alphanumeric string (including underscores) that describes the type of annotation followed by a strand indicator and an optional quality type indicator:

**Strand indicator:**

- **`+`**: Annotation is on the forward strand of the sequenced molecule
- **`-`**: Annotation is on the reverse strand of the sequenced molecule
- **`.`**: Strand information is not applicable or unknown

**Quality type indicator (optional):**

- **`P`**: Quality scores are phred-scaled
- **`Q`**: Quality scores are linearly scaled (like ML tag)
- **(omitted)**: No quality scores for this annotation type

e.g., `msp+P`, `nuc-`, `fire.Q`

When the quality type indicator is omitted, no quality values are stored for that annotation type in the AQ tag. This means the AQ tag only contains values for annotations that have `P` or `Q` specified.

#### Strand Convention

Strand information is relative to the sequenced molecule and follows the same convention as the MM and ML tags for base modifications in the SAM specification:

- **All coordinates are in read orientation**: Coordinates always refer to positions on the forward strand of the original read sequence (i.e. starting from the leftmost base of the unaligned read).
- **For reverse-strand annotations (`-`)**: The annotation feature is on the reverse/Crick strand of the DNA molecule; however, coordinates still start from the left of the read sequence on the forward strand.
- **Strand is independent of alignment**: The strand indicator describes the biology of the feature, not the alignment orientation. If a read aligns to the reverse strand of a reference, the MA tag strand indicators remain unchanged.

**Example of strand-specific annotations** on a read of length 10 (where `#` marks the annotated region, `-` marks unannotated positions, and coordinates are 0-based):

This shows CTCF annotations on both the forward strand (start 0 length 4) and reverse strand (start 5 length 3).

```
MA:Z:10;ctcf+Q:0;ctcf-Q:5
AL:B:I,4,3
AQ:B:C,200,180

Position: 0123456789
Forward:  ####------
Reverse:  -----###--
```

### Annotation Data

Each annotation is represented across the four tags with corresponding values:

1. **MA tag** - Read length (first value) followed by start positions in molecular coordinates (0-based, u32 arrays with annotation type prefixes)
2. **AL tag** - Length of the annotation in base pairs (u32)
3. **AQ tag** - Quality score (0-255, u8; only for annotations with P or Q specified)
4. **AN tag** - Optional name/label for the annotation (string)

The values in AL and AN tags correspond positionally to all annotations defined in the MA tag. The AQ tag only contains values for annotations whose type specifies `P` or `Q`; annotations without a quality indicator are skipped. For example, if the MA tag contains `1000;msp+P:100,200;nuc+:150`, the AQ tag would contain only 2 values (for the two MSP annotations), while AL would contain 3 values (for all annotations).

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
MA:Z:1000;msp+:100,200
AL:B:I,50,60
```

- Read length: 1000 bp
- Two MSP annotations on the forward strand (no quality scores)
- First MSP: start position 100, length 50
- Second MSP: start position 200, length 60
- No AQ tag since no quality indicator was specified

### Single Annotation Type (Linear Quality)

```
MA:Z:1000;msp+Q:100,200
AL:B:I,50,60
AQ:B:C,255,200
```

- Read length: 1000 bp
- Two MSP annotations on the forward strand with linear quality scores
- First MSP: start position 100, length 50, quality 255
- Second MSP: start position 200, length 60, quality 200

### Single Annotation Type with Phred Quality

```
MA:Z:1000;msp+P:100,200
AL:B:I,50,60
AQ:B:C,40,30
```

- Read length: 1000 bp
- Two MSP annotations on the forward strand with phred-scaled quality scores
- First MSP: start position 100, length 50, phred quality 40
- Second MSP: start position 200, length 60, phred quality 30

### Multiple Annotation Types with Mixed Quality

```
MA:Z:1000;msp+P:100,200;nuc+:150,300;fire.Q:500
AL:B:I,50,60,103,100,75
AQ:B:C,40,35,200
```

- Read length: 1000 bp
- 2 MSP annotations (forward strand, phred quality): positions 100 and 200, lengths 50 and 60, phred qualities 40 and 35
- 2 nucleosome annotations (forward strand, no quality): positions 150 and 300, lengths 103 and 100
- 1 FIRE annotation (no strand, linear quality): position 500, length 75, quality 200
- AQ tag contains 3 values: 2 for MSP (phred) + 1 for FIRE (linear); nucleosome annotations have no quality values

### With Partial Names

```
MA:Z:1000;msp+P:100,200;nuc+:150,300
AL:B:I,50,60,103,100
AQ:B:C,40,35
AN:Z:msp1,,,nuc2
```

- Read length: 1000 bp
- 2 MSP annotations with phred quality, 2 nucleosome annotations without quality
- AQ tag only contains 2 values (for the MSP annotations)
- Only the first MSP and second nucleosome have names
- Unnamed annotations use empty strings to maintain positional correspondence in AN
