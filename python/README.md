# molecular-annotation

Python bindings for molecular annotation tags in SAM/BAM/CRAM files.

This package provides a Python interface to parse and generate MA/AQ/AN tags according to the Molecular Annotation specification.

## Installation

```bash
pip install molecular-annotation
```

## Usage

```python
from molecular_annotation import MolecularAnnotations

# Parse from MA tag string (positions in the tag are 1-based per spec)
# If annotation type has quality (P, Q, PQ, etc.), must provide aq array
ma_string = "1000;msp+P:100-50,200-60"
aq_array = [40, 35]  # Quality scores for the 2 msp annotations
annotations = MolecularAnnotations.from_tags(ma_string, aq=aq_array)

# Parse without quality scores
ma_string = "1000;nuc+:100-147,250-147"
annotations = MolecularAnnotations.from_tags(ma_string)

# Access properties
print(annotations.read_length)  # 1000
print(annotations.total_annotation_count())  # 2

# Generate tag values
print(annotations.to_ma_string())  # "1000;nuc+:100-147,250-147"
print(annotations.to_aq_array())   # None (no quality for nuc type)
print(annotations.to_an_string())  # None (no names)
```

## Building Annotations

All API coordinates are **0-based half-open** `[start, end)`. The MA tag string uses 1-based positions per spec, and conversion is handled automatically.

```python
from molecular_annotation import MolecularAnnotations

# Create empty container
annotations = MolecularAnnotations(1000)  # read length

# Add annotations with a single quality value per annotation
annotations.add_annotations(
    'msp', '+', 'P',            # type_name, strand, quality_spec
    starts=[100, 200, 350],     # 0-based start positions
    lengths=[50, 60, 45],       # annotation lengths
    qualities=[40, 35, 38]      # one quality score per annotation
)

# Add annotations without quality
annotations.add_annotations(
    'nuc', '+', '',             # empty quality_spec = no quality values
    starts=[150, 400],
    lengths=[147, 147]
)

# Add annotations with multiple quality values per annotation
# "PQ" = 2 values per annotation (first phred-scaled, second linear-scaled)
annotations.add_annotations(
    'ctcf', '+', 'PQ',
    starts=[500, 700],
    lengths=[20, 30],
    qualities=[40, 255, 30, 200]  # 2 annotations x 2 values = 4 total
)

print(annotations.to_ma_string())
# "1000;msp+P:101-50,201-60,351-45;nuc+:151-147,401-147;ctcf+PQ:501-20,701-30"
```

You can also use `ends` instead of `lengths`:

```python
annotations.add_annotations(
    'msp', '+', 'P',
    starts=[100, 200],
    ends=[150, 260],            # 0-based exclusive end positions
    qualities=[40, 35]
)
```

## Iterating Over Annotations

```python
# Iterate over all annotations with full coordinate info
for type_name, strand, quality_spec, qs, qe, fs, fe, rs, re, quals, name in annotations.iter_full():
    print(f"{type_name} [{qs}, {qe}) quals={quals}")

# Iterate over a specific type
for qs, qe, fs, fe, rs, re, quals, name in annotations.iter_type("msp"):
    print(f"[{qs}, {qe}) quals={quals}")
```

## Working with pysam

```python
import pysam
from molecular_annotation import MolecularAnnotations, from_record, write_to_record

# Read annotations from a BAM record
with pysam.AlignmentFile("input.bam") as bam:
    for record in bam:
        try:
            annotations = from_record(record)
            print(f"{record.query_name}: {annotations.total_annotation_count()} annotations")
        except KeyError:
            pass  # No MA tag

# Write annotations to a BAM record
annotations = MolecularAnnotations(1000)
annotations.add_annotations('msp', '+', 'P', [100, 200], lengths=[50, 60], qualities=[40, 35])
write_to_record(annotations, record)
```

## Tag Format

See the [Molecular Annotation Specification](../README.md) for details on the tag formats and conventions.
