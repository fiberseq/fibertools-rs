# molecular-annotation

Python bindings for molecular annotation tags in SAM/BAM/CRAM files.

This package provides a Python interface to parse and generate MA/AL/AQ/AN tags according to the Molecular Annotation specification.

## Installation

```bash
pip install molecular-annotation
```

## Usage

```python
from molecular_annotation import MolecularAnnotations

# Parse from inline format (start-length pairs in MA string)
# Note: if annotation type has quality (P or Q), must provide aq array
ma_string = "1000;msp+P:100-50,200-60"
aq_array = [40, 35]  # Quality scores for the 2 msp annotations
annotations = MolecularAnnotations.from_tags(ma_string, [], aq=aq_array)

# Parse from inline format without quality scores
ma_string = "1000;nuc+:100-147,250-147"
annotations = MolecularAnnotations.from_tags(ma_string, [])

# Parse from separate format (lengths in AL array)
ma_string = "1000;msp+P:100,200"
al_array = [50, 60]
aq_array = [40, 35]  # Quality scores for msp annotations
annotations = MolecularAnnotations.from_tags(ma_string, al_array, aq=aq_array)

# Access properties
print(annotations.read_length)  # 1000
print(annotations.total_annotation_count())  # 2

# Generate tag values
print(annotations.to_ma_string())  # "1000;msp+P:100,200"
print(annotations.to_al_array())   # [50, 60]
print(annotations.to_aq_array())   # [40, 35]
print(annotations.to_an_string())  # None (no names)
```

## Building Annotations

```python
from molecular_annotation import MolecularAnnotations

# Create empty container
annotations = MolecularAnnotations(1000)  # read length

# Add annotations in batches (vectorized)
annotations.add_annotations(
    'msp', '+', 'P',            # type_name, strand, quality_type
    starts=[100, 200, 350],     # 1-based positions
    lengths=[50, 60, 45],       # annotation lengths
    qualities=[40, 35, 38]      # quality scores (required for 'P' types)
)

# Add more without quality
annotations.add_annotations(
    'nuc', '+', '',             # no quality type
    starts=[150, 400],
    lengths=[147, 147]
)

print(annotations.to_ma_string())
# "1000;msp+P:100-50,200-60,350-45;nuc+:150-147,400-147"
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
annotations.add_annotations('msp', '+', 'P', [100, 200], [50, 60], qualities=[40, 35])
write_to_record(annotations, record)
```

## Tag Format

See the [Molecular Annotation Specification](../README.md) for details on the tag formats and conventions.