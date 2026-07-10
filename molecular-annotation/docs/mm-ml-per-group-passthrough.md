# Follow-up: per-group MM/ML passthrough

Status: proposed (follow-up to PR #107). Not blocking the MA migration merge.

## Problem

Base modifications on disk (`MM`/`ML`) must round-trip byte-identically. The
in-memory model normalizes MM/ML on read (splits grouped codes, canonicalizes
skip flags), so reconstructing MM/ML *from the model* is lossy for spec-legal
but non-canonical encodings (grouped multi-code `C+mh`, `N+a` wildcards).

Today we work around this with two write functions, and the caller picks based
on intent:

- `write_record` / `to_record` — writes only MA-family tags (`MA`/`AL`/`AQ`/`AN`),
  leaves the record's MM/ML bytes untouched. Used by structural editors
  (fire, add-nucs, footprint, pileup, extract, center, convert-tags) via
  `FiberseqData::serialize_annotations`.
- `write_record_with_basemods` / `write_mm_ml` — destructively re-encodes MM/ML
  canonically from the model. Used by base-mod producers (predict-m6a,
  ddda-to-m6a, strip-basemods, mock-fire, fibertig synthesis).

The split is a proxy for one boolean: "did this code path change base mods?"
It works, but the decision lives in the caller's choice of function — a footgun
(pick the wrong one → silently wrong MM/ML) and the source of the confusing
`_with_basemods` naming.

## Proposal: track dirtiness per original MM group, re-encode only what changed

Re-encode *only* the base-mod groups that were actually added or modified; pass
every untouched group through byte-identically. This collapses the two write
functions into one `write_record` whose behavior is derived from the data, not
declared by the caller.

### The tracking unit is the MM *group*, not the model *type*

`parse_into` splits a multi-code group like `C+mh` into two separate annotation
types (`m` and `h`) that share positions with interleaved ML, and currently
keeps no record of which group a type came from. The unit that can be preserved
byte-identically is the original MM group substring (`C+mh,5,12;`) plus its ML
byte slice — not an individual model type.

### Mechanics

- On read (`parse_into`): retain, per original MM group, its exact substring +
  its ML byte slice + the set of model types it produced. (New machinery — this
  provenance is discarded today.)
- Dirty rule: a group is *clean* iff **every** model type derived from it is
  unmodified. Clean → emit retained original bytes. Dirty or newly-created →
  reconstruct canonically from the model.
- On write: emit clean groups first, in original recorded order with their
  original ML slices, then append re-encoded / new groups. If nothing is dirty,
  output == input exactly.

### Why this beats the current whole-record split

Shrinks the re-encode blast radius to exactly the modification that changed.
Concrete win (the case @mrvollger flagged): a BAM with `MM:Z:C+mh,…;A+a,…`
(5mC + m6a). Run predict-m6a, which regenerates only `A+a`.

- Whole-record re-encode: the untouched `C+mh` gets canonicalized to
  `C+m;C+h` → 5mC bytes collaterally damaged.
- Per-group: `A+a` re-encoded, `C+mh` spliced back verbatim → 5mC survives
  byte-identical.

It also makes non-canonical encodings byte-identical as long as nothing touches
them — strictly stronger than #107's current "canonical only" guarantee.

### Sharp edge (document, rarely bites)

A multi-code group's codes share one MM group and interleaved ML, so you can't
preserve one code while re-encoding its sibling. Modifying `C+h` but not `C+m`
(both from `C+mh`) forces re-encoding the whole group, which comes back as split
form `C+m;C+h`, not `C+mh`. Rule: preserve a group iff *all* its codes are clean.

In practice this ~never triggers: real fiberseq data keeps m6a (`A+a`) and 5mC
(`C+m`) as separate groups, and grouped `C+mh` isn't emitted in the wild.

## Scope / where it lives

Library-side change in `molecular-annotation`:

- `basemods/parse.rs` — retain per-group provenance on read.
- `serialize.rs` — splice clean (original bytes) vs. re-encoded (from model)
  groups; unify into a single write entry point.
- Callers in fibertools-rs collapse to one `write_record`; drop
  `write_record_with_basemods`.

### Tests to add

- Mixed `A+a` + `C+mh`: modify only the m6a group, assert the 5mC group is
  byte-identical.
- Sharp edge: modify one code of a grouped pair, assert documented split
  behavior.
- Structural-only edit (fire/add-nucs): assert MM/ML byte-identical (regression
  guard for the current guarantee).
- Producer (predict-m6a / strip-basemods): assert re-encode / removal happens.

## Interim (in #107)

Keep the two functions but rename to state intent and document the "which do I
call?" contract in one place. Behavior-preserving; safe to land with the merge.
