//! Type definitions for molecular annotations.
//!
//! - [`error`] - [`ParseError`]
//! - [`primitives`] - small shared value types ([`Strand`], [`QualityScaling`],
//!   [`QualitySpec`], [`MaEncoding`], [`SkipFlag`], [`Encoding`])
//! - [`annotation`] - [`Annotation`] and the [`Qualities`] alias
//! - [`annotation_type`] - [`AnnotationType`] plus its tag-emission helpers and
//!   the [`MaParts`] / [`MmGroup`] / [`MmMlParts`] output structs
//! - [`view`] - read-side views ([`AnnotationInfo`], [`ProjectedAnnotation`],
//!   [`LiftedCoords`])

mod annotation;
mod annotation_type;
mod error;
mod primitives;
mod view;

pub use annotation::{Annotation, Qualities};
pub use annotation_type::{AnnotationType, MaParts, MmGroup, MmMlParts};
pub use error::ParseError;
pub use primitives::{Encoding, MaEncoding, QualityScaling, QualitySpec, SkipFlag, Strand};
pub use view::{AnnotationInfo, LiftedCoords, ProjectedAnnotation};
