//! Small shared value types describing an annotation's strand, quality, and
//! on-disk encoding.
//!
//! - [`Strand`] - Strand orientation (+, -, .)
//! - [`QualityScaling`] - Scaling type for a single quality value (Phred or Linear)
//! - [`QualitySpec`] - Quality specification for an annotation type (zero or more scaling types)
//! - [`SkipFlag`] - MM "skip flag" (`.` / `?`)
//! - [`Encoding`] - How an `AnnotationType`'s data lives on disk

use std::fmt;
use std::str::FromStr;

use crate::ParseError;

/// Strand orientation of an annotation
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    /// Forward strand (+)
    Forward,
    /// Reverse strand (-)
    Reverse,
    /// Unknown or not applicable (.)
    Unknown,
}

impl Strand {
    /// Returns the character representation of the strand
    pub fn as_char(&self) -> char {
        match self {
            Strand::Forward => '+',
            Strand::Reverse => '-',
            Strand::Unknown => '.',
        }
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_char())
    }
}

impl FromStr for Strand {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" => Ok(Strand::Forward),
            "-" => Ok(Strand::Reverse),
            "." => Ok(Strand::Unknown),
            _ => Err(ParseError::InvalidFormat(format!("Invalid strand: {}", s))),
        }
    }
}

/// MM "skip flag" — the optional `.` or `?` after the mod code(s) in an MM group.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SkipFlag {
    #[default]
    Implicit,
    LowProbability,
    Unknown,
}

/// How an `AnnotationType`'s data lives on disk in the BAM record.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Encoding {
    #[default]
    Ma,
    MmMl {
        skip_flag: SkipFlag,
    },
}

impl Encoding {
    /// `MmMl` with the default (implicit) skip flag — the common case for
    /// producers synthesizing base modifications.
    pub fn mm_ml() -> Self {
        Encoding::MmMl {
            skip_flag: SkipFlag::Implicit,
        }
    }
}

/// Scaling type for a single quality value.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum QualityScaling {
    /// Phred-scaled quality score (P)
    Phred,
    /// Linear quality score 0-255 (Q)
    Linear,
}

impl QualityScaling {
    /// Returns the character representation
    pub fn as_char(&self) -> char {
        match self {
            QualityScaling::Phred => 'P',
            QualityScaling::Linear => 'Q',
        }
    }
}

impl fmt::Display for QualityScaling {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_char())
    }
}

/// Quality specification for an annotation type.
///
/// Defines the number and scaling type of quality values per annotation.
/// Empty means no quality values. The length determines how many quality
/// values each annotation in the type contributes to the AQ array.
///
/// # Examples
/// - `QualitySpec::none()` - no quality values
/// - `"P".parse::<QualitySpec>()` - one phred-scaled value per annotation
/// - `QualitySpec::from_str("PQ")` - two values per annotation (phred, then linear)
/// - `QualitySpec::from_str("PQQP")` - four values per annotation
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct QualitySpec(Vec<QualityScaling>);

impl QualitySpec {
    /// Create a quality spec from a list of scaling types.
    ///
    /// # Example
    /// ```
    /// use molecular_annotation::{QualitySpec, QualityScaling};
    /// let spec = QualitySpec::new(vec![
    ///     QualityScaling::Phred,
    ///     QualityScaling::Linear,
    ///     QualityScaling::Linear,
    ///     QualityScaling::Phred,
    /// ]);
    /// assert_eq!(spec.num_qualities(), 4);
    /// assert_eq!(spec.to_string(), "PQQP");
    /// ```
    pub fn new(scalings: Vec<QualityScaling>) -> Self {
        Self(scalings)
    }

    /// No quality values.
    pub fn none() -> Self {
        Self(Vec::new())
    }

    /// Number of quality values per annotation.
    pub fn num_qualities(&self) -> usize {
        self.0.len()
    }

    /// Returns true if this spec has any quality values.
    pub fn has_quality(&self) -> bool {
        !self.0.is_empty()
    }

    /// Access the individual scaling types.
    pub fn scalings(&self) -> &[QualityScaling] {
        &self.0
    }
}

impl fmt::Display for QualitySpec {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for s in &self.0 {
            write!(f, "{}", s.as_char())?;
        }
        Ok(())
    }
}

impl FromStr for QualitySpec {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut scalings = Vec::with_capacity(s.len());
        for c in s.chars() {
            match c {
                'P' => scalings.push(QualityScaling::Phred),
                'Q' => scalings.push(QualityScaling::Linear),
                _ => {
                    return Err(ParseError::InvalidFormat(format!(
                        "Invalid quality spec character '{}' in '{}'",
                        c, s
                    )))
                }
            }
        }
        Ok(Self(scalings))
    }
}
