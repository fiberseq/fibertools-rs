//! Error type for parsing molecular annotations.

use std::fmt;

/// Error types for parsing molecular annotations
#[derive(Debug, Clone, PartialEq)]
pub enum ParseError {
    /// The MA tag format is invalid
    InvalidFormat(String),
    /// A coordinate value is invalid
    InvalidCoordinate(String),
    /// The number of values in AL/AQ/AN doesn't match MA
    MismatchedArrayLengths { expected: usize, got: usize },
    /// Annotation type already exists with different strand/quality_type
    ConflictingAnnotationType {
        /// The annotation type name that conflicts
        name: String,
        /// Description of the conflict
        reason: String,
    },
    /// Coordinate arithmetic would overflow
    CoordinateOverflow {
        /// The operation that caused overflow
        operation: String,
    },
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseError::InvalidFormat(msg) => write!(f, "Invalid format: {}", msg),
            ParseError::InvalidCoordinate(msg) => write!(f, "Invalid coordinate: {}", msg),
            ParseError::MismatchedArrayLengths { expected, got } => {
                write!(
                    f,
                    "Mismatched array lengths: expected {}, got {}",
                    expected, got
                )
            }
            ParseError::ConflictingAnnotationType { name, reason } => {
                write!(f, "Conflicting annotation type '{}': {}", name, reason)
            }
            ParseError::CoordinateOverflow { operation } => {
                write!(f, "Coordinate overflow in: {}", operation)
            }
        }
    }
}

impl std::error::Error for ParseError {}
