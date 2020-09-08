#![crate_name = "squiggler"]
#![crate_type = "lib"]
//#![warn(missing_docs)]

//!A tiny implementation of convertion from fasta to squiggle for Rust language.
//!
extern crate rand;
mod squiggler;
pub use squiggler::*;
mod test;
