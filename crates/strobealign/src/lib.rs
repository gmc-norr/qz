#![allow(warnings, clippy::all)]
// Slimmed to the MAP-ONLY core qz uses (index build + seeding + chaining → NAMs).
// The upstream alignment/SAM stack (mapper, aligner, piecewisealigner, ssw, cigar,
// io::sam) and the CLI read I/O (io::reads, io::fastq, io::paf) were removed — qz
// does ungapped-Hamming placement itself and feeds raw read bytes directly. See the
// crate README / git history for the prune.
pub mod chainer;
pub mod details;
pub mod hit;
pub mod index;
pub mod indexer;
pub mod insertsize;
pub mod io;
pub mod maponly;
pub mod math;
pub mod mcsstrategy;
pub mod nam;
pub mod partition;
pub mod read;
pub mod revcomp;
pub mod seeding;
pub mod shuffle;
