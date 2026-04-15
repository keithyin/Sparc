# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build and Test

This is a Rust project that wraps a C++ consensus algorithm for long sequencing reads.

**Build commands:**
- `make` in `sparc-source-code/` compiles the C++ library (`libsparc.a`)
- `cargo build` or `cargo test` builds/tests the Rust bindings

**Test command:**
- `cargo test` - runs the consensus test in `src/lib.rs`

## Architecture

The project has two layers:

1. **C++ core** (`sparc-source-code/`): Implements the sparsity-based consensus algorithm with:
   - Graph construction from k-mers
   - Graph simplification
   - Best path finding
   - Main entry point: `SparcConsensus()` in `sparc.cpp`

2. **Rust bindings** (`src/lib.rs`): Provides safe Rust wrappers around the C++ library:
   - `SparcConfig` struct for algorithm parameters (k-mer size, coverage threshold, etc.)
   - `Query` struct representing aligned reads
   - `sparc_consensus()` function calling the C++ consensus engine
   - Memory-managed FFI types (`SparcQuery`, `SparcConsensusResult`)

## FFI Design

The Rust bindings use opaque C types (`CQuery`) to hide C++ implementation details. The `Query` struct in Rust fills a `CQuery` via FFI functions defined in `BasicDataStructure.h`.

Key C++ structures exposed to Rust: `Query`, `SparcConfig`, `SparcConsensusResult`.
