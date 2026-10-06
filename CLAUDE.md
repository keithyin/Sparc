# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build and Test

This is a Rust project that wraps a C++ consensus algorithm for long sequencing reads. Requires `g++` and `make`.

**Build commands:**
- `cargo build` / `cargo test` — `build.rs` copies `sparc-source-code/` into `OUT_DIR`, runs `make` there, and links `libsparc.a`
- `make` in `sparc-source-code/` compiles the C++ library standalone (`libsparc.a`); `make -C sparc-source-code clean` to reset

**Test command:**
- `cargo test` — unit tests in `src/lib.rs` (consensus behavior, input validation, m5 parsing, e2e against `testdata2/out.consensus.fasta`)

**Memory checking (leak regressions):**
```
cargo test --no-run
valgrind --leak-check=full --show-leak-kinds=definite \
  --errors-for-leak-kinds=definite --error-exitcode=99 \
  target/debug/deps/sparc-<hash>
```
CI runs the full test suite under valgrind; a `definitely lost` report is a failure. The C++ layer uses raw malloc/free paired with Rust-side `Drop` — any new code path that allocates `ConsensusNode`/`ConsensusEdgeNode` must be reachable from `SparcFreeInfo` (backbone right-subtrees **or** `Backbone::orphan_nodes`; orphan chains hang off nothing and rejoin backbone nodes via forward edges, so they must be freed **before** the backbone nodes).

## Architecture

The project has two layers:

1. **C++ core** (`sparc-source-code/`): Implements the sparsity-based consensus algorithm with:
   - Graph construction from k-mers (`GraphConstruction.cpp`: backbone k-mer chain + read branch paths)
   - Graph simplification / best path (`GraphSimplification.cpp`: `SparcFindBestPath`)
   - Library entry point: `SparcConsensus()` in `sparc.cpp` (refactored from `main.cpp`, which remains the original CLI)
   - C ABI wrappers for `Query` construction: `NewQuery`/`FreeQuery`/`QuerySet*` in `BasicDataStructure.cpp`

2. **Rust bindings** (`src/lib.rs`): Safe API around the C++ library:
   - `SparcConfig` — algorithm parameters; `#[repr(C)]` and must stay field-for-field identical to `struct SparcConfig` in `sparc.h` (change both sides together)
   - `Query` — one read alignment; build with `Query::new` / `Query::reverse_strand` / `Query::from_m5_row`, or parse an m5 file with `parse_m5`
   - `sparc_consensus(backbone, &queries, &config) -> Result<Consensus, SparcError>` — validates all inputs (k-mer range [1,16], backbone length/alphabet, alignment lengths, target span vs non-gap base count, coordinates) before crossing FFI, then returns the consensus sequence plus the best-path range on the backbone
   - Memory-managed FFI types: `SparcQuery` (per-query `FreeQuery` via `Drop`, panic-safe), `CConsensusResult` (frees the malloc'd sequence via `Drop`)

## FFI conventions

- Every C allocation crossing the boundary is owned by exactly one Rust type with a `Drop` impl (`SparcQuery`, `CConsensusResult`); never hand `raw pointers` to users.
- The C++ core must never unwind across the FFI boundary — keep input validation on the Rust side.
- Debug mode (`SparcConfig::debug = true`) makes the C++ layer write files (align_profile.txt, subgraph.dot, ...) into the process CWD; leave `debug` off in tests unless asserting on those artifacts.
