# Sparc: a sparsity-based consensus algorithm for long erroneous sequencing reads

rust binding of consensus

Ye C, Ma Z. (2016) Sparc: a sparsity-based consensus algorithm for long erroneous sequencing reads. PeerJ 4:e2016 https://doi.org/10.7717/peerj.2016  

test data:
https://sourceforge.net/projects/sparc-consensus/files/testdata/


## Rust binding

Requires `g++` and `make` (Linux/macOS). The C++ core is compiled by `build.rs`.

```rust
use sparc::{parse_m5, sparc_consensus, Query, SparcConfig};
use std::io::BufReader;
use std::ffi::c_int;

// backbone: reference-like sequence to polish
let backbone = "GATCGGGCTAA";

// each Query is one read alignment against the backbone
// (m5-style aligned strings with '-' gaps)
let queries = vec![
    Query::new(backbone.to_string(), "GATCGCGCTAA".to_string(), 0, backbone.len()),
    Query::new(backbone.to_string(), "GATCGCGCCAA".to_string(), 0, backbone.len()),
    Query::new(backbone.to_string(), "GCTCGGCCCAA".to_string(), 0, backbone.len()),
];

let mut config = SparcConfig::default();
config.kmer = 1;          // k-mer size, 1..=16 (suggested [1, 2])
config.coverage_threshold = 2;  // CLI "c", suggested [1, 5]
config.threshold = -0.1;  // CLI "t", adaptive threshold (<0 disables)

let consensus = sparc_consensus(backbone, &queries, &config).unwrap();
println!("consensus: {}", consensus.seq);
// consensus.start / consensus.end: best-path range on the backbone [start, end);
// end == 0 means "no trusted path, seq is the backbone as-is";
// start == None means the path head is not on the backbone.

// alignments can also be parsed from blasr m5 rows/files:
let m5 = std::fs::File::open("backbone-0.mapped.m5").unwrap();
let queries = parse_m5(BufReader::new(m5)).unwrap();
// negative-strand rows (tStrand '-') are handled automatically
```

All inputs are validated before entering the C++ core; invalid input (k-mer
out of range, backbone shorter than k, non-ACGT bases, inconsistent alignment
lengths/spans, out-of-range coordinates) returns a `SparcError` instead of
crashing.

## Original CLI

To compile the original command line tool from scratch, clone the directory and use the following command:

g++ -O3 -o Sparc *.cpp


Parameters:

b: backbone file.

m: the reads mapping files produced by blasr, using option -m 5. (A blasr example command: blasr -nproc 32 query.reads.fasta backbone.fasta -bestn 1 -m 5 -minMatch 19 -out backbone.mapped.m5)

k: k-mer size (suggested range: [1,2]).

c: coverage threshold (range: [1,5], suggest: 2).

t: adaptive threshold (suggested range[0.0,0.3]).

g: skip size, the larger the value, the more memory efficient the algorithm is (suggested range: [1,3]).

HQ_Prefix: Shared prefix of the high quality read names. (e.g. if the sec-gen sequences have names >Contig_xxx, then ‘Contig’ is a shared prefix of the high quality reads)

boost: boosting weight for the high quality reads (suggested range: [1,10]).

Example command: 

Using third-gen data only:

Sparc b Backbone.fa m backbone.mapped.m5 k 2 g 2 c 2 t 0.1 o ConsensusOutput

Using hybrid data:

Sparc b backbone.fasta m backbone.mapped.m5 k 2 g 2 c 2 t 0.1 HQ_Prefix Contig boost 5 o ConsensusOutput



