use std::ffi::{CStr, CString, c_char, c_int};

/// 对应 C++ 的 `struct Query`
/// Rust 不知道内部结构，只拿指针用
#[repr(C)]
struct CQuery {
    _private: [u8; 0],
}

#[repr(C)]
pub struct SparcConfig {
    pub debug: bool,
    pub kmer: c_int,
    pub converage_threshold: c_int,
    pub scoring_method: c_int,
    pub subgraph_begin: c_int,
    pub subgraph_end: c_int,
    pub cns_start: c_int,
    pub cns_end: c_int,
    pub report_begin: c_int,
    pub report_end: c_int,
    pub cov_radius: c_int,
}
impl Default for SparcConfig {
    fn default() -> Self {
        Self {
            debug: false,
            kmer: 1,
            converage_threshold: 2,
            scoring_method: 2,
            subgraph_begin: 0,
            subgraph_end: 0,
            cns_start: 0,
            cns_end: 0,
            report_begin: 0,
            report_end: 0,
            cov_radius: 2,
        }
    }
}

#[repr(C)]
pub struct SparcConsensusResult {
    seq: *mut u8,
}

impl SparcConsensusResult {
    pub fn into_string(&self) -> String {
        unsafe {
            CStr::from_ptr(self.seq as *const i8)
                .to_str()
                .unwrap()
                .to_string()
        }
    }
}

impl Drop for SparcConsensusResult {
    fn drop(&mut self) {
        unsafe {
            SparcFreeConsensusResult(self.seq);
        }
    }
}

#[allow(unused)]
unsafe extern "C" {

    unsafe fn SparcConsensus(
        backbone_c: *const c_char,
        queries: *mut *mut CQuery,
        n_queries: c_int,
        config: *const SparcConfig,
    ) -> SparcConsensusResult;

    unsafe fn NewQuery() -> *mut CQuery;
    unsafe fn FreeQuery(query: *mut CQuery);

    /* ---------- string fields ---------- */
    unsafe fn QuerySetQueryName(query: *mut CQuery, name: *const c_char);
    unsafe fn QuerySetTargetName(query: *mut CQuery, name: *const c_char);
    unsafe fn QuerySetQueryAlignedSeq(query: *mut CQuery, seq: *const c_char);
    unsafe fn QuerySetMatchPattern(query: *mut CQuery, pattern: *const c_char);
    unsafe fn QuerySetTargetAlignedSeq(query: *mut CQuery, seq: *const c_char);

    /* ---------- basic attributes ---------- */
    unsafe fn QuerySetQueryLength(query: *mut CQuery, queryLength: c_int);
    unsafe fn QuerySetQueryStart(query: *mut CQuery, qStart: c_int);
    unsafe fn QuerySetQueryEnd(query: *mut CQuery, qEnd: c_int);

    unsafe fn QuerySetTargetLength(query: *mut CQuery, targetLength: c_int);
    unsafe fn QuerySetTargetStart(query: *mut CQuery, tStart: c_int);
    unsafe fn QuerySetTargetEnd(query: *mut CQuery, tEnd: c_int);

    /* ---------- alignment stats ---------- */
    unsafe fn QuerySetScore(query: *mut CQuery, score: c_int);
    unsafe fn QuerySetNumMatch(query: *mut CQuery, numMatch: c_int);
    unsafe fn QuerySetNumMismatch(query: *mut CQuery, numMismatch: c_int);
    unsafe fn QuerySetNumIns(query: *mut CQuery, numIns: c_int);
    unsafe fn QuerySetNumDel(query: *mut CQuery, numDel: c_int);
    unsafe fn QuerySetMapQV(query: *mut CQuery, mapQV: c_int);

    /* ---------- strand / index ---------- */
    unsafe fn QuerySetQueryStrand(query: *mut CQuery, strand: c_char);
    unsafe fn QuerySetTargetStrand(query: *mut CQuery, strand: c_char);
    unsafe fn QuerySetReadIndex(query: *mut CQuery, read_idx: usize);

    /* ---------- report range ---------- */
    unsafe fn QuerySetReportBegin(query: *mut CQuery, report_b: c_int);
    unsafe fn QuerySetReportEnd(query: *mut CQuery, report_e: c_int);

    /* ---------- counters ---------- */
    unsafe fn QuerySetNumExist(query: *mut CQuery, n_exist: c_int);
    unsafe fn QuerySetNumNew(query: *mut CQuery, n_new: c_int);

    /* ---------- patch flags ---------- */
    unsafe fn QuerySetPatch(query: *mut CQuery, patch: bool);
    unsafe fn QuerySetFill(query: *mut CQuery, fill: bool);

    /* ---------- patch parameters ---------- */
    unsafe fn QuerySetPatchK(query: *mut CQuery, k: c_int);
    unsafe fn QuerySetPatchD(query: *mut CQuery, d: c_int);
    unsafe fn QuerySetPatchG(query: *mut CQuery, g: c_int);

    unsafe fn SparcFreeConsensusResult(seq: *mut u8);

}

#[allow(unused)]
pub struct Query {
    query_aligned_seq: String,
    target_aligned_seq: String,
    rev_strand: bool,
    query_start: usize,
    query_end: usize,
    target_start: usize,
    target_end: usize,
}

impl Query {
    pub fn new(
        target_aligned: String,
        query_aligned: String,
        target_start: usize,
        target_end: usize,
    ) -> Self {
        Self {
            query_aligned_seq: query_aligned,
            target_aligned_seq: target_aligned,
            rev_strand: false,
            query_start: 0,
            query_end: 0,
            target_start,
            target_end,
        }
    }

    fn fill_c_query(&self, c_query: *mut CQuery) {
        unsafe {
            let query_aligned_seq = CString::new(self.query_aligned_seq.as_bytes()).unwrap();
            QuerySetQueryAlignedSeq(c_query, query_aligned_seq.as_c_str().as_ptr());

            let target_aligned_seq = CString::new(self.target_aligned_seq.as_bytes()).unwrap();
            QuerySetTargetAlignedSeq(c_query, target_aligned_seq.as_c_str().as_ptr());

            QuerySetQueryStrand(c_query, '+' as i8);
            QuerySetTargetStrand(c_query, '+' as i8);

            QuerySetTargetLength(c_query, (self.target_end - self.target_start) as c_int);
            QuerySetTargetStart(c_query, self.target_start as c_int);
            QuerySetTargetEnd(c_query, self.target_end as c_int);
        }
    }
}

struct SparcQuery {
    c_query: *mut CQuery,
}

impl From<&Query> for SparcQuery {
    fn from(value: &Query) -> Self {
        let c_query = unsafe {
            let c_query = NewQuery();
            value.fill_c_query(c_query);
            c_query
        };
        Self { c_query }
    }
}

impl Drop for SparcQuery {
    fn drop(&mut self) {
        unsafe {
            FreeQuery(self.c_query);
        }
    }
}

pub fn sparc_consensus(backbone: &str, queries: &[Query], config: &SparcConfig) -> String {
    let queries = queries
        .iter()
        .map(|v| v.into())
        .collect::<Vec<SparcQuery>>();
    let mut c_queries = queries
        .iter()
        .map(|v| v.c_query)
        .collect::<Vec<*mut CQuery>>();

    let backbone_str = CString::new(backbone).unwrap();
    let result = unsafe {
        SparcConsensus(
            backbone_str.as_ptr(),
            c_queries.as_mut_ptr(),
            c_queries.len() as c_int,
            config as *const SparcConfig,
        )
    };

    result.into_string()
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_sparc_consensus() {
        let backbone = "GATCGGGCTAA";

        let mut config = SparcConfig::default();
        config.debug = false;
        config.report_end = backbone.as_bytes().len() as c_int;
        config.subgraph_end = backbone.as_bytes().len() as c_int;
        config.cns_end = backbone.as_bytes().len() as c_int;

        let queries = vec![
            Query {
                query_aligned_seq: "GATCGCGCTAA".to_string(),
                target_aligned_seq: "GATCGGGCTAA".to_string(),
                rev_strand: false,
                query_start: 0,
                query_end: 11,
                target_start: 0,
                target_end: 11,
            },
            Query {
                query_aligned_seq: "GATCGCGCCAA".to_string(),
                target_aligned_seq: "GATCGGGCTAA".to_string(),
                rev_strand: false,
                query_start: 0,
                query_end: 11,
                target_start: 0,
                target_end: 11,
            },
            Query {
                query_aligned_seq: "GCTCGGCCCAA".to_string(),
                target_aligned_seq: "GATCGGGCTAA".to_string(),
                rev_strand: false,
                query_start: 0,
                query_end: 11,
                target_start: 0,
                target_end: 11,
            },
            Query {
                query_aligned_seq: "GCTCGGCCCAA".to_string(),
                target_aligned_seq: "GATCGGGCTAA".to_string(),
                rev_strand: false,
                query_start: 0,
                query_end: 11,
                target_start: 0,
                target_end: 11,
            },
            Query {
                query_aligned_seq: "GCTCGGCCCAA".to_string(),
                target_aligned_seq: "GATCGGGCTAA".to_string(),
                rev_strand: false,
                query_start: 0,
                query_end: 11,
                target_start: 0,
                target_end: 11,
            },
            Query {
                query_aligned_seq: "GATCGCGCCAA".to_string(),
                target_aligned_seq: "GATCGGGCTAA".to_string(),
                rev_strand: false,
                query_start: 0,
                query_end: 11,
                target_start: 0,
                target_end: 11,
            },
            Query {
                query_aligned_seq: "GATCGCGCCAA".to_string(),
                target_aligned_seq: "GATCGGGCTAA".to_string(),
                rev_strand: false,
                query_start: 0,
                query_end: 11,
                target_start: 0,
                target_end: 11,
            },
        ];

        let seq = sparc_consensus(backbone, &queries, &config);
        println!("consensus_seq:{seq}");
    }
}
