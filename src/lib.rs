use std::ffi::{CString, c_char, c_int};

pub fn add(left: u64, right: u64) -> u64 {
    left + right
}

/// 对应 C++ 的 `struct Query`
/// Rust 不知道内部结构，只拿指针用
#[repr(C)]
pub struct CQuery {
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
        }
    }
}

unsafe extern "C" {

    pub unsafe fn SparcConsensus(
        backbone_c: *const c_char,
        queries: *mut *mut CQuery,
        n_queries: c_int,
        config: *const SparcConfig,
    ) -> *mut c_char;

    pub fn NewQuery() -> *mut CQuery;
    pub unsafe fn FreeQuery(query: *mut CQuery);

    /* ---------- string fields ---------- */
    pub unsafe fn QuerySetQueryName(query: *mut CQuery, name: *const c_char);
    pub unsafe fn QuerySetTargetName(query: *mut CQuery, name: *const c_char);
    pub unsafe fn QuerySetQueryAlignedSeq(query: *mut CQuery, seq: *const c_char);
    pub unsafe fn QuerySetMatchPattern(query: *mut CQuery, pattern: *const c_char);
    pub unsafe fn QuerySetTargetAlignedSeq(query: *mut CQuery, seq: *const c_char);

    /* ---------- basic attributes ---------- */
    pub unsafe fn QuerySetQueryLength(query: *mut CQuery, queryLength: c_int);
    pub unsafe fn QuerySetQueryStart(query: *mut CQuery, qStart: c_int);
    pub unsafe fn QuerySetQueryEnd(query: *mut CQuery, qEnd: c_int);

    pub unsafe fn QuerySetTargetLength(query: *mut CQuery, targetLength: c_int);
    pub unsafe fn QuerySetTargetStart(query: *mut CQuery, tStart: c_int);
    pub unsafe fn QuerySetTargetEnd(query: *mut CQuery, tEnd: c_int);

    /* ---------- alignment stats ---------- */
    pub unsafe fn QuerySetScore(query: *mut CQuery, score: c_int);
    pub unsafe fn QuerySetNumMatch(query: *mut CQuery, numMatch: c_int);
    pub unsafe fn QuerySetNumMismatch(query: *mut CQuery, numMismatch: c_int);
    pub unsafe fn QuerySetNumIns(query: *mut CQuery, numIns: c_int);
    pub unsafe fn QuerySetNumDel(query: *mut CQuery, numDel: c_int);
    pub unsafe fn QuerySetMapQV(query: *mut CQuery, mapQV: c_int);

    /* ---------- strand / index ---------- */
    pub unsafe fn QuerySetQueryStrand(query: *mut CQuery, strand: c_char);
    pub unsafe fn QuerySetTargetStrand(query: *mut CQuery, strand: c_char);
    pub unsafe fn QuerySetReadIndex(query: *mut CQuery, read_idx: usize);

    /* ---------- report range ---------- */
    pub unsafe fn QuerySetReportBegin(query: *mut CQuery, report_b: c_int);
    pub unsafe fn QuerySetReportEnd(query: *mut CQuery, report_e: c_int);

    /* ---------- counters ---------- */
    pub unsafe fn QuerySetNumExist(query: *mut CQuery, n_exist: c_int);
    pub unsafe fn QuerySetNumNew(query: *mut CQuery, n_new: c_int);

    /* ---------- patch flags ---------- */
    pub unsafe fn QuerySetPatch(query: *mut CQuery, patch: bool);
    pub unsafe fn QuerySetFill(query: *mut CQuery, fill: bool);

    /* ---------- patch parameters ---------- */
    pub unsafe fn QuerySetPatchK(query: *mut CQuery, k: c_int);
    pub unsafe fn QuerySetPatchD(query: *mut CQuery, d: c_int);
    pub unsafe fn QuerySetPatchG(query: *mut CQuery, g: c_int);

}

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
    pub fn fill_c_query(self, c_query: *mut CQuery) {
        unsafe {
            let query_aligned_seq = CString::new(self.query_aligned_seq).unwrap();
            QuerySetQueryAlignedSeq(c_query, query_aligned_seq.as_c_str().as_ptr());

            let target_aligned_seq = CString::new(self.target_aligned_seq).unwrap();
            QuerySetTargetAlignedSeq(c_query, target_aligned_seq.as_c_str().as_ptr());

            QuerySetQueryStrand(c_query, '+' as i8);
            QuerySetTargetStrand(c_query, '+' as i8);

            QuerySetTargetLength(c_query, (self.target_end - self.target_start) as c_int);
            QuerySetTargetStart(c_query, self.target_start as c_int);
            QuerySetTargetEnd(c_query, self.target_end as c_int);
        }
    }
}

#[cfg(test)]
mod tests {
    use std::ffi::{CStr, CString};

    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }

    #[test]
    fn test_sparc_consensus() {
        let backbone_str = CString::new("GATCGGGCTAA").unwrap();

        // backbone_str.as_c_str().as_ptr();

        let mut c_queries = vec![];
        unsafe {
            let c_query = NewQuery();

            let query = Query {
                query_aligned_seq: "GATCGCGCTAA".to_string(),
                target_aligned_seq: "GATCGGGCTAA".to_string(),
                rev_strand: false,
                query_start: 0,
                query_end: 11,
                target_start: 0,
                target_end: 11,
            };
            query.fill_c_query(c_query);
            c_queries.push(c_query);

            let c_query = NewQuery();

            let query = Query {
                query_aligned_seq: "GATCGCGCCAA".to_string(),
                target_aligned_seq: "GATCGGGCTAA".to_string(),
                rev_strand: false,
                query_start: 0,
                query_end: 11,
                target_start: 0,
                target_end: 11,
            };
            query.fill_c_query(c_query);
            c_queries.push(c_query);

            let query = Query {
                query_aligned_seq: "GCTCGGCCCAA".to_string(),
                target_aligned_seq: "GATCGGGCTAA".to_string(),
                rev_strand: false,
                query_start: 0,
                query_end: 11,
                target_start: 0,
                target_end: 11,
            };
            query.fill_c_query(c_query);
            c_queries.push(c_query);

            let query = Query {
                query_aligned_seq: "GCTCGGCCCAA".to_string(),
                target_aligned_seq: "GATCGGGCTAA".to_string(),
                rev_strand: false,
                query_start: 0,
                query_end: 11,
                target_start: 0,
                target_end: 11,
            };
            query.fill_c_query(c_query);
            c_queries.push(c_query);

            let query = Query {
                query_aligned_seq: "GCTCGGCCCAA".to_string(),
                target_aligned_seq: "GATCGGGCTAA".to_string(),
                rev_strand: false,
                query_start: 0,
                query_end: 11,
                target_start: 0,
                target_end: 11,
            };
            query.fill_c_query(c_query);
            c_queries.push(c_query);

            let query = Query {
                query_aligned_seq: "GATCGCGCCAA".to_string(),
                target_aligned_seq: "GATCGGGCTAA".to_string(),
                rev_strand: false,
                query_start: 0,
                query_end: 11,
                target_start: 0,
                target_end: 11,
            };
            query.fill_c_query(c_query);
            c_queries.push(c_query);

            let query = Query {
                query_aligned_seq: "GATCGCGCCAA".to_string(),
                target_aligned_seq: "GATCGGGCTAA".to_string(),
                rev_strand: false,
                query_start: 0,
                query_end: 11,
                target_start: 0,
                target_end: 11,
            };
            query.fill_c_query(c_query);
            c_queries.push(c_query);

        }

        let mut config = SparcConfig::default();
        config.debug = true;
        config.report_end = backbone_str.as_bytes().len() as c_int;
        config.subgraph_end = backbone_str.as_bytes().len() as c_int;
        config.cns_end = backbone_str.as_bytes().len() as c_int;
        unsafe {
            let c_str = SparcConsensus(
                backbone_str.as_ptr(),
                c_queries.as_mut_ptr(),
                c_queries.len() as c_int,
                &config as *const SparcConfig,
            );

            let c_str = CStr::from_ptr(c_str);
            println!("consensus:{}", c_str.to_str().unwrap());
        }
    }
}
