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
    start_pos: c_int,
    end_pos: c_int,
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
        // fill_c_query 在序列含内嵌 NUL 时会 panic，此时 guard 保证已创建的
        // C++ Query 仍被释放；成功路径 forget 掉 guard，交由 SparcQuery 的 Drop 管理
        struct FreeOnPanic(*mut CQuery);
        impl Drop for FreeOnPanic {
            fn drop(&mut self) {
                unsafe { FreeQuery(self.0) };
            }
        }

        let c_query = unsafe { NewQuery() };
        let guard = FreeOnPanic(c_query);
        value.fill_c_query(c_query);
        std::mem::forget(guard);
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

/// sparc_consensus
/// return: (cons_seq, start_in_backbone, end_in_backbone)
/// [start_in_backbone, end_in_backbone)
pub fn sparc_consensus(
    backbone: &str,
    queries: &[Query],
    config: &SparcConfig,
) -> (String, i32, i32) {
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

    let cons_seq = result.into_string();
    (cons_seq, result.start_pos, result.end_pos)
}

#[cfg(test)]
mod tests {

    use super::*;

    fn make_config(backbone_len: usize) -> SparcConfig {
        let mut config = SparcConfig::default();
        config.debug = false;
        config.report_end = backbone_len as c_int;
        config.subgraph_end = backbone_len as c_int;
        config.cns_end = backbone_len as c_int;
        config
    }

    fn make_queries(backbone: &str) -> Vec<Query> {
        let q = |query_aligned_seq: &str| Query {
            query_aligned_seq: query_aligned_seq.to_string(),
            target_aligned_seq: backbone.to_string(),
            rev_strand: false,
            query_start: 0,
            query_end: query_aligned_seq.len(),
            target_start: 0,
            target_end: backbone.len(),
        };
        vec![
            q("GATCGCGCTAA"),
            q("GATCGCGCCAA"),
            q("GCTCGGCCCAA"),
            q("GCTCGGCCCAA"),
            q("GCTCGGCCCAA"),
            q("GATCGCGCCAA"),
            q("GATCGCGCCAA"),
        ]
    }

    #[test]
    fn test_sparc_consensus() {
        let backbone = "GATCGGGCTAA";
        let config = make_config(backbone.len());
        let queries = make_queries(backbone);

        let seq = sparc_consensus(backbone, &queries, &config);
        println!("consensus_seq:{seq:?}");
    }

    /// 回归：max_score == 0 的早退路径必须先释放 k-mer 图和 ref.read_bits。
    /// 此前该路径在 SparcFreeInfo/free 之前 return，空 queries 即可触发整体泄漏。
    #[test]
    fn test_empty_queries_returns_backbone() {
        let backbone = "GATCGGGCTAA";
        let config = make_config(backbone.len());

        let (seq, start, end) = sparc_consensus(backbone, &[], &config);
        assert_eq!(seq, backbone);
        assert_eq!((start, end), (0, 0));
    }

    /// 回归：SparcFreeInfo 须覆盖最后一个 backbone 节点（每次调用都会走到该路径），
    /// 同时验证跨调用的 FFI 配对（NewQuery/FreeQuery、malloc/free）不累积泄漏。
    #[test]
    fn test_sparc_consensus_repeated() {
        let backbone = "GATCGGGCTAA";
        let config = make_config(backbone.len());
        let queries = make_queries(backbone);

        for _ in 0..50 {
            let (seq, start, end) = sparc_consensus(backbone, &queries, &config);
            assert!(!seq.is_empty());
            assert!(start < end);
        }
    }

    /// 回归：3' 端插入使非 backbone 节点挂在最后一个 backbone 节点上，
    /// SparcFreeInfo 现在会释放该子树；若释放逻辑存在双重释放会直接 abort。
    #[test]
    fn test_consensus_with_terminal_insertion() {
        let backbone = "GATCGGGCTAA";
        let config = make_config(backbone.len());
        let mut queries = make_queries(backbone);

        // tAligned 末尾为 '-'，即在 backbone 末端插入一个碱基
        let ins = |query_aligned_seq: &str| Query {
            query_aligned_seq: query_aligned_seq.to_string(),
            target_aligned_seq: format!("{backbone}-"),
            rev_strand: false,
            query_start: 0,
            query_end: query_aligned_seq.len(),
            target_start: 0,
            target_end: backbone.len(),
        };
        queries.push(ins("GATCGGGCTAAA"));
        queries.push(ins("GATCGGGCTAAG"));

        let (seq, start, end) = sparc_consensus(backbone, &queries, &config);
        assert!(!seq.is_empty());
        assert!(start < end);
    }
}
