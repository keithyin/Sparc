//! Safe Rust bindings for [Sparc](https://github.com/keithyin/Sparc), a
//! sparsity-based consensus algorithm for long erroneous sequencing reads
//! (Ye C, Ma Z. 2016, PeerJ 4:e2016).
//!
//! The C++ core in `sparc-source-code/` is compiled and linked by `build.rs`.
//! All interaction goes through [`sparc_consensus`]: give it a backbone
//! sequence plus alignments ([`Query`]) of reads against that backbone, and it
//! returns the polished consensus.
//!
//! # Example
//!
//! ```
//! use sparc::{Query, SparcConfig, sparc_consensus};
//!
//! let backbone = "GATCGGGCTAA";
//! let queries = ["GATCGCGCTAA", "GATCGCGCCAA", "GATCGCGCCAA", "GATCGCGCCAA"]
//!     .iter()
//!     .map(|seq| Query::new(backbone.to_string(), seq.to_string(), 0, backbone.len()))
//!     .collect::<Vec<_>>();
//!
//! let consensus = sparc_consensus(backbone, &queries, &SparcConfig::default()).unwrap();
//! assert!(!consensus.seq.is_empty());
//! ```
//!
//! Inputs are validated before crossing the FFI boundary; invalid input
//! returns [`SparcError`] instead of crashing inside the C++ core.

use std::ffi::{CStr, CString, c_char, c_int};
use std::fmt;
use std::io::BufRead;

/// 对应 C++ 的 `struct Query`
/// Rust 不知道内部结构，只拿指针用
#[repr(C)]
struct CQuery {
    _private: [u8; 0],
}

/// 算法参数，与 C++ `SparcConfig`（`sparc-source-code/sparc.h`）逐字段对应。
/// 增删/移动字段必须两端同步修改。
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SparcConfig {
    /// debug=true 时 C++ 层会向进程当前目录写调试文件
    /// （align_profile.txt / subgraph.dot / DEBUG.consensus.fasta 等）。
    pub debug: bool,
    /// k-mer 大小，范围 [1, 16]（k-mer 以 uint32 位编码，2 bit/碱基）。建议 [1, 2]。
    pub kmer: c_int,
    /// 覆盖度阈值（CLI 的 c），建议 [1, 5]。
    pub coverage_threshold: c_int,
    /// 打分方法：1 = 对数比例法，2 = 默认的线性减法。
    pub scoring_method: c_int,
    /// debug 模式输出子图 dot 文件的区间起点（其余情况下不生效）。
    pub subgraph_begin: c_int,
    /// debug 模式输出子图 dot 文件的区间终点（其余情况下不生效）。
    pub subgraph_end: c_int,
    /// 保留字段，当前实现未使用。
    pub cns_start: c_int,
    /// 保留字段，当前实现未使用。
    pub cns_end: c_int,
    /// 保留字段：仅对原 CLI 读入的 m5 行生效，不影响 FFI 传入的 query。
    pub report_begin: c_int,
    /// 保留字段：仅对原 CLI 读入的 m5 行生效，不影响 FFI 传入的 query。
    pub report_end: c_int,
    /// 覆盖度滑动窗口半径（原 CLI 固定 200，绑定默认 2）。
    pub cov_radius: c_int,
    /// 自适应阈值（CLI 的 t）。<0 关闭自适应（CLI 默认 -0.1），建议 [0.0, 0.3]。
    pub threshold: f64,
}

impl Default for SparcConfig {
    fn default() -> Self {
        Self {
            debug: false,
            kmer: 1,
            coverage_threshold: 2,
            scoring_method: 2,
            subgraph_begin: 0,
            subgraph_end: 0,
            cns_start: 0,
            cns_end: 0,
            report_begin: 0,
            report_end: 0,
            cov_radius: 2,
            threshold: 0.2,
        }
    }
}

/// consensus 结果。
///
/// `start`/`end` 是最优路径在 backbone 上的节点下标区间 `[start, end)`：
/// - `start == None`：路径头部不在 backbone 上（例如起点是一个插入分支节点）；
/// - `end == 0`：无可信路径（如无 query 输入），此时 `seq` 为原 backbone。
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Consensus {
    pub seq: String,
    pub start: Option<u32>,
    pub end: u32,
}

/// [`sparc_consensus`] 与 m5 解析的错误类型。
#[derive(Clone, Debug, PartialEq)]
pub enum SparcError {
    /// kmer 必须 >= 1，且 <= 16（k-mer 以 uint32 位编码）。
    InvalidKmer { kmer: i32 },
    /// backbone 长度不足以切出 k-mer。
    BackboneTooShort { backbone_len: usize, kmer: i32 },
    /// backbone 含非法碱基（仅允许 ACGTacgt；含 N 时上游会静默出错，这里显式拒绝）。
    InvalidBackboneBase { base: char, position: usize },
    /// cov_radius 为负会使 C++ 滑动窗口下标越界。
    InvalidCovRadius { cov_radius: i32 },
    /// threshold 为 NaN。
    InvalidThreshold { threshold: f64 },
    /// 比对串为空。
    EmptyAlignment,
    /// 两条比对串必须等长。
    AlignedLengthMismatch { query_len: usize, target_len: usize },
    /// 比对串含非法字符（仅允许 ACGTacgt 与 gap '-'）。
    InvalidAlignedBase { base: char, position: usize },
    /// target 坐标超出 backbone 范围或 start > end。
    InvalidTargetSpan {
        target_start: usize,
        target_end: usize,
        backbone_len: usize,
    },
    /// target 比对串中非 gap 碱基数与声明的 span 不一致。
    SpanMismatch { span: usize, target_bases: usize },
    /// m5 行解析失败。
    M5Parse(String),
}

impl fmt::Display for SparcError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SparcError::InvalidKmer { kmer } => {
                write!(f, "invalid kmer {kmer}: expected 1..=16")
            }
            SparcError::BackboneTooShort { backbone_len, kmer } => {
                write!(f, "backbone length {backbone_len} is shorter than kmer {kmer}")
            }
            SparcError::InvalidBackboneBase { base, position } => {
                write!(f, "invalid base {base:?} in backbone at position {position}: expected ACGT")
            }
            SparcError::InvalidCovRadius { cov_radius } => {
                write!(f, "invalid cov_radius {cov_radius}: expected >= 0")
            }
            SparcError::InvalidThreshold { threshold } => {
                write!(f, "invalid threshold {threshold}: must not be NaN")
            }
            SparcError::EmptyAlignment => write!(f, "aligned sequences must not be empty"),
            SparcError::AlignedLengthMismatch { query_len, target_len } => {
                write!(
                    f,
                    "aligned sequence length mismatch: query {query_len} != target {target_len}"
                )
            }
            SparcError::InvalidAlignedBase { base, position } => {
                write!(
                    f,
                    "invalid base {base:?} in aligned sequence at position {position}: expected ACGT or '-'"
                )
            }
            SparcError::InvalidTargetSpan {
                target_start,
                target_end,
                backbone_len,
            } => {
                write!(
                    f,
                    "invalid target span [{target_start}, {target_end}): backbone length is {backbone_len}"
                )
            }
            SparcError::SpanMismatch { span, target_bases } => {
                write!(
                    f,
                    "target aligned sequence has {target_bases} non-gap bases, but span is {span}"
                )
            }
            SparcError::M5Parse(msg) => write!(f, "m5 parse error: {msg}"),
        }
    }
}

impl std::error::Error for SparcError {}

/// ConsensusNode.kmer 为 uint32_t，2 bit/碱基，最多编码 16 个碱基。
const MAX_KMER: i32 = 16;

fn is_base(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
}

fn validate_inputs(
    backbone: &str,
    queries: &[Query],
    config: &SparcConfig,
) -> Result<(), SparcError> {
    if config.kmer < 1 || config.kmer > MAX_KMER {
        return Err(SparcError::InvalidKmer { kmer: config.kmer });
    }
    if config.cov_radius < 0 {
        return Err(SparcError::InvalidCovRadius {
            cov_radius: config.cov_radius,
        });
    }
    if config.threshold.is_nan() {
        return Err(SparcError::InvalidThreshold {
            threshold: config.threshold,
        });
    }
    if backbone.len() < config.kmer as usize {
        return Err(SparcError::BackboneTooShort {
            backbone_len: backbone.len(),
            kmer: config.kmer,
        });
    }
    for (position, &b) in backbone.as_bytes().iter().enumerate() {
        if !is_base(b) {
            return Err(SparcError::InvalidBackboneBase {
                base: b as char,
                position,
            });
        }
    }
    for query in queries {
        validate_query(query, backbone.len())?;
    }
    Ok(())
}

fn validate_query(query: &Query, backbone_len: usize) -> Result<(), SparcError> {
    let q = query.query_aligned_seq.as_bytes();
    let t = query.target_aligned_seq.as_bytes();
    if q.is_empty() || t.is_empty() {
        return Err(SparcError::EmptyAlignment);
    }
    if q.len() != t.len() {
        return Err(SparcError::AlignedLengthMismatch {
            query_len: q.len(),
            target_len: t.len(),
        });
    }
    // 字母表校验同时排除内嵌 NUL，保证后续 CString::new 不会失败
    for (position, &b) in q.iter().enumerate() {
        if !is_base(b) && b != b'-' {
            return Err(SparcError::InvalidAlignedBase {
                base: b as char,
                position,
            });
        }
    }
    let mut target_bases = 0usize;
    for (position, &b) in t.iter().enumerate() {
        if !is_base(b) && b != b'-' {
            return Err(SparcError::InvalidAlignedBase {
                base: b as char,
                position,
            });
        }
        if b != b'-' {
            target_bases += 1;
        }
    }
    if query.target_start > query.target_end || query.target_end > backbone_len {
        return Err(SparcError::InvalidTargetSpan {
            target_start: query.target_start,
            target_end: query.target_end,
            backbone_len,
        });
    }
    if target_bases != query.target_end - query.target_start {
        return Err(SparcError::SpanMismatch {
            span: query.target_end - query.target_start,
            target_bases,
        });
    }
    Ok(())
}

#[repr(C)]
struct CConsensusResult {
    seq: *mut u8,
    start_pos: c_int,
    end_pos: c_int,
}

impl Drop for CConsensusResult {
    fn drop(&mut self) {
        unsafe {
            SparcFreeConsensusResult(self.seq);
        }
    }
}

impl CConsensusResult {
    fn into_consensus(self) -> Consensus {
        let seq = unsafe { CStr::from_ptr(self.seq as *const c_char) }
            .to_string_lossy()
            .into_owned();
        Consensus {
            seq,
            start: if self.start_pos < 0 {
                None
            } else {
                Some(self.start_pos.max(0) as u32)
            },
            end: self.end_pos.max(0) as u32,
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
    ) -> CConsensusResult;

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

/// 一条 read 相对 backbone 的比对，对应一行 blasr m5 记录。
///
/// `query_aligned_seq` / `target_aligned_seq` 是带 `-` gap 列的比对串，
/// 两者必须等长；`target_aligned_seq` 中非 `-` 的碱基数必须等于
/// `target_end - target_start`。以上约束由 [`sparc_consensus`] 在调用时校验。
#[derive(Clone, Debug, PartialEq, Eq)]
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
    /// 由 target/query 两条比对串构造（参数顺序：先 target 后 query）。
    /// `target_aligned` 为 backbone 方向的比对串，坐标为 backbone 上的
    /// `[target_start, target_end)`。默认按正链（tStrand = '+'）处理。
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

    /// 标记该比对来自负链（m5 的 tStrand == '-'）。
    ///
    /// 传入的两条比对串保持 m5 负链行的原始方向（read 方向），
    /// C++ 层会对其做反向互补。
    pub fn reverse_strand(mut self) -> Self {
        self.rev_strand = true;
        self
    }

    /// 解析一行 blasr m5 记录（19 个空白分隔字段，见 README）。
    ///
    /// 只取绑定层用到的字段：qStart/qEnd、tStart/tEnd、tStrand、
    /// qAlignedSeq、tAlignedSeq（matchPattern 等忽略）。
    /// 字母表/坐标一致性由 [`sparc_consensus`] 校验。
    pub fn from_m5_row(row: &str) -> Result<Self, SparcError> {
        let fields: Vec<&str> = row.split_whitespace().collect();
        if fields.len() != 19 {
            return Err(SparcError::M5Parse(format!(
                "expected 19 whitespace-separated fields, got {}",
                fields.len()
            )));
        }
        // qName qLength qStart qEnd qStrand tName tLength tStart tEnd tStrand
        // score numMatch numMismatch numIns numDel mapQV
        // qAlignedSeq matchPattern tAlignedSeq
        let parse_usize = |s: &str, name: &str| {
            s.parse::<usize>()
                .map_err(|_| SparcError::M5Parse(format!("invalid {name}: {s:?}")))
        };
        let query_start = parse_usize(fields[2], "qStart")?;
        let query_end = parse_usize(fields[3], "qEnd")?;
        let target_start = parse_usize(fields[7], "tStart")?;
        let target_end = parse_usize(fields[8], "tEnd")?;
        let rev_strand = match fields[9] {
            "+" => false,
            "-" => true,
            other => return Err(SparcError::M5Parse(format!("invalid tStrand: {other:?}"))),
        };
        Ok(Self {
            query_aligned_seq: fields[16].to_string(),
            target_aligned_seq: fields[18].to_string(),
            rev_strand,
            query_start,
            query_end,
            target_start,
            target_end,
        })
    }

    fn fill_c_query(&self, c_query: *mut CQuery) {
        unsafe {
            let query_aligned_seq = CString::new(self.query_aligned_seq.as_bytes()).unwrap();
            QuerySetQueryAlignedSeq(c_query, query_aligned_seq.as_c_str().as_ptr());

            let target_aligned_seq = CString::new(self.target_aligned_seq.as_bytes()).unwrap();
            QuerySetTargetAlignedSeq(c_query, target_aligned_seq.as_c_str().as_ptr());

            QuerySetQueryStrand(c_query, b'+' as c_char);
            // C++ 层按 tStrand == '-' 对两条比对串做反向互补
            QuerySetTargetStrand(
                c_query,
                if self.rev_strand { b'-' } else { b'+' } as c_char,
            );

            QuerySetTargetLength(c_query, (self.target_end - self.target_start) as c_int);
            QuerySetTargetStart(c_query, self.target_start as c_int);
            QuerySetTargetEnd(c_query, self.target_end as c_int);
        }
    }
}

/// 逐行解析 m5 比对文件，跳过空行。
pub fn parse_m5<R: BufRead>(mut reader: R) -> Result<Vec<Query>, SparcError> {
    let mut queries = Vec::new();
    let mut line = String::new();
    loop {
        line.clear();
        let n = reader
            .read_line(&mut line)
            .map_err(|e| SparcError::M5Parse(e.to_string()))?;
        if n == 0 {
            return Ok(queries);
        }
        let row = line.trim_end();
        if row.is_empty() {
            continue;
        }
        queries.push(Query::from_m5_row(row)?);
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

/// 由 backbone 和一组比对计算 consensus。
///
/// 返回 [`Consensus`]，其 `start`/`end` 为最优路径在 backbone 上的
/// `[start, end)` 节点下标区间（`end == 0` 表示回退为原 backbone，
/// `start == None` 表示路径头部不在 backbone 上）。
///
/// 所有输入在进入 C++ 层前完成校验；不满足约束时返回 [`SparcError`]。
pub fn sparc_consensus(
    backbone: &str,
    queries: &[Query],
    config: &SparcConfig,
) -> Result<Consensus, SparcError> {
    validate_inputs(backbone, queries, config)?;

    let sparc_queries = queries
        .iter()
        .map(SparcQuery::from)
        .collect::<Vec<SparcQuery>>();
    let mut c_queries = sparc_queries
        .iter()
        .map(|v| v.c_query)
        .collect::<Vec<*mut CQuery>>();

    // 字母表校验已排除内嵌 NUL，CString 不会失败
    let backbone_str = CString::new(backbone).expect("NUL is rejected by validation");
    let result = unsafe {
        SparcConsensus(
            backbone_str.as_ptr(),
            c_queries.as_mut_ptr(),
            c_queries.len() as c_int,
            config as *const SparcConfig,
        )
    };

    Ok(result.into_consensus())
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

    fn revcomp(s: &str) -> String {
        s.bytes()
            .rev()
            .map(|b| match b {
                b'A' => 'T',
                b'T' => 'A',
                b'C' => 'G',
                b'G' => 'C',
                _ => '-',
            })
            .collect()
    }

    #[test]
    fn test_sparc_consensus() {
        let backbone = "GATCGGGCTAA";
        let config = make_config(backbone.len());
        let queries = make_queries(backbone);

        let consensus = sparc_consensus(backbone, &queries, &config).unwrap();
        assert!(!consensus.seq.is_empty());
        assert_eq!(consensus.start, Some(0));
        assert!(consensus.end > 0);
    }

    /// 回归：max_score == 0 的早退路径必须先释放 k-mer 图和 ref.read_bits。
    /// 此前该路径在 SparcFreeInfo/free 之前 return，空 queries 即可触发整体泄漏。
    #[test]
    fn test_empty_queries_returns_backbone() {
        let backbone = "GATCGGGCTAA";
        let config = make_config(backbone.len());

        let consensus = sparc_consensus(backbone, &[], &config).unwrap();
        assert_eq!(consensus.seq, backbone);
        assert_eq!(consensus.start, None);
        assert_eq!(consensus.end, 0);
    }

    /// 回归：SparcFreeInfo 须覆盖最后一个 backbone 节点（每次调用都会走到该路径），
    /// 同时验证跨调用的 FFI 配对（NewQuery/FreeQuery、malloc/free）不累积泄漏。
    #[test]
    fn test_sparc_consensus_repeated() {
        let backbone = "GATCGGGCTAA";
        let config = make_config(backbone.len());
        let queries = make_queries(backbone);

        for _ in 0..50 {
            let consensus = sparc_consensus(backbone, &queries, &config).unwrap();
            assert!(!consensus.seq.is_empty());
            assert!(consensus.start.expect("start should be set") < consensus.end);
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

        let consensus = sparc_consensus(backbone, &queries, &config).unwrap();
        assert!(!consensus.seq.is_empty());
        assert!(consensus.start.expect("start should be set") < consensus.end);
    }

    /// 回归：比对起始处即 mismatch 的 read 会在链首产生不在任何 backbone
    /// 右子图内的孤儿节点；SparcFreeInfo 需经由 orphan_nodes 释放
    /// （内存层面由 CI 的 valgrind job 兜底验证）。
    #[test]
    fn test_leading_mismatch_queries() {
        let backbone = "GATCGGGCTAA";
        let config = make_config(backbone.len());
        let queries: Vec<Query> = (0..1000)
            .map(|_| Query {
                query_aligned_seq: "TATCGGGCTAA".to_string(),
                target_aligned_seq: backbone.to_string(),
                rev_strand: false,
                query_start: 0,
                query_end: 11,
                target_start: 0,
                target_end: backbone.len(),
            })
            .collect();

        let consensus = sparc_consensus(backbone, &queries, &config).unwrap();
        assert!(!consensus.seq.is_empty());
    }

    /// 负链：m5 负链行的两条比对串按 read 方向给出，C++ 层反向互补后
    /// 应与直接给正向比对得到完全一致的结果。
    #[test]
    fn test_reverse_strand_matches_forward() {
        let backbone = "GATCGGGCTAA";
        let config = make_config(backbone.len());

        let forward = make_queries(backbone);
        let reversed: Vec<Query> = forward
            .iter()
            .map(|q| Query {
                query_aligned_seq: revcomp(&q.query_aligned_seq),
                target_aligned_seq: revcomp(&q.target_aligned_seq),
                rev_strand: true,
                query_start: q.query_start,
                query_end: q.query_end,
                target_start: q.target_start,
                target_end: q.target_end,
            })
            .collect();

        let a = sparc_consensus(backbone, &forward, &config).unwrap();
        let b = sparc_consensus(backbone, &reversed, &config).unwrap();
        assert_eq!(a, b);
    }

    /// FFI 布局回归：Rust `SparcConfig` 必须与 C++ `struct SparcConfig`
    /// 逐字段同偏移（已用 C++ offsetof 打点核对），错位会让参数静默串位。
    /// C++ 侧由 sparc.h 中的 static_assert 同步守护。
    #[test]
    fn test_sparc_config_ffi_layout() {
        assert_eq!(std::mem::size_of::<SparcConfig>(), 56);
        assert_eq!(std::mem::offset_of!(SparcConfig, debug), 0);
        assert_eq!(std::mem::offset_of!(SparcConfig, kmer), 4);
        assert_eq!(std::mem::offset_of!(SparcConfig, coverage_threshold), 8);
        assert_eq!(std::mem::offset_of!(SparcConfig, scoring_method), 12);
        assert_eq!(std::mem::offset_of!(SparcConfig, subgraph_begin), 16);
        assert_eq!(std::mem::offset_of!(SparcConfig, subgraph_end), 20);
        assert_eq!(std::mem::offset_of!(SparcConfig, cns_start), 24);
        assert_eq!(std::mem::offset_of!(SparcConfig, cns_end), 28);
        assert_eq!(std::mem::offset_of!(SparcConfig, report_begin), 32);
        assert_eq!(std::mem::offset_of!(SparcConfig, report_end), 36);
        assert_eq!(std::mem::offset_of!(SparcConfig, cov_radius), 40);
        assert_eq!(std::mem::offset_of!(SparcConfig, threshold), 48);
    }

    /// threshold 行为回归：自适应阈值必须真正从 FFI 传到 C++。
    /// threshold=100 使所有边的 new_score 为负、从不更新，全部节点 score=0，
    /// 触发 fallback（原样返回 backbone）；threshold=-0.1 则正常产出 consensus。
    #[test]
    fn test_threshold_reaches_cpp() {
        let backbone = "GATCGGGCTAA";
        let queries = make_queries(backbone);

        let mut config = make_config(backbone.len());
        config.kmer = 1;
        config.threshold = -0.1;
        let consensus = sparc_consensus(backbone, &queries, &config).unwrap();
        assert_ne!(consensus.seq, backbone);
        assert!(consensus.end > 0);

        let mut config = make_config(backbone.len());
        config.kmer = 1;
        config.threshold = 100.0;
        let consensus = sparc_consensus(backbone, &queries, &config).unwrap();
        assert_eq!(consensus.seq, backbone);
        assert_eq!(consensus.end, 0);
        assert_eq!(consensus.start, None);
    }

    /// 回归：k=3 且 cov_radius 大于 node_vec 大小时，C++ 滑动窗口初始化循环
    /// 会越界读 node_vec（上游既有 bug，k=1/2 恰好不触发）。修复后该组合
    /// 必须良定义且 valgrind 干净（CI 兜底），输出确定性。
    #[test]
    fn test_k3_large_radius() {
        let bb = "GATCGGGCTAAACGTACGATCGATCGATCGGATCCGATACGTACGTACGATCGTACGATC";
        let mut config = make_config(bb.len());
        config.kmer = 3;
        config.cov_radius = 500;
        let exact = Query {
            query_aligned_seq: bb.to_string(),
            target_aligned_seq: bb.to_string(),
            rev_strand: false,
            query_start: 0,
            query_end: bb.len(),
            target_start: 0,
            target_end: bb.len(),
        };
        let mut queries = vec![exact.clone(), exact.clone(), exact];
        // 首列替换：产生孤儿链首；末端插入：产生末端分支
        queries.push(Query {
            query_aligned_seq: format!("T{}", &bb[1..]),
            target_aligned_seq: bb.to_string(),
            rev_strand: false,
            query_start: 0,
            query_end: bb.len(),
            target_start: 0,
            target_end: bb.len(),
        });
        queries.push(Query {
            query_aligned_seq: format!("{bb}A"),
            target_aligned_seq: format!("{bb}-"),
            rev_strand: false,
            query_start: 0,
            query_end: bb.len() + 1,
            target_start: 0,
            target_end: bb.len(),
        });

        let first = sparc_consensus(bb, &queries, &config).unwrap();
        assert!(!first.seq.is_empty());
        // 确定性：同输入重复调用结果一致
        for _ in 0..3 {
            let again = sparc_consensus(bb, &queries, &config).unwrap();
            assert_eq!(again, first);
        }
    }

    #[test]
    fn test_validation_errors() {
        let backbone = "GATCGGGCTAA";
        let config = make_config(backbone.len());
        let q = |qs: &str, ts: &str, s: usize, e: usize| Query {
            query_aligned_seq: qs.to_string(),
            target_aligned_seq: ts.to_string(),
            rev_strand: false,
            query_start: 0,
            query_end: qs.len(),
            target_start: s,
            target_end: e,
        };

        // kmer 越界（下界 / 上界）
        for kmer in [0, 17] {
            let mut cfg = make_config(backbone.len());
            cfg.kmer = kmer;
            assert_eq!(
                sparc_consensus(backbone, &[], &cfg),
                Err(SparcError::InvalidKmer { kmer })
            );
        }

        // backbone 长度不足
        let mut cfg = make_config(3);
        cfg.kmer = 5;
        assert_eq!(
            sparc_consensus("GAT", &[], &cfg),
            Err(SparcError::BackboneTooShort {
                backbone_len: 3,
                kmer: 5
            })
        );

        // backbone 含 N
        assert_eq!(
            sparc_consensus("GATCGNGCTAA", &[], &config),
            Err(SparcError::InvalidBackboneBase {
                base: 'N',
                position: 5
            })
        );

        // cov_radius 为负
        let mut cfg = make_config(backbone.len());
        cfg.cov_radius = -1;
        assert_eq!(
            sparc_consensus(backbone, &[], &cfg),
            Err(SparcError::InvalidCovRadius { cov_radius: -1 })
        );

        // threshold 为 NaN（NaN != NaN，用 matches! 而非 assert_eq! 比较）
        let mut cfg = make_config(backbone.len());
        cfg.threshold = f64::NAN;
        assert!(matches!(
            sparc_consensus(backbone, &[], &cfg),
            Err(SparcError::InvalidThreshold { .. })
        ));

        // 空比对
        assert_eq!(
            sparc_consensus(backbone, &[q("", "", 0, 11)], &config),
            Err(SparcError::EmptyAlignment)
        );

        // 长度不齐
        assert_eq!(
            sparc_consensus(
                backbone,
                &[q("GATC", "GATCG", 0, 5)],
                &config
            ),
            Err(SparcError::AlignedLengthMismatch {
                query_len: 4,
                target_len: 5
            })
        );

        // 非法碱基
        assert_eq!(
            sparc_consensus(backbone, &[q("GATCN", "GATCN", 0, 5)], &config),
            Err(SparcError::InvalidAlignedBase {
                base: 'N',
                position: 4
            })
        );

        // span 越界
        assert_eq!(
            sparc_consensus(
                backbone,
                &[q("GATCG", "GATCG", 0, 99)],
                &config
            ),
            Err(SparcError::InvalidTargetSpan {
                target_start: 0,
                target_end: 99,
                backbone_len: backbone.len(),
            })
        );

        // span 与非 gap 碱基数不一致（含一个 gap 列）
        assert_eq!(
            sparc_consensus(
                backbone,
                &[q("GATC-GGCTAA", "GATC-GGCTAA", 0, 11)],
                &config
            ),
            Err(SparcError::SpanMismatch {
                span: 11,
                target_bases: 10
            })
        );
    }

    #[test]
    fn test_parse_m5_row() {
        let row = "read_24/0_10 10 0 12 + ref_template 11 0 11 + 100 10 0 0 0 60 \
                   GATCGCGCTAA 10M GATCGGGCTAA";
        let q = Query::from_m5_row(row).unwrap();
        assert_eq!(q.query_aligned_seq, "GATCGCGCTAA");
        assert_eq!(q.target_aligned_seq, "GATCGGGCTAA");
        assert_eq!(q.query_start, 0);
        assert_eq!(q.query_end, 12);
        assert_eq!(q.target_start, 0);
        assert_eq!(q.target_end, 11);
        assert!(!q.rev_strand);

        let neg_row = "read_24/0_10 10 0 12 + ref_template 11 0 11 - 100 10 0 0 0 60 \
                       GATCGCGCTAA 10M GATCGGGCTAA";
        let q = Query::from_m5_row(neg_row).unwrap();
        assert!(q.rev_strand);

        let bad = Query::from_m5_row("only three fields");
        assert!(matches!(bad, Err(SparcError::M5Parse(_))));

        let bad_strand = Query::from_m5_row(
            "r 10 0 12 + t 11 0 11 ? 100 10 0 0 0 60 GATCGCGCTAA 10M GATCGGGCTAA",
        );
        assert!(matches!(bad_strand, Err(SparcError::M5Parse(_))));
    }

    #[test]
    fn test_parse_m5_testdata() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/testdata2/backbone-0.mapped.m5");
        let queries = parse_m5(std::io::BufReader::new(std::fs::File::open(path).unwrap())).unwrap();
        assert_eq!(queries.len(), 7);
        assert!(queries.iter().all(|q| q.target_start == 0 && !q.rev_strand));
    }

    /// 端到端：testdata2 的 backbone + m5 应复现原版 CLI 的 out.consensus.fasta。
    #[test]
    fn test_e2e_testdata2() {
        let dir = env!("CARGO_MANIFEST_DIR");
        let fasta =
            std::fs::read_to_string(format!("{dir}/testdata2/backbone-0.fasta")).unwrap();
        let backbone: String = fasta.lines().skip(1).collect();
        assert_eq!(backbone, "GATCGGGCTAA");

        let f = std::fs::File::open(format!("{dir}/testdata2/backbone-0.mapped.m5")).unwrap();
        let queries = parse_m5(std::io::BufReader::new(f)).unwrap();

        // 对齐原 CLI 生成 out.consensus.fasta 时的参数：
        // k=1（默认），threshold=-0.1 关闭自适应，滑动窗口半径 200
        let mut config = make_config(backbone.len());
        config.kmer = 1;
        config.threshold = -0.1;
        config.cov_radius = 200;

        let consensus = sparc_consensus(&backbone, &queries, &config).unwrap();
        let expected =
            std::fs::read_to_string(format!("{dir}/testdata2/out.consensus.fasta")).unwrap();
        let expected_seq: String = expected.lines().skip(1).collect();
        assert_eq!(consensus.seq, expected_seq);
    }
}
