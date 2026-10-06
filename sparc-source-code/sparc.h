
#include <cstddef>
#include <string>
#include "BasicDataStructure.h"

struct SparcConfig
{
    bool debug;
    int kmer;
    int coverage_threshold;
    int scoring_method;
    int subgraph_begin;
    int subgraph_end;
    int cns_start;
    int cns_end;
    int report_begin;
    int report_end;
    int cov_radius;
    // 自适应阈值（对应 CLI 的 t）。<0 关闭自适应，建议范围 [0.0, 0.3]
    double threshold;
};

// Rust 侧 (src/lib.rs) 的 #[repr(C)] SparcConfig 依赖此布局，
// 两端由 static_assert / offset_of! 测试同步守护，改动字段必须两端一致。
static_assert(sizeof(SparcConfig) == 56, "SparcConfig FFI layout changed");
static_assert(offsetof(SparcConfig, threshold) == 48, "SparcConfig FFI layout changed");
static_assert(offsetof(SparcConfig, cov_radius) == 40, "SparcConfig FFI layout changed");
static_assert(offsetof(SparcConfig, kmer) == 4, "SparcConfig FFI layout changed");

struct SparcConsensusResult {
    char* seq;
    int start_pos;
    int end_pos; // exclusive
};

// std::string SparcConsensus();
extern "C"
{
    SparcConsensusResult SparcConsensus(char *backbone_c, Query **queries, int n_queries, SparcConfig *config);

    void SparcFreeConsensusResult(char *consensus);
}
