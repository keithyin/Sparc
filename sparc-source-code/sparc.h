

#include <string>
#include "BasicDataStructure.h"

struct SparcConfig
{
    bool debug;
    int kmer;
    int converage_threshold;
    int scoring_method;
    int subgraph_begin;
    int subgraph_end;
    int cns_start;
    int cns_end;
    int report_begin;
    int report_end;
    int cov_radius;
};

struct SparcConsensusResult {
    char* seq;
};

// std::string SparcConsensus();
extern "C"
{
    SparcConsensusResult SparcConsensus(char *backbone_c, Query **queries, int n_queries, SparcConfig *config);

    void SparcFreeConsensusResult(char *consensus);
}
