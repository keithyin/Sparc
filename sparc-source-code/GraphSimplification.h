#pragma once
#include <iostream>
#include <string>
#include <string.h>
#include <stdint.h>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include <fstream>
#include "time.h"
#include "BasicDataStructure.h"
#include "GraphConstruction.h"
using namespace std;

char SparcKmer2base(uint32_t kmer);

void SparcMergeNodes(struct Backbone *backbone_info);

void SparcOutputPathsFromANode(ConsensusNode *begin_node, string filename, map<ConsensusNode *, bool> &Visited);

void SparcOutputSubGraph(struct Backbone *backbone_info, int begin, int end, string filename);

char *SparcMultiply(const char *a_in, const char *b_in, char *mul);

void SparcBFSFindBestPath(struct Backbone *backbone_info, int node_idx);

void SparcFindBestPath(struct Backbone *backbone_info);

void SparcBFSFree(struct Backbone *backbone_info, int node_idx);

void SparcBFSClear(struct Backbone *backbone_info, int node_idx);

void SparcClearInfo(struct Backbone *backbone_info);

void SparcFreeInfo(struct Backbone *backbone_info);
