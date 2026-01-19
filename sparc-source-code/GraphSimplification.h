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

char kmer2base(uint32_t kmer);

void MergeNodes(struct Backbone *backbone_info);

void OutputPathsFromANode(ConsensusNode *begin_node, string filename, map<ConsensusNode *, bool> &Visited);

void OutputSubGraph(struct Backbone *backbone_info, int begin, int end, string filename);

char *multiply(const char *a_in, const char *b_in, char *mul);

void BFSFindBestPath(struct Backbone *backbone_info, int node_idx);

void FindBestPath(struct Backbone *backbone_info);

void BFSFree(struct Backbone *backbone_info, int node_idx);

void BFSClear(struct Backbone *backbone_info, int node_idx);

void ClearInfo(struct Backbone *backbone_info);

void FreeInfo(struct Backbone *backbone_info);
