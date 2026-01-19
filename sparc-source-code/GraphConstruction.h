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
#include <sstream>
#include "time.h"
#include "Align.h"
#include "BasicDataStructure.h"
using namespace std;

void Consensus_Kmer_Graph_Construction(struct RefRead *read, struct Backbone *backbone_info, int K_size);

// 将 mismatch 拆成了两个 GAP
// gap 右对齐
void NormalizeAlignment(Query *query_info);

void PatchGaps(Query *query_info);

void FillGaps(Query *query_info);

void Add_Path_To_Backbone(struct Backbone *backbone_info, struct Query *query_info, int K_size);