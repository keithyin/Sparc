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

void SparcConsensusKmerGraphConstruction(struct RefRead *read, struct Backbone *backbone_info, int K_size);

// 将 mismatch 拆成了两个 GAP
// gap 右对齐
void SparcNormalizeAlignment(Query *query_info);

void SparcPatchGaps(Query *query_info);

void SparcFillGaps(Query *query_info);

void SparcAddPathToBackbone(struct Backbone *backbone_info, struct Query *query_info, int K_size);