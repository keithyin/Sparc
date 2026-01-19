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
using namespace std;

struct aln_t
{
	const static int MaxLen = 301;
	int score[MaxLen][MaxLen];
	uint8_t b[MaxLen][MaxLen];
	int insertions, deletions, substitutions;
};
void str2int(char *s, int *out);

int GlobalAlign(struct aln_t *aln_t, char *A, char *B, int match, int mismatch, int gap_cost, int band_width);

void printAlign(struct aln_t *aln_t, char *A, char *B, string &A_aln, string &B_aln);