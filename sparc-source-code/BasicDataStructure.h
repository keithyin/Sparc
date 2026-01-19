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
// typedef unsigned __int64 uint64_t;
// typedef __int64 int64_t;
// typedef unsigned __int8 uint8_t;
// typedef unsigned __int16 uint16_t;
// typedef unsigned __int32 uint32_t;
// typedef __int32 int32_t;

// These are the structures to save the k-mers.32 bases, 64 bases, 96 bases, 128 bases.
struct kmer_t
{
	uint64_t kmer;
};

struct kmer_t2
{
	uint64_t kmer[2];
};

struct kmer_t3
{
	uint64_t kmer[3];
};

struct kmer_t4
{
	uint64_t kmer[4];
};

// read structure

struct Read
{
	char tag[1000];
	bool error_nt[1000];
	char c_seq[10000]; // char representation
	// char *c_seq;//char representation

	// uint64_t read_bits[10000];//bit representation
	uint64_t *read_bits;
	// char read[1000];//char representation
	int readLen; // read length
	int read_idx;
};

struct RefRead
{
	char tag[1000];
	uint64_t *read_bits; // bit representation
	size_t read_idx;
	int alloc_sz;
	int readLen; // read length
	int contig_no;
};

struct AlignProfile
{
	vector<int> match_vec, mismatch_vec, deletion_vec, insertion_vec;
};

struct Query
{
	string qName, tName, qAlignedSeq, matchPattern, tAlignedSeq;
	int qLength, qStart, qEnd, tLength, tStart, tEnd, score, numMatch, numMismatch, numIns, numDel, mapQV;
	char qStrand, tStrand;
	size_t read_idx;
	int report_b, report_e;
	int n_exist;
	int n_new;
	bool Patch, Fill;
	int Patch_K;
	int Patch_D;
	int Patch_G;
};

extern "C"
{
	Query *NewQuery();
	void FreeQuery(Query *query);

	/* ---------- string fields ---------- */
	void QuerySetQueryName(Query &query, const char *name);
	void QuerySetTargetName(Query &query, const char *name);
	void QuerySetQueryAlignedSeq(Query &query, const char *seq);
	void QuerySetMatchPattern(Query &query, const char *pattern);
	void QuerySetTargetAlignedSeq(Query &query, const char *seq);

	/* ---------- basic attributes ---------- */
	void QuerySetQueryLength(Query &query, int queryLength);
	void QuerySetQueryStart(Query &query, int qStart);
	void QuerySetQueryEnd(Query &query, int qEnd);

	void QuerySetTargetLength(Query &query, int targetLength);
	void QuerySetTargetStart(Query &query, int tStart);
	void QuerySetTargetEnd(Query &query, int tEnd);

	/* ---------- alignment stats ---------- */
	void QuerySetScore(Query &query, int score);
	void QuerySetNumMatch(Query &query, int numMatch);
	void QuerySetNumMismatch(Query &query, int numMismatch);
	void QuerySetNumIns(Query &query, int numIns);
	void QuerySetNumDel(Query &query, int numDel);
	void QuerySetMapQV(Query &query, int mapQV);

	/* ---------- strand / index ---------- */
	void QuerySetQueryStrand(Query &query, char strand);
	void QuerySetTargetStrand(Query &query, char strand);
	void QuerySetReadIndex(Query &query, size_t read_idx);

	/* ---------- report range ---------- */
	void QuerySetReportBegin(Query &query, int report_b);
	void QuerySetReportEnd(Query &query, int report_e);

	/* ---------- counters ---------- */
	void QuerySetNumExist(Query &query, int n_exist);
	void QuerySetNumNew(Query &query, int n_new);

	/* ---------- patch flags ---------- */
	void QuerySetPatch(Query &query, bool patch);
	void QuerySetFill(Query &query, bool fill);

	/* ---------- patch parameters ---------- */
	void QuerySetPatchK(Query &query, int k);
	void QuerySetPatchD(Query &query, int d);
	void QuerySetPatchG(Query &query, int g);
}

struct ConsensusEdgeNode
{
	int edge_cov;
	struct ConsensusNode *node_ptr;
	struct ConsensusEdgeNode *nxt_edge;
};

struct sparse_consensus_edge_node
{
	uint32_t edge;
	int32_t edge_cov : 24, len : 8;
	struct sparse_consensus_node *node_ptr;
	struct sparse_consensus_edge_node *nxt_edge;
};

struct ConsensusNode
{
	uint32_t kmer;
	uint32_t coord;
	int64_t score;
	uint32_t cns_coord;
	uint32_t cov : 28, used : 1, in_backbone : 1, in_cns_backbone : 1;
	bool selected;
	ConsensusEdgeNode *left;
	ConsensusEdgeNode *right;
	ConsensusNode *last_node;
};

struct sparse_consensus_node
{
	uint32_t kmer;
	int64_t score;
	uint32_t coord;
	uint32_t cns_coord;
	uint32_t cov : 28, used : 1, in_backbone : 1, in_cns_backbone : 1;
	sparse_consensus_edge_node *left;
	sparse_consensus_edge_node *right;
	sparse_consensus_node *last_node;
};

struct Backbone
{
	vector<ConsensusNode *> node_vec;
	vector<sparse_consensus_node *> sparse_node_vec;
	vector<uint64_t> cov_vec;
	vector<uint64_t> cnt_vec;
	string backbone;
	int64_t n_nodes;
	int64_t n_edges;
	int CovTh;
	int ScoringMethod;
	size_t ref_matched, ref_mismatch;
	int boost;
	int gap;
	double threshold;
	string ContigPrefix;
};

struct read_index
{
	vector<map<uint64_t, bool>> repeat_maps;
	uint64_t repeat_cnt;
	int MaxMatch;
	vector<int> read_len_vt;
};

struct reads_table
{
	bool in_use;
	list<uint64_t *> pblocks;
	int BytesPerBlock;
	int current_block;
	int current_byte;
	int current_read;
	map<int, uint64_t *> read_ptr;
	vector<int> read_len_vt;
};

// contig graph

struct contigs_info
{
	int total_contigs;
	int K_size;
	vector<int> contig_sz_vt, kmer_cnt_vt, comp_vt;
	vector<int> contigs_hp_b, contigs_hp_e;
	vector<string> contigs_str;
	map<int, vector<int>> Length_ID;
	map<int, vector<int>> Cov_Length;
	// map<int, vector<int> > ctg_in_scf;
	// vector<vector<int> > scaffolds;
	// vector<vector<int> > gaps_in_scaffolds;

	vector<vector<int>::iterator> LengthRank;
	vector<int> cov_vt;
	vector<struct c_info> c_info_vt;
	vector<map<int, struct scaffold_contig_info>> scaffold_adjacency_left, scaffold_adjacency_right;
	vector<map<int, struct adjacent_contig_info>> contig_adjacency_left, contig_adjacency_right;
};

// path information in the BFS bubble removal

bool get_a_fasta_read(ifstream &fasta_in, string &tag, string &str, string &n_tag);

bool get_a_fastq_read(ifstream &fastq_in, string &tag, string &seq, string &quality);

// left shift and right shift of shift_sz bits of the whole bit array, arr_sz is the array length
static inline void L_shift_NB(uint64_t *bitsarr, int shift_sz, int arr_sz)
{

	uint64_t temp_arr[100];

	/*
	for (int i=0;i<arr_sz;++i)
	{
		temp_arr[i]=0;
	}
	memset(temp_arr,0,sizeof(uint64_t)*arr_sz);
	*/

	int jmp = shift_sz / 64;
	int offset = shift_sz % 64;

	for (int i = 0; i < arr_sz; ++i)
	{
		if (i + jmp + 1 < arr_sz)
		{

			uint64_t tt = 0;
			if (offset == 0)
			{
				tt = 0;
			}
			else
			{
				tt = (bitsarr[i + jmp + 1] >> (64 - offset));
			}
			temp_arr[i] = ((bitsarr[i + jmp] << offset) | tt);
		}
		if (i + jmp + 1 == arr_sz)
		{
			temp_arr[i] = bitsarr[i + jmp] << offset;
		}
		if (i + jmp + 1 > arr_sz)
		{
			temp_arr[i] = 0;
		}
	}

	memcpy(bitsarr, temp_arr, sizeof(uint64_t) * arr_sz);
	/*
		for (int i=0;i<arr_sz;++i)
		{
			bitsarr[i]=temp_arr[i];
		}
		*/
}

static inline void R_shift_NB(uint64_t *bitsarr, int shift_sz, int arr_sz)
{
	uint64_t temp_arr[100];
	/*
	for (int i=0;i<arr_sz;++i)
	{
		temp_arr[i]=0;
	}
	memset(temp_arr,0,sizeof(uint64_t)*arr_sz);
	*/
	int jmp = shift_sz / 64;
	int offset = shift_sz % 64;

	for (int i = arr_sz - 1; i >= 0; --i)
	{
		if (i - jmp > 0)
		{
			if (offset > 0)
			{
				temp_arr[i] = (bitsarr[i - jmp] >> offset) | (bitsarr[i - jmp - 1] << (64 - offset));
			}
			else
			{
				temp_arr[i] = bitsarr[i - jmp];
			}
		}
		if (i - jmp == 0)
		{
			if (offset > 0)
			{
				temp_arr[i] = (bitsarr[i - jmp] >> offset);
			}
			else
			{
				temp_arr[i] = bitsarr[i - jmp];
			}
		}
		if (i - jmp < 0)
		{
			temp_arr[i] = 0;
		}
	}
	memcpy(bitsarr, temp_arr, sizeof(uint64_t) * arr_sz);
	/*
	for (int i=0;i<arr_sz;++i)
	{
		bitsarr[i]=temp_arr[i];
	}
	*/
}

// get reverse complement of a k-mer.
static inline uint64_t get_rev_comp_seq(uint64_t seq, int seq_size)
{
	seq = ~seq;

	seq = ((seq & 0x3333333333333333) << 2) | ((seq & 0xCCCCCCCCCCCCCCCC) >> 2);
	seq = ((seq & 0x0F0F0F0F0F0F0F0F) << 4) | ((seq & 0xF0F0F0F0F0F0F0F0) >> 4);
	seq = ((seq & 0x00FF00FF00FF00FF) << 8) | ((seq & 0xFF00FF00FF00FF00) >> 8);
	seq = ((seq & 0x0000FFFF0000FFFF) << 16) | ((seq & 0xFFFF0000FFFF0000) >> 16);
	seq = ((seq & 0x00000000FFFFFFFF) << 32) | ((seq & 0xFFFFFFFF00000000) >> 32);

	return seq >> (64 - (seq_size * 2));
}

static inline uint64_t *get_rev_comp_seq_arr(uint64_t *seq_arr, int seq_size, int arr_sz)
{

	int tot_bits = arr_sz * 64;
	for (int i = 0; i < arr_sz; ++i)
	{
		seq_arr[i] = ~seq_arr[i];
		seq_arr[i] = ((seq_arr[i] & 0x3333333333333333) << 2) | ((seq_arr[i] & 0xCCCCCCCCCCCCCCCC) >> 2);
		seq_arr[i] = ((seq_arr[i] & 0x0F0F0F0F0F0F0F0F) << 4) | ((seq_arr[i] & 0xF0F0F0F0F0F0F0F0) >> 4);
		seq_arr[i] = ((seq_arr[i] & 0x00FF00FF00FF00FF) << 8) | ((seq_arr[i] & 0xFF00FF00FF00FF00) >> 8);
		seq_arr[i] = ((seq_arr[i] & 0x0000FFFF0000FFFF) << 16) | ((seq_arr[i] & 0xFFFF0000FFFF0000) >> 16);
		seq_arr[i] = ((seq_arr[i] & 0x00000000FFFFFFFF) << 32) | ((seq_arr[i] & 0xFFFFFFFF00000000) >> 32);
	}

	int j = 0, k = arr_sz - 1;
	for (; j < k; ++j, --k)
	{
		uint64_t temp;
		temp = seq_arr[j];
		seq_arr[j] = seq_arr[k];
		seq_arr[k] = temp;
	}

	R_shift_NB(seq_arr, tot_bits - (seq_size * 2), arr_sz);
	return seq_arr;
	// return seq >> (64 - (seq_size*2));
}

// get sub bit array from a bit array.
inline void get_sub_arr(uint64_t *bitsarr_in, int bitsarr_len, int begin_pos, int sub_sz, uint64_t *bitsarr_out)
{
	if (bitsarr_len < sub_sz)
	{
		cout << "Error! Input kmer too short." << bitsarr_len << " " << sub_sz << endl;
		return;
	}
	int arr_sz_in = bitsarr_len / 32 + 1;
	int rem = bitsarr_len % 32;
	if (rem == 0)
	{
		arr_sz_in--;
	}

	int arr_sz_out = sub_sz / 32 + 1;
	if (sub_sz % 32 == 0)
	{
		arr_sz_out--;
	}

	uint64_t temp_arr[10];
	memset(temp_arr, 0, sizeof(temp_arr));

	memset(bitsarr_out, 0, sizeof(uint64_t) * arr_sz_out);

	int rem2 = (32 - rem + begin_pos) % 32;
	int block_beg = (32 - rem + begin_pos) / 32;
	if (rem == 0)
	{
		block_beg--;
	}

	int rem3 = (32 - rem + begin_pos + sub_sz) % 32;
	int block_end = (32 - rem + begin_pos + sub_sz) / 32;
	if (rem3 != 0)
	{
		rem3 = 32 - rem3;
	}
	else
	{
		block_end--;
	}
	if (rem == 0)
	{
		block_end--;
	}

	int orig_sz = (block_end - block_beg + 1);
	memcpy(temp_arr, &bitsarr_in[block_beg], orig_sz * sizeof(uint64_t));
	L_shift_NB(temp_arr, rem2 * 2, orig_sz);
	R_shift_NB(temp_arr, (rem2 + rem3) % 32 * 2, arr_sz_out);
	memcpy(bitsarr_out, temp_arr, arr_sz_out * sizeof(uint64_t));
}

uint64_t *str2bitsarr(const char *c_str, int len, uint64_t *b_str, int arr_sz);

// hash functions
uint64_t MurmurHash64A(const void *key, int len, unsigned int seed);

uint64_t MurmurHash64B(const void *key, int len, unsigned int seed);

// convert a string of nucleotide bases into bit array.
void Init_Read(string &seq, struct Read &read);

void Init_Ref_Read(string &seq, struct RefRead &read);

void reverse_complement_str(string &str);

// get the complement of a string of nucleotide bases
void complement_str(string &str);

// convert a bit array into a string
char *bitsarr2str(uint64_t *b_seq, int len, char *c_str, int arr_sz);
