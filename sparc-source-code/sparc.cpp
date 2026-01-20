#include "iostream"
#include "stdio.h"
#include "string"
#include "vector"
#include "cstdlib"
#include "bitset"
#include <map>
#include <math.h>
#include "memory"
#include <algorithm>
#include "fstream"
#include "sstream"
#include "list"
#include "stdlib.h"
#include "time.h"
#include "stdint.h"
#include "BasicDataStructure.h"
#include "GraphConstruction.h"
#include "GraphSimplification.h"
#include "sparc.h"

char *SparcConsensus(char *backbone_c, Query **queries, int n_queries, SparcConfig *config)
{
	std::string backbone = backbone_c;
	bool HELP = 0;
	int kmer = config->kmer;
	int CovTh = config->converage_threshold;
	int ScoringMethod = config->scoring_method;
	int subgraph_begin = config->subgraph_begin;
	int subgraph_end = config->subgraph_end;
	int cns_begin = config->cns_start;
	int cns_end = config->cns_end;
	int report_begin = config->report_begin;
	int report_end = config->report_end;
	int boost = 1;
	bool debug = config->debug;
	string ContigPrefix;
	bool Patch = 0, Fill = 0;
	int Patch_K = 5;
	int Patch_D = 30;
	int Patch_G = 2;
	double threshold = 0.2;
	string str;

	if (debug)
	{
		std::cout << "BackboneSize=" << backbone.size() << std::endl;
	}

	RefRead ref;
	ref.read_bits = (uint64_t *)malloc((size_t)(backbone.size() / 4) + 100);
	ref.alloc_sz = (size_t)(backbone.size() / 4 + 100);
	Init_Ref_Read(backbone, ref);

	int64_t bucket_count = 0, edge_cnt = 0;
	struct Backbone backbone_info;
	backbone_info.threshold = threshold;
	backbone_info.gap = 1;
	backbone_info.ContigPrefix = ContigPrefix;
	backbone_info.ScoringMethod = ScoringMethod;
	backbone_info.boost = boost;
	backbone_info.backbone = backbone;
	backbone_info.n_edges = 0;
	backbone_info.n_nodes = 0;
	backbone_info.CovTh = CovTh;
	// cout << "constructing backbone graph." << endl;
	struct Backbone backbone_info_org = backbone_info;

	SparcConsensusKmerGraphConstruction(&ref, &backbone_info_org, kmer);

	ofstream o_profile;
	if (debug)
	{
		o_profile.open("align_profile.txt");
	}

	if (debug)
	{
		std::cout << "Nodes: " << backbone_info.n_nodes << " nodes." << std::endl;

		std::cout << "Edges: " << backbone_info.n_edges << " edges." << std::endl;
		std::cout << "BackboneSize: " << backbone_info.backbone.size() << std::endl;

		std::cout << "Patch_K=" << Patch_K << std::endl;
		std::cout << "Patch_D=" << Patch_D << std::endl;
		std::cout << "Patch_G=" << Patch_G << std::endl;
		std::cout << "Patch=" << Patch << std::endl;
		std::cout << "Fill=" << Fill << std::endl;
	}

	// cout << "adding query branches." << endl;

	size_t n_reads = 0;
	backbone_info.ref_matched = 0;
	backbone_info.ref_mismatch = 0;
	backbone_info_org.ref_matched = 0;
	backbone_info_org.ref_mismatch = 0;
	struct Query query_info, query_info_org;
	query_info.Patch = Patch;
	query_info.Fill = Fill;
	query_info.Patch_K = Patch_K;
	query_info.Patch_D = Patch_D;
	query_info.Patch_G = Patch_G;
	query_info.report_b = report_begin;
	query_info.report_e = report_end;

	if (report_end > report_begin)
	{
		std::ofstream o_report_align("subalign.txt");
	}

	for (int i = 0; i < n_queries; i++)
	{
		Query *query_info = *(queries + i);
		query_info->n_exist = 0;
		query_info->n_new = 0;
		n_reads++;
		query_info->read_idx = n_reads;
		SparcAddPathToBackbone(&backbone_info_org, query_info, kmer);
	}

	// cout << "done." << endl;;

	if (debug)
	{
		std::cout << "Nodes: " << backbone_info_org.n_nodes << "." << std::endl;
		std::cout << "Edges: " << backbone_info_org.n_edges << "." << std::endl;
	}
	uint64_t cum_sum = 0;

	// coverage 值所 对应的 计数
	map<int, int> cov_cnt;

	backbone_info_org.cov_vec.resize(backbone_info_org.node_vec.size());
	if (debug)
	{

		std::cout << "backbone.size=" << backbone.size() << std::endl;
		std::cout << "backbone_info_org.node_vec.size=" << backbone_info_org.node_vec.size() << std::endl;
	}

	// 这个半径需要调整一下？ TODO
	int radius = 200;

	if (radius > backbone.size())
	{
		radius = backbone.size() - 1;
	}

	// 这个 radius 的 cov 的 count ++ 是要干什么？
	// 滑动窗口最大覆盖度值的 计算
	for (int i = 0; i < radius; ++i)
	{
		cov_cnt[backbone_info_org.node_vec[i]->cov]++;
	}

	for (int i = 0; i < backbone_info_org.node_vec.size(); ++i)
	{
		if (i >= radius)
		{
			// 上边++，这边--，不知道在搞些是什么
			cov_cnt[backbone_info_org.node_vec[i - radius]->cov]--;
			if (cov_cnt[backbone_info_org.node_vec[i - radius]->cov] <= 0)
			{
				// 会删除一些 coverage
				cov_cnt.erase(backbone_info_org.node_vec[i - radius]->cov);
			}
		}
		if (i + radius < backbone_info_org.node_vec.size())
		{
			cov_cnt[backbone_info_org.node_vec[i + radius]->cov]++;
		}

		if (cov_cnt.size() > 0)
		{
			// std::cout << "cov_cnt_max=" << cov_cnt.rbegin()->first << std::endl;
			// ;
			backbone_info_org.cov_vec[i] = cov_cnt.rbegin()->first;
		}
		else
		{
			backbone_info_org.cov_vec[i] = 0;
		}
	}

	SparcFindBestPath(&backbone_info_org);
	std::string filename = "subgraph.dot";

	if (subgraph_end > subgraph_begin && debug)
	{
		SparcOutputSubGraph(&backbone_info_org, subgraph_begin, subgraph_end, filename);
	}

	std::string consensus;

	uint64_t max_score = 0;
	size_t position = 0;

	// 这里只是遍历了 backbone 上的最大分。但是也有可能最大分不在 backbone上？
	for (int i = 0; i < backbone_info_org.node_vec.size(); ++i)
	{
		if (backbone_info_org.node_vec[i]->score > max_score)
		{
			max_score = backbone_info_org.node_vec[i]->score;
			position = i;
		}
		// backbone_info_org.node_vec[i]->in_backbone = 0;
	}

	if (max_score == 0)
	{
		char *result = (char *)malloc(backbone.size() + 1);
		strcpy(result, backbone.c_str());
		std::cout << "Empty ouput. Backbone copied." << std::endl;
		return result;
	}

	ConsensusNode *current_node = backbone_info_org.node_vec[position];

	char KmerStr[100];
	if (current_node != NULL)
	{
		uint64_t kmer_uint64 = (current_node->kmer);
		bitsarr2str(&kmer_uint64, kmer, KmerStr, 1);
		consensus = KmerStr;
		reverse(consensus.begin(), consensus.end());
		current_node = current_node->last_node;
		while (current_node != NULL)
		{
			kmer_uint64 = (current_node->kmer);
			bitsarr2str(&kmer_uint64, kmer, KmerStr, 1);

			consensus.push_back(KmerStr[0]);
			current_node = current_node->last_node;
		}

		reverse(consensus.begin(), consensus.end());
	}

	if (debug)
	{
		std::string OutputFilename2 = "DEBUG.consensus.fasta";
		std::ofstream o_cns(OutputFilename2.c_str());
		o_cns << ">Debug" << endl;
		o_cns << consensus << endl;
	}

	current_node = backbone_info_org.node_vec[position];
	// backbone_info_org.node_vec.clear();
	// backbone_info_org.node_vec.push_back(current_node);
	if (current_node != NULL)
	{
		current_node->selected = true;
		current_node = current_node->last_node;
		while (current_node != NULL)
		{
			current_node->selected = true;
			// backbone_info_org.node_vec.push_back(current_node);
			current_node = current_node->last_node;
		}
	}
	// reverse(backbone_info_org.node_vec.begin(), backbone_info_org.node_vec.end());

	// for (int ii = 0; ii < backbone_info_org.node_vec.size(); ++ii)
	// {
	// 	backbone_info_org.node_vec[ii]->in_backbone = 1;
	// 	backbone_info_org.node_vec[ii]->cns_coord = ii + 1;
	// }

	filename = "subgraph_cns.dot";
	if (subgraph_end > subgraph_begin && debug)
	{
		SparcOutputSubGraph(&backbone_info_org, subgraph_begin, subgraph_end, filename);
	}
	if (debug)
	{
		for (int i = 0; i < backbone_info_org.node_vec.size(); ++i)
		{

			o_profile << i << " " << backbone_info_org.node_vec[i]->score;

			if (backbone_info_org.node_vec[i]->coord > 0)
			{
				o_profile << " bb_coord: " << backbone_info_org.node_vec[i]->coord;
			}
			o_profile << " cns_coord: " << backbone_info_org.node_vec[i]->cns_coord << endl;
		}
	}

	std::cout << "Finished." << std::endl;

	SparcFreeInfo(&backbone_info_org);
	free(ref.read_bits);

	char *result = (char *)malloc(consensus.size() + 1);
	strcpy(result, consensus.c_str());

	return result;
}

void SparcFreeConsensusResult(char *consensus)
{
	free(consensus);
}