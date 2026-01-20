
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

void SparcConsensusKmerGraphConstruction(struct RefRead *read, struct Backbone *backbone_info, int K_size)
{
	int readLen = read->readLen;
	int OverlappingKmers = readLen - K_size + 1;

	int Read_arr_sz = readLen / 32 + 1;
	int rem = readLen % 32;
	if (rem == 0)
	{
		Read_arr_sz--;
	}

	int Kmer_arr_sz = K_size / 32 + 1;
	int rem1 = K_size % 32;
	if (rem1 == 0)
	{
		Kmer_arr_sz--;
	}

	int tot_bits = Read_arr_sz * 64;

	// cout << "ReadLen=" << readLen << endl;
	// cout << "OverlappingKmers=" << OverlappingKmers << endl;
	// cout << "Read_arr_sz=" << Read_arr_sz << endl;

	uint64_t seq;
	// check the read to see if there is a saved kmer in the hashtable or bloom filter
	ConsensusNode *previous_node = NULL, *current_node = NULL;
	backbone_info->node_vec.resize(OverlappingKmers);

	for (int j = 0; j < OverlappingKmers; j++)
	{

		previous_node = current_node;
		get_sub_arr(read->read_bits, read->readLen, j, K_size, &seq);

		current_node = (ConsensusNode *)malloc(sizeof(ConsensusNode));
		backbone_info->node_vec[j] = current_node;
		memset(current_node, 0, sizeof(ConsensusNode));
		memcpy(&(current_node->kmer), &seq, sizeof(uint64_t));
		current_node->in_backbone = 1;
		current_node->cov++;
		backbone_info->n_nodes++;
		if (j >= 1)
		{
			// left edge,right edge
			previous_node->right = (ConsensusEdgeNode *)malloc(sizeof(ConsensusEdgeNode));
			memset(previous_node->right, 0, sizeof(ConsensusEdgeNode));
			previous_node->right->node_ptr = current_node;
			previous_node->right->edge_cov++;
			previous_node->coord = j - 1;
			backbone_info->n_edges++;
			current_node->left = (ConsensusEdgeNode *)malloc(sizeof(ConsensusEdgeNode));
			memset(current_node->left, 0, sizeof(ConsensusEdgeNode));
			current_node->left->node_ptr = previous_node;
			current_node->left->edge_cov++;
			backbone_info->n_edges++;
		}
	}
}

// 将 mismatch 拆成了两个 GAP
// gap 右对齐
void SparcNormalizeAlignment(Query *query_info)
{
	cout << "query_info->tAlignedSeq=" << query_info->tAlignedSeq << endl;
	cout << "query_info->qAlignedSeq=" << query_info->qAlignedSeq << endl;

	string qAlignedSeq_new, tAlignedSeq_new;
	size_t seq_sz = query_info->qAlignedSeq.size();
	qAlignedSeq_new.resize(2 * seq_sz);
	tAlignedSeq_new.resize(2 * seq_sz);
	int target_position = query_info->tStart;
	string target_crop, query_crop;
	int n_char = 0;
	for (int i = 0; i < seq_sz; ++i)
	{

		if (query_info->tAlignedSeq[i] != '-')
		{
			target_position++;
		}

		if (query_info->qAlignedSeq[i] == query_info->tAlignedSeq[i])
		{
			// qAlignedSeq_new.push_back(query_info->qAlignedSeq[i]);
			// tAlignedSeq_new.push_back(query_info->tAlignedSeq[i]);
			qAlignedSeq_new[n_char] = query_info->qAlignedSeq[i];
			tAlignedSeq_new[n_char] = query_info->tAlignedSeq[i];
			n_char++;
		}
		else
		{
			if (query_info->qAlignedSeq[i] != '-' && query_info->tAlignedSeq[i] != '-')
			{
				// tAlignedSeq_new.push_back(query_info->tAlignedSeq[i]);
				//  mismatch 的处理逻辑。target 先占住位置
				tAlignedSeq_new[n_char] = query_info->tAlignedSeq[i];
				n_char++;

				bool replace = 0;
				for (int j = i + 1; j < query_info->tAlignedSeq.size(); ++j)
				{
					// 找到 第一个不是 - 的 base
					// ATC   ATC
					// ACC
					if (query_info->tAlignedSeq[j] != '-')
					{
						if (query_info->tAlignedSeq[j] == query_info->qAlignedSeq[i])
						{
							replace = 1;
							query_info->tAlignedSeq[j] = '-';
							// tAlignedSeq_new.push_back(query_info->qAlignedSeq[i]);
							tAlignedSeq_new[n_char] = query_info->qAlignedSeq[i];
						}
						break;
					}
				}
				if (replace == 0)
				{
					// tAlignedSeq_new.push_back('-');
					tAlignedSeq_new[n_char] = '-';
				}
				n_char--;

				// qAlignedSeq_new.push_back('-');
				// qAlignedSeq_new.push_back(query_info->qAlignedSeq[i]);
				qAlignedSeq_new[n_char] = '-';
				n_char++;
				qAlignedSeq_new[n_char] = query_info->qAlignedSeq[i];

				n_char++;
			}
			else
			{
				if (query_info->qAlignedSeq[i] == '-')
				{

					for (int j = i + 1; j < query_info->qAlignedSeq.size(); ++j)
					{
						if (query_info->qAlignedSeq[j] != '-')
						{
							if (query_info->qAlignedSeq[j] == query_info->tAlignedSeq[i])
							{
								query_info->qAlignedSeq[i] = query_info->qAlignedSeq[j];
								query_info->qAlignedSeq[j] = '-';
							}
							break;
						}
					}
				}
				else
				{

					for (int j = i + 1; j < query_info->qAlignedSeq.size(); ++j)
					{
						if (query_info->tAlignedSeq[j] != '-')
						{
							if (query_info->tAlignedSeq[j] == query_info->qAlignedSeq[i])
							{
								query_info->tAlignedSeq[i] = query_info->tAlignedSeq[j];
								query_info->tAlignedSeq[j] = '-';
							}
							break;
						}
					}
				}

				// qAlignedSeq_new.push_back(query_info->qAlignedSeq[i]);
				// tAlignedSeq_new.push_back(query_info->tAlignedSeq[i]);
				qAlignedSeq_new[n_char] = query_info->qAlignedSeq[i];
				tAlignedSeq_new[n_char] = query_info->tAlignedSeq[i];
				n_char++;
			}
		}
	}
	qAlignedSeq_new.resize(n_char);
	tAlignedSeq_new.resize(n_char);

	cout << "tAlignedSeq_new        =" << tAlignedSeq_new << endl;
	cout << "qAlignedSeq_new        =" << qAlignedSeq_new << endl;

	query_info->qAlignedSeq = qAlignedSeq_new;
	query_info->tAlignedSeq = tAlignedSeq_new;
}

void SparcPatchGaps(Query *query_info)
{
	bool Align = 1;
	string qAlignedSeq_new, tAlignedSeq_new;
	size_t seq_sz = query_info->qAlignedSeq.size();
	qAlignedSeq_new.resize(2 * seq_sz);
	tAlignedSeq_new.resize(2 * seq_sz);
	int target_position = query_info->tStart;
	string target_crop, query_crop;
	int MinGapSize = 2, K_size = query_info->Patch_K, SearchDepth = query_info->Patch_D;
	int n_char = 0;

	for (int i = 0; i < seq_sz; ++i)
	{

		if (query_info->tAlignedSeq[i] != '-')
		{
			target_position++;
		}

		bool Patch = 0;

		if (query_info->qAlignedSeq[i] == '-')
		{
			int gap_bases = 1;
			for (int j = i + 1; j < i + MinGapSize; ++j)
			{
				if ((query_info->qAlignedSeq[j] != '-') || (j + 1 == seq_sz))
				{
					gap_bases = 0;
					break;
				}
				gap_bases++;
			}
			if (gap_bases > 1)
			{
				Patch = 1;
			}
		}

		if (!Patch && Align)
		{
			qAlignedSeq_new[n_char] = query_info->qAlignedSeq[i];
			tAlignedSeq_new[n_char] = query_info->tAlignedSeq[i];
			n_char++;
		}

		if (Patch)
		{
			int shift = -1;
			map<string, int> target_position, query_position;
			vector<int> target_index, query_index;
			target_index.resize(300);
			query_index.resize(300);
			string target_crop, query_crop;
			int target_bases = 0, query_bases = 0;
			for (int j = i; j < seq_sz; ++j)
			{
				if (query_info->tAlignedSeq[j] != '-')
				{
					target_crop.push_back(query_info->tAlignedSeq[j]);
					target_index[target_bases] = j;
					target_bases++;
				}
				if (target_bases >= SearchDepth + K_size)
				{
					break;
				}
			}

			for (int j = i; j < seq_sz; ++j)
			{
				if (query_info->qAlignedSeq[j] != '-')
				{
					query_crop.push_back(query_info->qAlignedSeq[j]);
					query_index[query_bases] = j;
					query_bases++;
				}
				if (query_bases >= SearchDepth + K_size)
				{
					break;
				}
			}

			if (!Align)
			{
				shift = MinGapSize;
				for (int j = 0; j + K_size < target_crop.size(); ++j)
				{
					string kmer = target_crop.substr(j, K_size);
					if (target_position.count(kmer) == 0)
					{
						target_position[kmer] = j;
					}
				}
				int target_matched = -1, query_matched = -1;
				for (int j = 0; j + K_size < query_crop.size(); ++j)
				{
					string kmer = query_crop.substr(j, K_size);
					if (target_position.count(kmer))
					{
						target_matched = target_position[kmer];
						query_matched = j;
						break;
					}
				}

				if (target_matched >= 0)
				{
					bool Debug = 0;
					if (Debug)
					{
						cout << "before: " << endl;
						cout << query_info->tAlignedSeq.substr(i, 100) << endl;
						cout << query_info->qAlignedSeq.substr(i, 100) << endl;
					}
					// clear the query bases
					for (int j = 0; j < query_matched + K_size; ++j)
					{
						query_info->qAlignedSeq[query_index[j]] = '-';
					}
					for (int j = 0; j < K_size; ++j)
					{
						query_info->qAlignedSeq[target_index[target_matched + j]] = query_info->tAlignedSeq[target_index[target_matched + j]];
					}

					if (Debug)
					{
						cout << "after: " << endl;
						cout << query_info->tAlignedSeq.substr(i, 100) << endl;
						cout << query_info->qAlignedSeq.substr(i, 100) << endl;
					}
				}
			}
			else
			{

				int match = 100, mismatch = -100, gap_cost = -10, band_width = 50;
				string A_aln, B_aln;
				struct aln_t aln_t;
				int score = 0;
				char qry_char[300];
				char ref_char[300];
				strcpy(ref_char, target_crop.c_str());
				strcpy(qry_char, query_crop.c_str());

				score = GlobalAlign(&aln_t, qry_char, ref_char, match, mismatch, gap_cost, band_width);
				PrintAlign(&aln_t, qry_char, ref_char, A_aln, B_aln);

				bool Debug = 0;

				if (Debug)
				{
					cout << i << endl;

					if (i == 27)
					{
						cout << "";
					}
					cout << "before: " << endl;
					cout << query_info->tAlignedSeq.substr(i, 70) << endl;
					cout << query_info->qAlignedSeq.substr(i, 70) << endl;
				}
				if (Debug)
				{
					cout << "after: " << endl;

					cout << B_aln << endl;
					cout << A_aln << endl;
				}
				int A_cnt = 0, B_cnt = 0;
				for (int jj = 0; jj != A_aln.size(); ++jj)
				{
					if (A_aln[jj] != '-')
					{
						A_cnt++;
					}
				}

				for (int jj = 0; jj != B_aln.size(); ++jj)
				{
					if (B_aln[jj] != '-')
					{
						B_cnt++;
					}
				}

				int cnt = 0;
				int match_position = -1;
				for (int jj = (int)A_aln.size() - 1; jj >= 0; --jj)
				{
					if (A_aln[jj] == B_aln[jj])
					{
						cnt++;
					}
					else
					{
						cnt = 0;
					}
					if (cnt == K_size)
					{
						match_position = jj + K_size - 1;
						break;
					}
				}
				// append align

				if (match_position == -1)
				{
					shift = MinGapSize - 1;
					qAlignedSeq_new[n_char] = query_info->qAlignedSeq[i];
					tAlignedSeq_new[n_char] = query_info->tAlignedSeq[i];
					n_char++;

					i = i + shift - 1;
					continue;
				}
				for (int jj = 0; jj <= match_position; ++jj)
				{
					qAlignedSeq_new[n_char] = A_aln[jj];
					tAlignedSeq_new[n_char] = B_aln[jj];
					n_char++;
				}

				int target_shift = 0, query_shift = 0;
				for (int jj = 0; jj <= match_position; ++jj)
				{
					if (B_aln[jj] != '-')
					{
						target_shift++;
					}
					if (A_aln[jj] != '-')
					{
						query_shift++;
					}
				}

				// clear the query bases
				int q_base = 0;
				for (int j = i; j < seq_sz; ++j)
				{
					if (query_info->qAlignedSeq[j] != '-')
					{
						q_base++;
						query_info->qAlignedSeq[j] = '-';
					}
					if (q_base >= query_shift)
					{
						break;
					}
				}

				int t_base = 0;
				for (int j = i; j < seq_sz; ++j)
				{

					if (query_info->tAlignedSeq[j] != '-')
					{
						t_base++;
						query_info->tAlignedSeq[j] = '-';
					}
					if (t_base >= target_shift)
					{
						shift = j - i;
						break;
					}
				}
			}

			i = i + shift - 1;
		}
	}

	if (!Align)
	{
		for (int i = 0; i < seq_sz; ++i)
		{
			if (query_info->tAlignedSeq[i] != '-' || query_info->qAlignedSeq[i] != '-')
			{
				qAlignedSeq_new[n_char] = query_info->qAlignedSeq[i];
				tAlignedSeq_new[n_char] = query_info->tAlignedSeq[i];
				n_char++;
			}
		}
	}

	qAlignedSeq_new.resize(n_char);
	tAlignedSeq_new.resize(n_char);
	query_info->qAlignedSeq = qAlignedSeq_new;

	query_info->tAlignedSeq = tAlignedSeq_new;
}

void SparcFillGaps(Query *query_info)
{
	bool Align = 1;
	string qAlignedSeq_new, tAlignedSeq_new;
	size_t seq_sz = query_info->qAlignedSeq.size();
	qAlignedSeq_new.resize(2 * seq_sz);
	tAlignedSeq_new.resize(2 * seq_sz);
	int target_position = query_info->tStart;
	string target_crop, query_crop;
	int MinGapSize = query_info->Patch_G, K_size = query_info->Patch_K, SearchDepth = query_info->Patch_D;
	int n_char = 0;

	for (int i = 0; i < seq_sz; ++i)
	{

		if (query_info->tAlignedSeq[i] != '-')
		{
			target_position++;
		}

		bool Patch = 0;
		// detect gaps in query
		if (query_info->qAlignedSeq[i] == '-')
		{
			int gap_bases = 1;
			for (int j = i + 1; j < i + MinGapSize; ++j)
			{
				if ((query_info->qAlignedSeq[j] != '-') || (j + 1 == seq_sz))
				{
					gap_bases = 0;
					break;
				}
				gap_bases++;
			}
			if (gap_bases > 1)
			{
				Patch = 1;
			}
		}

		// detect gaps in target
		if (query_info->tAlignedSeq[i] == '-')
		{
			int gap_bases = 1;
			for (int j = i + 1; j < i + MinGapSize; ++j)
			{
				if ((query_info->tAlignedSeq[j] != '-') || (j + 1 == seq_sz))
				{
					gap_bases = 0;
					break;
				}
				gap_bases++;
			}
			if (gap_bases > 1)
			{
				Patch = 1;
			}
		}

		if (!Patch)
		{
			qAlignedSeq_new[n_char] = query_info->qAlignedSeq[i];
			tAlignedSeq_new[n_char] = query_info->tAlignedSeq[i];
			n_char++;
			continue;
		}

		int shift = -1;
		vector<int> target_index, query_index;
		target_index.resize(300);
		query_index.resize(300);
		string target_crop, query_crop;
		target_crop.resize(300);
		query_crop.resize(300);
		int target_bases = 0, query_bases = 0;
		for (int j = i; j < seq_sz; ++j)
		{
			if (query_info->tAlignedSeq[j] != '-')
			{
				target_crop[target_bases] = query_info->tAlignedSeq[j];
				target_index[target_bases] = j;
				target_bases++;
			}
			if (target_bases >= SearchDepth + K_size)
			{
				break;
			}
		}

		for (int j = i; j < seq_sz; ++j)
		{
			if (query_info->qAlignedSeq[j] != '-')
			{
				query_crop[query_bases] = query_info->qAlignedSeq[j];
				query_index[query_bases] = j;
				query_bases++;
			}
			if (query_bases >= SearchDepth + K_size)
			{
				break;
			}
		}
		target_crop.resize(target_bases);
		query_crop.resize(query_bases);

		int match = 100, mismatch = -100, gap_cost = -10, band_width = 50;
		string Query_aln, Target_aln;
		struct aln_t aln_t;
		int score = 0;
		char qry_char[300];
		char ref_char[300];
		strcpy(ref_char, target_crop.c_str());
		strcpy(qry_char, query_crop.c_str());

		score = GlobalAlign(&aln_t, qry_char, ref_char, match, mismatch, gap_cost, band_width);
		PrintAlign(&aln_t, qry_char, ref_char, Query_aln, Target_aln);

		bool Debug = 0;
		if (Debug)
		{
			cout << i << endl;

			cout << "before: " << endl;
			cout << query_info->tAlignedSeq.substr(i, 70) << endl;
			cout << query_info->qAlignedSeq.substr(i, 70) << endl;
		}
		if (Debug)
		{
			cout << "after: " << endl;

			cout << Target_aln << endl;
			cout << Query_aln << endl;
		}

		int cnt = 0;
		int match_position = -1;
		for (int jj = (int)Query_aln.size() - 1; jj >= 0; --jj)
		{
			if (Query_aln[jj] == Target_aln[jj])
			{
				cnt++;
			}
			else
			{
				cnt = 0;
			}
			if (cnt == K_size)
			{
				match_position = jj + K_size - 1;
				break;
			}
		}
		// append align

		int shift_end;
		shift_end = i + 1;
		if (query_crop.size() == target_crop.size() && target_crop.size() == SearchDepth + K_size)
		{
			shift_end = i + SearchDepth; // min(target_index[SearchDepth + K_size - 1], query_index[SearchDepth + K_size - 1]);
		}

		if (match_position == -1)
		{
			int j;
			for (j = i; j < shift_end; ++j)
			{
				qAlignedSeq_new[n_char] = query_info->qAlignedSeq[j];
				tAlignedSeq_new[n_char] = query_info->tAlignedSeq[j];
				n_char++;
			}
			i = j - 1;

			continue;
		}

		for (int jj = 0; jj <= match_position; ++jj)
		{
			qAlignedSeq_new[n_char] = Query_aln[jj];
			tAlignedSeq_new[n_char] = Target_aln[jj];
			n_char++;
		}

		int target_shift = 0, query_shift = 0;
		for (int jj = 0; jj <= match_position; ++jj)
		{
			if (Target_aln[jj] != '-')
			{
				target_shift++;
			}
			if (Query_aln[jj] != '-')
			{
				query_shift++;
			}
		}

		// clear the aligned bases
		int q_base = 0;
		for (int j = i; j < seq_sz; ++j)
		{
			if (query_info->qAlignedSeq[j] != '-')
			{
				q_base++;
				query_info->qAlignedSeq[j] = '-';
			}
			if (q_base >= query_shift)
			{
				shift = query_index[q_base - 1];
				break;
			}
		}

		int t_base = 0;
		for (int j = i; j < seq_sz; ++j)
		{

			if (query_info->tAlignedSeq[j] != '-')
			{
				t_base++;
				query_info->tAlignedSeq[j] = '-';
			}
			if (t_base >= target_shift)
			{
				shift = min(shift, target_index[t_base - 1]);
				break;
			}
		}

		i = shift;
	}

	qAlignedSeq_new.resize(n_char);
	tAlignedSeq_new.resize(n_char);
	query_info->qAlignedSeq = qAlignedSeq_new;
	query_info->tAlignedSeq = tAlignedSeq_new;
}

void SparcAddPathToBackbone(struct Backbone *backbone_info, struct Query *query_info, int K_size)
{
	cout << "Add_Path_To_Backbone" << endl;
	ofstream o_report_align;
	bool boost_edges = 0;
	string ContigPrefix = backbone_info->ContigPrefix;
	cout << "ContigPrefix=" << ContigPrefix << endl;
	if (ContigPrefix.size() > 0)
	{
		if (query_info->qName.substr(0, ContigPrefix.size()) == ContigPrefix)
		{
			boost_edges = 1;
		}
	}
	bool DEBUG = 1;
	cout << "query_info->tStrand=" << query_info->tStrand << endl;

	if (query_info->tStrand == '-')
	{
		reverse_complement_str(query_info->qAlignedSeq);
		reverse_complement_str(query_info->tAlignedSeq);
		reverse(query_info->matchPattern.begin(), query_info->matchPattern.end());
	}

	SparcNormalizeAlignment(query_info);
	if (query_info->Patch)
	{
		SparcPatchGaps(query_info);
	}
	else
	{
		if (query_info->Fill)
		{
			SparcFillGaps(query_info);
		}
	}
	//		FillGaps(query_info);

	for (int i = 0; i < query_info->tAlignedSeq.size(); ++i)
	{
		if (query_info->tAlignedSeq[i] != '-')
		{
			if (query_info->qAlignedSeq[i] != '-')
			{
				backbone_info->ref_matched++;
			}
			else
			{
				backbone_info->ref_mismatch++;
			}
		}
	}

	if ((query_info->report_e > query_info->report_b) && (query_info->tStart < query_info->report_b) && (query_info->tEnd > query_info->report_e))
	{
		o_report_align.open("subalign.txt", ios_base::app);
		int target_position = query_info->tStart;
		string target_crop, query_crop;
		for (int k = 0; k < query_info->tAlignedSeq.size(); ++k)
		{
			if (query_info->tAlignedSeq[k] != '-')
			{
				target_position++;
			}
			if (target_position > query_info->report_b)
			{
				target_crop.push_back(query_info->tAlignedSeq[k]);
				query_crop.push_back(query_info->qAlignedSeq[k]);
				//
			}

			if (target_position > query_info->report_e)
			{
				o_report_align << ">target_" << query_info->read_idx << endl
							   << target_crop << endl; //
				o_report_align << ">query_" << query_info->read_idx << endl
							   << query_crop << endl; //
				break;
			}
		}
	}

	// string backbone_seg = backbone_info->backbone.substr(query_info->tStart, query_info->tEnd - query_info->tStart);
	// cout << backbone_seg << endl;
	// cout << target_seg << endl;
	int target_position = query_info->tStart;
	ConsensusNode *previous_node = NULL, *current_node = NULL;
	int MatchPosition = -100, PreviousMatchPosition = -100;

	// 开始构图
	cout << "Start Build Graph" << endl;

	for (int i = 0; i + K_size <= query_info->qAlignedSeq.size(); ++i)
	{

		// 如果是 deletion 的话, 调整一下 target position。就直接 continue 了
		// 没有新的 kmer 产生，cov 也不会发生变化
		if (query_info->qAlignedSeq[i] == '-')
		{
			MatchPosition = -1;
			PreviousMatchPosition = -1;
			if (query_info->tAlignedSeq[i] != '-')
			{
				target_position++;
			}
			continue;
		}

		// 如果是 insertion, match, mismatch 会继续
		PreviousMatchPosition = MatchPosition;
		MatchPosition = target_position;

		// 判断当前position 是不是 match position
		for (int j = 0; j < K_size; ++j)
		{
			if (query_info->qAlignedSeq[i + j] != query_info->tAlignedSeq[i + j])
			{
				MatchPosition = -1;
				break;
			}
		}

		// 如果是 insertion, match, mismatch 会继续
		string KmerStr;
		KmerStr.resize(K_size);
		int KmerSize = 0;
		// 确定当前位置， query 的 kmer 是啥

		// 在这个位置, 可能是 insertion 区域、match 区域。
		for (int j = i; j < query_info->qAlignedSeq.size(); ++j)
		{
			if (query_info->qAlignedSeq[j] != '-')
			{
				KmerStr[KmerSize] = query_info->qAlignedSeq[j];
				KmerSize++;
			}

			if (KmerSize == K_size)
			{
				break;
			}
		}

		if (KmerSize != K_size)
		{
			break;
		}
		uint64_t seq = 0;
		str2bitsarr(KmerStr.c_str(), K_size, &seq, 1);

		// 如果是 match。 这里的 match 是 equal
		// match 是对应了 backbone 的 match。 在 backbone 是一定能找到 锚点的
		if (MatchPosition >= 0)
		{

			query_info->n_exist++;
			current_node = backbone_info->node_vec[MatchPosition];

			if (current_node != NULL) // TODO: 为啥要判断 previous_node != NULL ????
			{
				current_node->cov++;
			}

			// link previous node to the current one

			if (previous_node != NULL)
			{

				ConsensusEdgeNode **edge_p2p = &(previous_node->right);

				while (*edge_p2p != NULL)
				{
					// 这个 current node 一定是 backbone 上的 node
					// 如果 previous_node 在 backbone 上 ？
					// 如果 previsou_node 不在 backnone 上？
					if ((*edge_p2p)->node_ptr == current_node)
					{
						(*edge_p2p)->edge_cov++;
						if (boost_edges)
						{
							(*edge_p2p)->edge_cov += backbone_info->boost;
						}

						break;
					}

					edge_p2p = &((*edge_p2p)->nxt_edge);
				}

				if (*edge_p2p == NULL)
				{

					(*edge_p2p) = (ConsensusEdgeNode *)malloc(sizeof(ConsensusEdgeNode));
					memset((*edge_p2p), 0, sizeof(ConsensusEdgeNode));
					(*edge_p2p)->nxt_edge = NULL;
					(*edge_p2p)->node_ptr = current_node;
					(*edge_p2p)->edge_cov++;
					if (boost_edges)
					{
						(*edge_p2p)->edge_cov += backbone_info->boost;
					}
					backbone_info->n_edges++;
				}

				edge_p2p = &(current_node->left);
				while (*edge_p2p != NULL)
				{
					if ((*edge_p2p)->node_ptr == previous_node)
					{
						(*edge_p2p)->edge_cov++;
						if (boost_edges)
						{
							(*edge_p2p)->edge_cov += backbone_info->boost;
						}

						break;
					}
					edge_p2p = &((*edge_p2p)->nxt_edge);
				}
				if (*edge_p2p == NULL)
				{
					(*edge_p2p) = (ConsensusEdgeNode *)malloc(sizeof(ConsensusEdgeNode));
					memset((*edge_p2p), 0, sizeof(ConsensusEdgeNode));
					(*edge_p2p)->nxt_edge = NULL;
					(*edge_p2p)->node_ptr = previous_node;
					(*edge_p2p)->edge_cov++;
					if (boost_edges)
					{
						(*edge_p2p)->edge_cov += backbone_info->boost;
					}
					backbone_info->n_edges++;
				}
			}
		}
		else
		{
			// 在这个位置, 可能是 insertion 区域、mismatch 区域。
			// link previous node to the current one

			if (previous_node == NULL)
			{
				current_node = (ConsensusNode *)malloc(sizeof(ConsensusNode));
				memset(current_node, 0, sizeof(ConsensusNode));

				current_node->kmer = (uint32_t)seq;
				current_node->cov++;
				backbone_info->n_nodes++;
				previous_node = current_node;
				query_info->n_new++;
				// cout << i << ", ";
			}
			else
			{

				ConsensusEdgeNode **edge_p2p = &(previous_node->right);
				while (*edge_p2p != NULL)
				{
					// 还不明白为什么会有 in_backbone == 0 的判断。
					// 因为这里处理的不是 equal，所以一定不能指向 backbone
					if ((*edge_p2p)->node_ptr->kmer == (uint32_t)seq && (*edge_p2p)->node_ptr->in_backbone == 0)
					{
						(*edge_p2p)->edge_cov++;
						if (boost_edges)
						{
							(*edge_p2p)->edge_cov += backbone_info->boost;
						}
						current_node = (*edge_p2p)->node_ptr;
						current_node->cov++;

						break;
					}
					edge_p2p = &((*edge_p2p)->nxt_edge);
				}

				if (*edge_p2p == NULL) // 创建一个边指向新的 kmer
				{
					(*edge_p2p) = (ConsensusEdgeNode *)malloc(sizeof(ConsensusEdgeNode));
					memset((*edge_p2p), 0, sizeof(ConsensusEdgeNode));
					(*edge_p2p)->nxt_edge = NULL;
					current_node = (ConsensusNode *)malloc(sizeof(ConsensusNode));
					query_info->n_new++;
					memset(current_node, 0, sizeof(ConsensusNode));
					current_node->kmer = (uint32_t)seq;
					current_node->cov++;
					(*edge_p2p)->node_ptr = current_node;

					(*edge_p2p)->edge_cov++;
					if (boost_edges)
					{
						(*edge_p2p)->edge_cov += backbone_info->boost;
					}
					backbone_info->n_edges++;
					backbone_info->n_nodes++;
					// cout << i <<", ";
				}
				// link current node to previous,
				edge_p2p = &(current_node->left);
				while (*edge_p2p != NULL)
				{
					if ((*edge_p2p)->node_ptr->kmer == previous_node->kmer) // 这里直接判断 node_ptr 不就完了？为啥还要判断 kmer
					{
						(*edge_p2p)->edge_cov++;
						if (boost_edges)
						{
							(*edge_p2p)->edge_cov += backbone_info->boost;
						}

						break;
					}
					edge_p2p = &((*edge_p2p)->nxt_edge);
				}
				if (*edge_p2p == NULL) // 创建一个边指向来时路
				{
					(*edge_p2p) = (ConsensusEdgeNode *)malloc(sizeof(ConsensusEdgeNode));
					memset((*edge_p2p), 0, sizeof(ConsensusEdgeNode));
					(*edge_p2p)->nxt_edge = NULL;
					(*edge_p2p)->node_ptr = previous_node;
					(*edge_p2p)->edge_cov++;
					if (boost_edges)
					{
						(*edge_p2p)->edge_cov += backbone_info->boost;
					}
					backbone_info->n_edges++;
				}
			}
		}

		previous_node = current_node;

		if (query_info->tAlignedSeq[i] != '-')
		{
			target_position++;
		}
	}
}
