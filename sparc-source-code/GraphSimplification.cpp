#include <iostream>
#include <string>
#include <string.h>
#include <stdint.h>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include <fstream>
#include <cmath>
#include "time.h"
#include "BasicDataStructure.h"
#include "GraphConstruction.h"
#include "GraphSimplification.h"
using namespace std;

char SparcKmer2base(uint32_t kmer)
{
	if (kmer == 0)
	{
		return 'A';
	}
	if (kmer == 1)
	{
		return 'C';
	}

	if (kmer == 2)
	{
		return 'G';
	}

	if (kmer == 3)
	{
		return 'T';
	}

	return 'N';
}

void SparcMergeNodes(struct Backbone *backbone_info)
{
	// break the bubble links in the bfs.
	ofstream o_debug("debug_merge.txt");
	int n_merged = 0, n_replacements = 0;
	for (int i = 0; i + 1 < backbone_info->node_vec.size(); ++i)
	{
		// cout << i << " ";
		ConsensusEdgeNode *edge_ptr = backbone_info->node_vec[i]->right;
		map<uint32_t, ConsensusNode *> node_map;
		map<ConsensusNode *, int> node_cov;
		uint32_t kmer = edge_ptr->node_ptr->kmer;
		int backbone_edge_cov = backbone_info->node_vec[i]->right->edge_cov;
		// edge_ptr = edge_ptr->nxt_edge;
		int max_cov = 0;
		while (edge_ptr != NULL)
		{
			if (node_map.count(edge_ptr->node_ptr->kmer) == 0)
			{
				node_map[edge_ptr->node_ptr->kmer] = edge_ptr->node_ptr;
				node_cov[node_map[edge_ptr->node_ptr->kmer]] = edge_ptr->edge_cov;
			}
			else
			{
				node_cov[node_map[edge_ptr->node_ptr->kmer]] += edge_ptr->edge_cov;
				n_merged++;
				if (edge_ptr->node_ptr->used)
				{
					// cout << "err" << endl;
				}
				edge_ptr->node_ptr->used = 1;
			}
			edge_ptr = edge_ptr->nxt_edge;
		}

		edge_ptr = backbone_info->node_vec[i]->right;

		map<ConsensusNode *, int>::iterator it;

		for (it = node_cov.begin(); it != node_cov.end(); ++it)
		{
			if (it->first->kmer == kmer)
			{
				backbone_edge_cov = node_cov[node_map[kmer]];
			}
			else
			{
				if (max_cov < it->second)
				{
					max_cov = it->second;
				}
			}
		}

		if (max_cov >= backbone_edge_cov)
		{
			n_replacements++;

			o_debug << "position:" << i << " " << endl;
			for (it = node_cov.begin(); it != node_cov.end(); ++it)
			{
				if (it->first->kmer == kmer)
				{
					o_debug << "bb " << it->first->kmer << " " << it->second << endl;
				}
				else
				{
					o_debug << it->first->kmer << " " << it->second << endl;
					;
				}
			}
		}
	}

	cout << backbone_info->n_nodes << " nodes." << endl;
	cout << backbone_info->n_edges << " edges." << endl;
	cout << n_merged << " merges." << endl;
	cout << n_replacements << " replacements." << endl;
}

void SparcOutputPathsFromANode(ConsensusNode *begin_node, string filename, map<ConsensusNode *, bool> &Visited)
{

	ofstream o_graph(filename.c_str(), ios_base::app);
	map<uint64_t, ConsensusNode *> node_map;

	list<ConsensusNode *> node_list;
	node_list.push_back(begin_node);
	if (begin_node->coord == 14240)
	{
		cout << "";
	}
	// o_graph << "\"" << begin_node << "\"" << " [label=\"" << begin_node->kmer.kmer << "(" << begin_node->coord << ")" "\"];" << endl;
	while (node_list.size() > 0)
	{
		ConsensusNode *current_node = node_list.front();
		node_list.pop_front();
		ConsensusEdgeNode *edge_ptr = current_node->right;
		while (edge_ptr != NULL)
		{

			o_graph << "\"" << current_node << "\"" << " -> " << "\"" << edge_ptr->node_ptr << "\"";
			o_graph << " [label=\"" << edge_ptr->edge_cov << "\"];";
			o_graph << endl;

			if ((!edge_ptr->node_ptr->in_backbone))
			{
				if (Visited.count(edge_ptr->node_ptr) == 0)
				{
					Visited[edge_ptr->node_ptr] = 1;
					node_list.push_back(edge_ptr->node_ptr);
				}

				std::string color = edge_ptr->node_ptr->selected ? "color=red," : "";
				o_graph << "\"" << edge_ptr->node_ptr << "\"" << " [" << color << "label=\""
						<< SparcKmer2base(edge_ptr->node_ptr->kmer)
						<< "[" << edge_ptr->node_ptr->coord << "]"
						<< "[" << edge_ptr->node_ptr->cov << "]"
						<< "\"];"
						<< endl;
			}
			else
			{
				std::string color = edge_ptr->node_ptr->selected ? "color=red," : "";

				o_graph << "\"" << edge_ptr->node_ptr << "\"" << " [" << color << "label=\"" << SparcKmer2base(edge_ptr->node_ptr->kmer)
						<< "[" << edge_ptr->node_ptr->coord << "]"
						<< "[" << edge_ptr->node_ptr->cov << "]"
						<< ":B"
						<< "\"];"
						<< endl;
			}

			edge_ptr = edge_ptr->nxt_edge;
		}
	}

	o_graph.close();
}

void SparcOutputSubGraph(struct Backbone *backbone_info, int begin, int end, string filename)
{

	ofstream o_graph(filename.c_str());
	o_graph << "digraph G {" << endl;
	std::string color = backbone_info->node_vec[begin]->selected ? "color=red," : "";

	o_graph << "\"" << backbone_info->node_vec[begin]
			<< "\"" << " [" << color << "label=\""
			<< SparcKmer2base(backbone_info->node_vec[begin]->kmer)
			<< "[" << backbone_info->node_vec[begin]->coord << "]"
			<< "[" << backbone_info->node_vec[begin]->cov << "]"
			<< ":B"
			<< "\"];"
			<< endl;

	o_graph.close();
	uint32_t n_Rbranched = 0, n_Lbranched = 0;
	map<ConsensusNode *, bool> Visited;
	for (int i = begin; i + 1 < end; ++i)
	{
		SparcOutputPathsFromANode(backbone_info->node_vec[i], filename, Visited);
	}

	o_graph.open(filename.c_str(), ios_base::app);
	o_graph << "}" << endl;
	o_graph.close();
	o_graph.clear();
}

char *SparcMultiply(const char *a_in, const char *b_in, char *mul)
{
	char c[500], a[500], b[500];
	strcpy(a, a_in);
	strcpy(b, b_in);
	char temp[500];
	int la, lb;
	int i, j, k = 0, x = 0, y;
	long int r = 0;
	long sum = 0;
	la = strlen(a) - 1;
	lb = strlen(b) - 1;

	for (i = 0; i <= la; i++)
	{
		a[i] = a[i] - 48;
	}

	for (i = 0; i <= lb; i++)
	{
		b[i] = b[i] - 48;
	}

	for (i = lb; i >= 0; i--)
	{
		r = 0;
		for (j = la; j >= 0; j--)
		{
			temp[k++] = (b[i] * a[j] + r) % 10;
			r = (b[i] * a[j] + r) / 10;
		}
		temp[k++] = r;
		x++;
		for (y = 0; y < x; y++)
		{
			temp[k++] = 0;
		}
	}

	k = 0;
	r = 0;
	for (i = 0; i < la + lb + 2; i++)
	{
		sum = 0;
		y = 0;
		for (j = 1; j <= lb + 1; j++)
		{
			if (i <= la + j)
			{
				sum = sum + temp[y + i];
			}
			y += j + la + 1;
		}
		c[k++] = (sum + r) % 10;
		r = (sum + r) / 10;
	}
	c[k] = r;
	j = 0;
	for (i = k - 1; i >= 0; i--)
	{
		mul[j++] = c[i] + 48;
	}
	mul[j] = '\0';
	return mul;
}

void SparcBFSFindBestPath(struct Backbone *backbone_info, int node_idx)
{
	ConsensusNode *begin_node = backbone_info->node_vec[node_idx];
	list<ConsensusNode *> node_list;
	node_list.push_back(begin_node);
	int max_cov = backbone_info->cov_vec[node_idx];
	if (node_idx == 511)
	{
		cout << "";
	}
	while (node_list.size() > 0)
	{
		ConsensusNode *current_node = node_list.front();
		node_list.pop_front();
		ConsensusEdgeNode *edge_ptr = current_node->right;
		while (edge_ptr != NULL)
		{
			// if (edge_ptr->node_ptr->score < (current_node->score + edge_ptr->edge_cov))
			bool update = 0;
			if (backbone_info->ScoringMethod == 1) // 默认用的 ScoringMethod = 2
			{
				double score_org = 0.0;
				double score_new = 0.0;
				if (current_node->score == 0)
				{
					cout << "";
				}

				stringstream itoa_str1, itoa_str2, itoa_str3, itoa_str4;
				string score_str1, score_str2, degree_str1, degree_str2;

				itoa_str1 << edge_ptr->node_ptr->score;
				score_str1 = itoa_str1.str();
				itoa_str2 << edge_ptr->node_ptr->cov;
				degree_str1 = itoa_str2.str();

				itoa_str3 << current_node->score + edge_ptr->edge_cov;
				score_str2 = itoa_str3.str();
				itoa_str4 << current_node->cov + 1;
				degree_str2 = itoa_str4.str();
				char score_org_str[600], score_new_str[600];

				SparcMultiply(score_str1.c_str(), degree_str2.c_str(), score_org_str);
				SparcMultiply(score_str2.c_str(), degree_str1.c_str(), score_new_str);

				if (edge_ptr->node_ptr->score == 0)
				{
					update = 1;
				}
				else
				{
					string score_org2, score_new2;
					bool output = 0;
					for (int ii = 0; ii < strlen(score_org_str); ++ii)
					{
						if (score_org_str[ii] != '0')
						{
							output = 1;
						}
						if (output)
						{
							score_org2.push_back(score_org_str[ii]);
						}
					}
					output = 0;
					for (int ii = 0; ii < strlen(score_new_str); ++ii)
					{
						if (score_new_str[ii] != '0')
						{
							output = 1;
						}
						if (output)
						{
							score_new2.push_back(score_new_str[ii]);
						}
					}
					if (score_org2.size() < score_new2.size() || (score_org2.size() == score_new2.size() && score_org2 < score_new2))
					{
						update = 1;
					}
				}
			}

			int new_score = 0;
			if (backbone_info->ScoringMethod == 2)
			{
				if (backbone_info->threshold < 0.0) // 一般用的是 0.2
				{
					new_score = (current_node->score + max(edge_ptr->edge_cov - backbone_info->CovTh, -2));
					if (edge_ptr->node_ptr->score < new_score)
					{
						update = 1;
					}
				}
				else
				{
					int threshold = round((max_cov * backbone_info->threshold));
					if (threshold < backbone_info->CovTh) // 一般是2
					{
						threshold = backbone_info->CovTh;
					}
					new_score = (current_node->score + edge_ptr->edge_cov - threshold);
					if (edge_ptr->node_ptr->score < new_score)
					{
						update = 1;
					}
				}
			}

			// if (score_org <score_new)
			//
			if (update)
			{
				// edge_ptr->node_ptr->cov = current_node->cov + 1;
				// edge_ptr->node_ptr->score = (current_node->score + max(edge_ptr->edge_cov - backbone_info->CovTh,-2));

				edge_ptr->node_ptr->score = new_score;

				edge_ptr->node_ptr->last_node = current_node;
				if (!edge_ptr->node_ptr->in_backbone)
				{
					node_list.push_back(edge_ptr->node_ptr);
				}
			}

			edge_ptr = edge_ptr->nxt_edge;
		}
	}
}

void SparcFindBestPath(struct Backbone *backbone_info)
{

	uint32_t n_Rbranched = 0, n_Lbranched = 0;
	// backbone_info->node_vec[0]->cov = 1; // depth
	for (int i = 0; i + 1 < backbone_info->node_vec.size(); ++i)
	{

		SparcBFSFindBestPath(backbone_info, i);

		// cout << i << " ";
		ConsensusEdgeNode *edge_ptr = backbone_info->node_vec[i]->right;
		map<uint64_t, ConsensusNode *> node_map;
		map<ConsensusNode *, int> node_cov;
		int RBranch = 0;
		while (edge_ptr != NULL)
		{
			/*
			if (node_map.count(edge_ptr->node_ptr->kmer.kmer) == 0)
			{
				node_map[edge_ptr->node_ptr->kmer.kmer] = edge_ptr->node_ptr;
				node_cov[node_map[edge_ptr->node_ptr->kmer.kmer]] = edge_ptr->edge_cov;
			}
			else
			{
				node_cov[node_map[edge_ptr->node_ptr->kmer.kmer]] += edge_ptr->edge_cov;

			}
			*/
			RBranch++;
			edge_ptr = edge_ptr->nxt_edge;
		}
		if (RBranch > 1)
		{
			n_Rbranched++;
		}

		node_map.clear();
		node_cov.clear();

		edge_ptr = backbone_info->node_vec[i]->left;
		int LBranch = 0;
		while (edge_ptr != NULL)
		{
			/*
			if (node_map.count(edge_ptr->node_ptr->kmer.kmer) == 0)
			{
				node_map[edge_ptr->node_ptr->kmer.kmer] = edge_ptr->node_ptr;
				node_cov[node_map[edge_ptr->node_ptr->kmer.kmer]] = edge_ptr->edge_cov;
			}
			else
			{
				node_cov[node_map[edge_ptr->node_ptr->kmer.kmer]] += edge_ptr->edge_cov;

			}
			*/
			LBranch++;
			edge_ptr = edge_ptr->nxt_edge;
		}

		if (LBranch > 1)
		{
			n_Lbranched++;
		}
	}

	// cout << backbone_info->backbone.size() << " backbone size." << endl;
	// cout << backbone_info->n_nodes << " nodes." << endl;
	// cout << backbone_info->n_edges << " edges." << endl;
	// cout << n_Lbranched << " left branches." << endl;
	// cout << n_Rbranched << " right branches." << endl;
}

void SparcBFSFree(struct Backbone *backbone_info, int node_idx)
{
	ConsensusNode *begin_node = backbone_info->node_vec[node_idx];
	list<ConsensusNode *> node_list;
	node_list.push_back(begin_node);

	while (node_list.size() > 0)
	{
		ConsensusNode *current_node = node_list.front();
		node_list.pop_front();
		ConsensusEdgeNode *edge_ptr = current_node->right;
		while (edge_ptr != NULL)
		{
			edge_ptr->node_ptr->coord = 0;

			if (!edge_ptr->node_ptr->in_backbone)
			{
				node_list.push_back(edge_ptr->node_ptr);
			}
			ConsensusEdgeNode *nxt_edge_ptr = edge_ptr->nxt_edge;
			free(edge_ptr);

			edge_ptr = nxt_edge_ptr;
		}

		ConsensusEdgeNode *prev_edge_ptr = current_node->left;
		while (prev_edge_ptr != NULL)
		{
			ConsensusEdgeNode *edge_ptr = prev_edge_ptr->nxt_edge;
			free(prev_edge_ptr);
			prev_edge_ptr = edge_ptr;
		}
		free(current_node);
	}
}

void SparcBFSClear(struct Backbone *backbone_info, int node_idx)
{
	ConsensusNode *begin_node = backbone_info->node_vec[node_idx];
	list<ConsensusNode *> node_list;
	node_list.push_back(begin_node);
	begin_node->in_backbone = 0;
	begin_node->coord = 0;
	while (node_list.size() > 0)
	{
		ConsensusNode *current_node = node_list.front();
		node_list.pop_front();
		ConsensusEdgeNode *edge_ptr = current_node->right;
		while (edge_ptr != NULL)
		{
			edge_ptr->node_ptr->coord = 0;

			if (!edge_ptr->node_ptr->in_backbone)
			{
				node_list.push_back(edge_ptr->node_ptr);
			}

			edge_ptr = edge_ptr->nxt_edge;
		}
	}
}

void SparcClearInfo(struct Backbone *backbone_info)
{

	uint32_t n_Rbranched = 0, n_Lbranched = 0;
	backbone_info->node_vec[0]->cov = 1; // depth
	for (int i = 0; i + 1 < backbone_info->node_vec.size(); ++i)
	{

		SparcBFSClear(backbone_info, i);
	}

	cout << "Nodes information clear." << endl;
}

void SparcFreeInfo(struct Backbone *backbone_info)
{

	for (int i = 0; i + 1 < backbone_info->node_vec.size(); ++i)
	{

		SparcBFSFree(backbone_info, i);
	}

}
