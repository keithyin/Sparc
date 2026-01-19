
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
using namespace std;
// typedef unsigned __int64 uint64_t;
// typedef __int64 int64_t;
// typedef unsigned __int8 uint8_t;
// typedef unsigned __int16 uint16_t;
// typedef unsigned __int32 uint32_t;
// typedef __int32 int32_t;

// These are the structures to save the k-mers.32 bases, 64 bases, 96 bases, 128 bases.

// path information in the BFS bubble removal

Query *NewQuery()
{
	Query *query = new Query();
	query->Fill = false;
	query->Patch = false;

	return query;
}
void FreeQuery(Query *query)
{
	if (query != nullptr)
	{
		delete query;
	}
};

void QuerySetQueryName(Query &query, const char *name)
{
	query.qName = name ? name : "";
}

void QuerySetTargetName(Query &query, const char *name)
{
	query.tName = name ? name : "";
}

void QuerySetQueryAlignedSeq(Query &query, const char *seq)
{
	query.qAlignedSeq = seq ? seq : "";
}

void QuerySetMatchPattern(Query &query, const char *pattern)
{
	query.matchPattern = pattern ? pattern : "";
}

void QuerySetTargetAlignedSeq(Query &query, const char *seq)
{
	query.tAlignedSeq = seq ? seq : "";
}

/* ---------- basic attributes ---------- */
void QuerySetQueryLength(Query &query, int queryLength)
{
	query.qLength = queryLength;
}

void QuerySetQueryStart(Query &query, int qStart)
{
	query.qStart = qStart;
}

void QuerySetQueryEnd(Query &query, int qEnd)
{
	query.qEnd = qEnd;
}

void QuerySetTargetLength(Query &query, int targetLength)
{
	query.tLength = targetLength;
}

void QuerySetTargetStart(Query &query, int tStart)
{
	query.tStart = tStart;
}

void QuerySetTargetEnd(Query &query, int tEnd)
{
	query.tEnd = tEnd;
}

/* ---------- alignment stats ---------- */
void QuerySetScore(Query &query, int score)
{
	query.score = score;
}

void QuerySetNumMatch(Query &query, int numMatch)
{
	query.numMatch = numMatch;
}

void QuerySetNumMismatch(Query &query, int numMismatch)
{
	query.numMismatch = numMismatch;
}

void QuerySetNumIns(Query &query, int numIns)
{
	query.numIns = numIns;
}

void QuerySetNumDel(Query &query, int numDel)
{
	query.numDel = numDel;
}

void QuerySetMapQV(Query &query, int mapQV)
{
	query.mapQV = mapQV;
}

/* ---------- strand / index ---------- */
void QuerySetQueryStrand(Query &query, char strand)
{
	query.qStrand = strand;
}

void QuerySetTargetStrand(Query &query, char strand)
{
	query.tStrand = strand;
}

void QuerySetReadIndex(Query &query, size_t read_idx)
{
	query.read_idx = read_idx;
}

/* ---------- report range ---------- */
void QuerySetReportBegin(Query &query, int report_b)
{
	query.report_b = report_b;
}

void QuerySetReportEnd(Query &query, int report_e)
{
	query.report_e = report_e;
}

/* ---------- counters ---------- */
void QuerySetNumExist(Query &query, int n_exist)
{
	query.n_exist = n_exist;
}

void QuerySetNumNew(Query &query, int n_new)
{
	query.n_new = n_new;
}

/* ---------- patch flags ---------- */
void QuerySetPatch(Query &query, bool patch)
{
	query.Patch = patch;
}

void QuerySetFill(Query &query, bool fill)
{
	query.Fill = fill;
}

/* ---------- patch parameters ---------- */
void QuerySetPatchK(Query &query, int k)
{
	query.Patch_K = k;
}

void QuerySetPatchD(Query &query, int d)
{
	query.Patch_D = d;
}

void QuerySetPatchG(Query &query, int g)
{
	query.Patch_G = g;
}

bool get_a_fasta_read(ifstream &fasta_in, string &tag, string &str, string &n_tag)
{

	ifstream tmp_ifstream;
	string temp;
	if (!getline(fasta_in, temp))
	{
		return 0;
	}
	if (temp[temp.size() - 1] == '\n' || temp[temp.size() - 1] == '\r')
	{
		temp.resize(temp.size() - 1);
	}

	str.clear();
	if (temp[0] == '>')
	{
		tag = temp;
	}
	else
	{
		tag = n_tag;
		str = temp;
	}

	while (getline(fasta_in, temp))
	{

		if (temp[temp.size() - 1] == '\n' || temp[temp.size() - 1] == '\r')
		{
			temp.resize(temp.size() - 1);
		}

		if ((temp.size() > 0 && (temp[0] == '>' || temp[0] == '\n' || temp[0] == '\r')))
		{
			n_tag = temp;
			return 1;
		}
		else
		{
			str += temp;
		}
	}
	return 1;
}

bool get_a_fastq_read(ifstream &fastq_in, string &tag, string &seq, string &quality)
{

	ifstream tmp_ifstream;
	string temp;
	if (!getline(fastq_in, temp))
	{
		return 0;
	}
	seq.clear();
	if (temp[0] == '@')
	{
		tag = temp;
	}
	else
	{
		return 0;
	}
	getline(fastq_in, seq); // seq
	if (seq[seq.size() - 1] == '\n' || seq[seq.size() - 1] == '\r')
	{
		seq.resize(seq.size() - 1);
	}
	getline(fastq_in, temp); //'+'
	getline(fastq_in, quality);
	if (quality[quality.size() - 1] == '\n' || quality[quality.size() - 1] == '\r')
	{
		quality.resize(quality.size() - 1);
	}

	return 1;
}

uint64_t *str2bitsarr(const char *c_str, int len, uint64_t *b_str, int arr_sz)
{

	for (int k = 0; k < arr_sz; ++k)
	{
		b_str[k] = 0;
	}
	int arr_sz_needed = len / 32 + 1;
	int rem = len % 32;
	if (rem == 0)
	{
		arr_sz_needed--;
	}

	int beg_arr_idx = arr_sz - arr_sz_needed;
	if (rem == 0 && arr_sz_needed > 0)
	{
		rem = 32;
	}
	for (int k = 0; k < len; k++)
	{
		if (rem == 0)
		{
			beg_arr_idx++;
			rem = 32;
		}

		switch (c_str[k])
		{
		case ('A'):
		case ('a'):
		case ('0'):
			b_str[beg_arr_idx] <<= 2; // L_shift_NB(b_str,2,arr_sz);
			rem--;
			// b_str<<=2;
			break;

		case ('C'):
		case ('c'):
		case ('1'):
			b_str[beg_arr_idx] <<= 2; // L_shift_NB(b_str,2,arr_sz);
			++b_str[beg_arr_idx];	  //++b_str[arr_sz-1];
			rem--;
			//			++(b_str<<=2);
			break;

		case 'G':
		case 'g':
		case '2':
			b_str[beg_arr_idx] <<= 1; // L_shift_NB(b_str,1,arr_sz);
			++b_str[beg_arr_idx];	  //++b_str[arr_sz-1];
			b_str[beg_arr_idx] <<= 1; // L_shift_NB(b_str,1,arr_sz);
			rem--;					  //(++(b_str<<=1))<<=1;
			break;

		case 'T':
		case 't':
		case '3':
			b_str[beg_arr_idx] <<= 1; // L_shift_NB(b_str,1,arr_sz);
			++b_str[beg_arr_idx];
			b_str[beg_arr_idx] <<= 1; // L_shift_NB(b_str,1,arr_sz);
			++b_str[beg_arr_idx];
			rem--;
			//++((++(b_str<<=1))<<=1);
			break;
		default:
			return b_str;
		}

		//	cout<<b_str<<endl;
	}
	return b_str;
}

// hash functions
uint64_t MurmurHash64A(const void *key, int len, unsigned int seed)
{
	const uint64_t m = 0xc6a4a7935bd1e995;
	const int r = 47;

	uint64_t h = seed ^ (len * m);

	const uint64_t *data = (const uint64_t *)key;
	const uint64_t *end = data + (len / 8);

	while (data != end)
	{
		uint64_t k = *data++;

		k *= m;
		k ^= k >> r;
		k *= m;

		h ^= k;
		h *= m;
	}

	const unsigned char *data2 = (const unsigned char *)data;

	switch (len & 7)
	{
	case 7:
		h ^= uint64_t(data2[6]) << 48;
	case 6:
		h ^= uint64_t(data2[5]) << 40;
	case 5:
		h ^= uint64_t(data2[4]) << 32;
	case 4:
		h ^= uint64_t(data2[3]) << 24;
	case 3:
		h ^= uint64_t(data2[2]) << 16;
	case 2:
		h ^= uint64_t(data2[1]) << 8;
	case 1:
		h ^= uint64_t(data2[0]);
		h *= m;
	};

	h ^= h >> r;
	h *= m;
	h ^= h >> r;

	return h;
}

uint64_t MurmurHash64B(const void *key, int len, unsigned int seed)
{
	const unsigned int m = 0x5bd1e995;
	const int r = 24;

	unsigned int h1 = seed ^ len;
	unsigned int h2 = 0;

	const unsigned int *data = (const unsigned int *)key;

	while (len >= 8)
	{
		unsigned int k1 = *data++;
		k1 *= m;
		k1 ^= k1 >> r;
		k1 *= m;
		h1 *= m;
		h1 ^= k1;
		len -= 4;

		unsigned int k2 = *data++;
		k2 *= m;
		k2 ^= k2 >> r;
		k2 *= m;
		h2 *= m;
		h2 ^= k2;
		len -= 4;
	}

	if (len >= 4)
	{
		unsigned int k1 = *data++;
		k1 *= m;
		k1 ^= k1 >> r;
		k1 *= m;
		h1 *= m;
		h1 ^= k1;
		len -= 4;
	}

	switch (len)
	{
	case 3:
		h2 ^= ((unsigned char *)data)[2] << 16;
	case 2:
		h2 ^= ((unsigned char *)data)[1] << 8;
	case 1:
		h2 ^= ((unsigned char *)data)[0];
		h2 *= m;
	};

	h1 ^= h2 >> 18;
	h1 *= m;
	h2 ^= h1 >> 22;
	h2 *= m;
	h1 ^= h2 >> 17;
	h1 *= m;
	h2 ^= h1 >> 19;
	h2 *= m;

	uint64_t h = h1;

	h = (h << 32) | h2;

	return h;
}

// convert a string of nucleotide bases into bit array.
void Init_Read(string &seq, struct Read &read)
{
	read.readLen = strlen(seq.c_str());
	int Read_arr_sz = read.readLen / 32 + 1;
	int rem = read.readLen % 32;
	if (rem == 0)
	{
		Read_arr_sz--;
	}

	str2bitsarr(seq.c_str(), (int)seq.size(), read.read_bits, Read_arr_sz);
}

void Init_Ref_Read(string &seq, struct RefRead &read)
{
	read.readLen = strlen(seq.c_str());
	int Read_arr_sz = read.readLen / 32 + 1;
	int rem = read.readLen % 32;
	if (rem == 0)
	{
		Read_arr_sz--;
	}

	str2bitsarr(seq.c_str(), (int)seq.size(), read.read_bits, Read_arr_sz);
}

void reverse_complement_str(string &str)
{
	if (str.size() == 0)
	{
		return;
	}

	reverse(str.begin(), str.end());
	for (size_t i = 0; i != str.size(); ++i)
	{
		switch (str[i])
		{
		case 'A':
		case 'a':
			str[i] = 'T';
			break;
		case 'C':
		case 'c':
			str[i] = 'G';
			break;
		case 'G':
		case 'g':
			str[i] = 'C';
			break;
		case 'T':
		case 't':
			str[i] = 'A';
			break;
		case 'N':
		case 'n':
			str[i] = 'N';
			break;
		case '-':
			str[i] = '-';
			break;
		default:
			cout << "error: complement_str" << str[i] << endl;
			return;
		}
	}
}

// get the complement of a string of nucleotide bases
void complement_str(string &str)
{
	for (size_t i = 0; i != str.size(); ++i)
	{
		switch (str[i])
		{
		case 'A':
		case 'a':
			str[i] = 'T';
			break;
		case 'C':
		case 'c':
			str[i] = 'G';
			break;
		case 'G':
		case 'g':
			str[i] = 'C';
			break;
		case 'T':
		case 't':
			str[i] = 'A';
			break;
		case 'N':
		case 'n':
			str[i] = 'N';
			break;
		case '-':
			str[i] = '-';
			break;
		default:
			cout << "error: complement_str" << str[i] << endl;
			return;
		}
	}
}

// convert a bit array into a string
char *bitsarr2str(uint64_t *b_seq, int len, char *c_str, int arr_sz)
{

	int tot_bits = arr_sz * 64;
	// char *c_str;
	// c_str=(char*) malloc(sizeof(char)*(len+1));
	// #pragma omp parallel for
	for (int i = 0; i < len; ++i)
	{
		uint64_t temp, temp2[100]; /////////////////////////
		for (int k = 0; k < arr_sz; ++k)
		{
			temp2[k] = b_seq[k];
		}
		L_shift_NB(temp2, tot_bits - (len - i) * 2, arr_sz);
		R_shift_NB(temp2, tot_bits - 2, arr_sz);
		// uint64_t temp=(b_seq<<(64-(len-i)*2))>>62;
		temp = temp2[arr_sz - 1];
		switch (temp)
		{
		case 0:
			c_str[i] = 'A';
			break;
		case 1:
			c_str[i] = 'C';
			break;
		case 2:
			c_str[i] = 'G';
			break;
		case 3:
			c_str[i] = 'T';
			break;
		}
	}
	c_str[len] = '\0';
	return c_str;
}
