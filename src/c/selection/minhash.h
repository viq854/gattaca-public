#ifndef MH_H_
#define MH_H_

#pragma once

#include <emmintrin.h>
#include <smmintrin.h>
#include "../seq/types.h"
#include <stdlib.h>
#include <math.h>
#include <unordered_map>

#if defined(_OPENMP)
#include <omp.h>
#endif

typedef uint64 hash_size_t;
struct minhash_fp_t {
	std::vector<hash_size_t> v;
};

struct minhash_rand_range_generator_t {
	int rand_in_range(int n) {
		int r, rand_max = RAND_MAX - (RAND_MAX % n);
		while ((r = rand()) >= rand_max);
		return r / (rand_max / n);
	}
};

typedef uint64 proj_hash_t;
struct bucket_entry {
	int sample_id;
	proj_hash_t hash;
};

struct sample_support_t {
	int id;
	int support;
	sample_support_t(int _id, int _support): id(_id), support(_support) {};
};

bool sort_by_support(const sample_support_t& lhs, const sample_support_t& rhs) {
   return lhs.support > rhs.support;
}

class minhash_t {
public:
	int k;
	int FP_len; // sample fingerprint length
	int B;
	
	// kmer hashing, multiply-shift hash functions (Dietzfelbinger et al)
	std::vector<hash_size_t> fp_hash_funcs;
	std::vector<std::vector<std::vector<bucket_entry>>> tables; // buckets of sample ids
	
	minhash_t(): minhash_t(0, 0) {}
	
	minhash_t(const int kmer_len, const int fp_len) {
		B = 131072;
		k = kmer_len;
		FP_len = fp_len;
		init_hash_funcs();
	}
	
	void init_hash_funcs() {
		fp_hash_funcs.resize(FP_len);
		for(int i = 0; i < FP_len; i++) {
			hash_size_t a = 2ULL*rand() + 1ULL; // odd multipliers
			fp_hash_funcs[i] = a;
		}
	}
	
	void init_index_tables() {
		tables.resize(FP_len);
		for(int i = 0; i < FP_len; i++) {
			tables[i].resize(B);
		}
	}

	bool compute_fp(const std::vector<kmer_2bit_t>& keys, minhash_fp_t& fp) {
		if(keys.size() == 0) {
			std::cout << "Empty key set!\n";
			return false;
		}
		fp.v.resize(FP_len);
		#pragma omp parallel for
		for(int h = 0; h < FP_len; h++) {
			const hash_size_t s = fp_hash_funcs[h];
			fp.v[h] = s*keys[0];
			for(uint64 i = 1; i < keys.size(); i++) {
				const hash_size_t p = keys[i]*s;
				if(p < fp.v[h]) {
					fp.v[h] = p;
				}
			}
		}
		return true;
	}
	
	// insert the sample into the appropriate buckets based on its fingerprint projections
	void insert(const minhash_fp_t& fp, const int sample_id)  { 
		#pragma omp parallel for
		for(int t = 0; t < FP_len; t++) { // for each hash table
			const proj_hash_t key = fp.v[t];
			const int bucket_id = key % B; 
			bucket_entry e;
			e.sample_id = sample_id;
			e.hash = key;
			tables[t][bucket_id].push_back(e);
		}
	}
	
	// lookup the sample ids in all the buckets given by the fingerprint
	void lookup(const minhash_fp_t& query_fp, std::vector<sample_support_t>& id_counts) {
		std::vector<int> sample_ids;
		for(int t = 0; t < FP_len; t++) { // for each hash table
			const proj_hash_t key = query_fp.v[t];
			const int bucket_id = key % B;
			for(unsigned int i = 0; i < tables[t][bucket_id].size(); i++) {
				bucket_entry e = tables[t][bucket_id][i];
				if(e.hash != key) continue;
				sample_ids.push_back(e.sample_id);
			}			
		}
		std::sort(sample_ids.begin(), sample_ids.end());
		if(sample_ids.size() == 0) return;
	
		int c = 1;
		for(uint64 i = 1; i < sample_ids.size(); i++) {
			if(sample_ids[i] == sample_ids[i-1]) {
				c++;
			} else {
				id_counts.push_back(sample_support_t(sample_ids[i-1], c));
				c = 1;
			}
		}
		id_counts.push_back(sample_support_t(sample_ids[sample_ids.size()-1], c));
		std::sort(id_counts.begin(), id_counts.end(), sort_by_support);
		//std::vector<int>::iterator it = std::unique(sample_ids.begin(), sample_ids.end());
		//sample_ids.erase(it, sample_ids.end());
	}

	static void save_fp_to_file(std::string& fname, const int k, const minhash_fp_t& fp) {
		fname += std::string(".k");
		fname += std::to_string(k);
		fname += std::string(".L");
		fname += std::to_string(fp.v.size());
		fname += std::string(".gaf");
		std::ofstream file;
		file.open(fname.c_str(), std::ios::out | std::ios::binary);
		file.write(reinterpret_cast<const char*>(&(fp.v[0])), fp.v.size()*sizeof(hash_size_t));
		file.close();
	}

	static void load_fp_from_file(const std::string& fname, const int fp_len, minhash_fp_t& fp) {
		fp.v.resize(fp_len);
		std::ifstream file;
		file.open(fname.c_str(), std::ios::in | std::ios::binary);
		file.read(reinterpret_cast<char*>(&(fp.v[0])), fp_len*sizeof(hash_size_t));
		file.close();
	}

	static int compare_fps(const minhash_fp_t& fp1, const minhash_fp_t& fp2) {
		int count = 0;
		for(unsigned int i = 0; i < fp1.v.size(); i++) {
			if(fp1.v[i] == fp2.v[i]) {
				count++;
			}
		}
		return count;
	}

	void save_index_to_file(std::string& fname) {
		fname += std::string(".k");
		fname += std::to_string(k);
		fname += std::string(".L");
		fname += std::to_string(FP_len);
		fname += std::string(".gaf.idx");
		std::ofstream file;
		file.open(fname.c_str(), std::ios::out | std::ios::binary);
		for(int i = 0; i < FP_len; i++) {
			for(int j = 0; j < B; j++) {
				int table_size = tables[i][j].size();
				file.write(reinterpret_cast<char*>(&table_size), sizeof(table_size));
				file.write(reinterpret_cast<const char*>(&(tables[i][j][0])), table_size*sizeof(bucket_entry));
			}
		}
		file.close();
	}

	void load_index_from_file(const std::string& fname) {
		std::ifstream file;
		file.open(fname.c_str(), std::ios::in | std::ios::binary);
		for(int i = 0; i < FP_len; i++) {
			for(int j = 0; j < B; j++) {
				int table_size;
				file.read(reinterpret_cast<char*>(&table_size), sizeof(table_size));
				tables[i][j].resize(table_size);
				file.read(reinterpret_cast<char*>(&(tables[i][j][0])), table_size*sizeof(bucket_entry));
			}
		}
	}
};

#endif



