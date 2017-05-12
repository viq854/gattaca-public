#ifndef MHSKETCH_H_
#define MHSKETCH_H_

#pragma once

#include <emmintrin.h>
#include <smmintrin.h>
#include "../seq/types.h"
#include "counts_table.h"
#include "basic_map.h"
#include "mph.h"
#include <stdlib.h>
#include <math.h>
#include <unordered_map>

#if defined(_OPENMP)
#include <omp.h>
#endif

static const int BITS_IN_KMER_HASH = 64;

typedef uint64 hash_size_t;
struct minhash_fp_t {
	std::vector<hash_size_t> v;
};


struct rand_range_generator_t {
	int rand_in_range(int n) {
		int r, rand_max = RAND_MAX - (RAND_MAX % n);
		while ((r = rand()) >= rand_max);
		return r / (rand_max / n);
	}
};

// multi-stream minhash index
class minhash_table_t : public counts_table_t  {
public:
	int k;
	int n_streams;
	int FP_len;
	int FP_proj_len;
	int n_tables;
	
	// kmer hashing, multiply-shift hash functions (Dietzfelbinger et al)
	std::vector<hash_size_t> fp_hash_funcs;
	std::vector<std::vector<hash_size_t> > fp_proj_funcs;
	std::vector<std::vector<int> > fp_proj_ind;
	typedef uint64 proj_hash_t;
	std::vector<basic_table_t> tables;
	std::vector<mphf_table_t> static_tables;

	minhash_table_t(): minhash_table_t(0,0,0,0,0) {}
	
	minhash_table_t(const int kmer_len, const int fp_len, const int n_proj, const int proj_len, const int n_input_streams) {
		k = kmer_len;
		FP_len = fp_len;
		n_tables = n_proj;
		FP_proj_len = proj_len;
		n_streams  = n_input_streams;
		init_hash_funcs();
		tables.resize(n_tables);
		static_tables.resize(n_tables);
	}
	
	void init_hash_funcs() {
		fp_hash_funcs.resize(FP_len);
		for(int i = 0; i < FP_len; i++) {
			hash_size_t a = 2ULL*rand() + 1ULL; // odd multipliers
			fp_hash_funcs[i] = a;
		} 
	
		fp_proj_funcs.resize(n_tables);
		fp_proj_ind.resize(n_tables);
		rand_range_generator_t rgen;
		std::vector<int> idx(FP_len); 
		for(int i = 0; i < n_tables; i++) {
			for(int k = 0; k < FP_len; k++) {
				idx[k] = k;
			}
			// pick random indices from the sketch
			fp_proj_funcs[i].resize(FP_proj_len);
			fp_proj_ind[i].resize(FP_proj_len);
			int cnt = 0;
			int len = FP_len;
			while(cnt < FP_proj_len) {
				int j = rgen.rand_in_range(len); // exclude 0
				fp_proj_ind[i][cnt] = idx[j];
				fp_proj_funcs[i][cnt] = 2ULL * rand() + 1ULL;
				idx[j] = idx[len-1];
				cnt++;
				len--;
			}
		}
	}
	
	bool compute_fp(const std::string seq, minhash_fp_t& fp) {
		const int n_kmers = seq.size() - k + 1;
		kmer_2bit_t v[n_kmers]  __attribute__((aligned(16)));;
		int n_valid_kmers = 0;
		kmer_parser_t seq_parser;
		seq_parser.init(seq, k);
		kmer_t kmer;
		while(seq_parser.get_next_kmer(kmer)) {
			if(!kmer.valid) continue;
			v[n_valid_kmers] = kmer.packed_rep; 
			n_valid_kmers++;
		}
		if(n_valid_kmers <= k) {
				return false;
		}
		fp.v.resize(FP_len);
		for(int h = 0; h < FP_len; h++) {
			const hash_size_t s = fp_hash_funcs[h];
			hash_size_t curr_min = s*v[0];	
			for(int i = 1; i < n_valid_kmers; i++) {
				hash_size_t p = v[i]*s;
				if(p < curr_min) {
					curr_min = p;
				}
			}
			fp.v[h] = curr_min;
		}
		return true;
	}
	
	virtual proj_hash_t compute_fp_proj(const minhash_fp_t& fp, const int proj_id) const {
		proj_hash_t s = 0;
		for(int i = 0; i < FP_proj_len; i++) {
			s += fp_proj_funcs[proj_id][i]*fp.v[fp_proj_ind[proj_id][i]];
		}
		return s;
	}
	
	// ---- interface ---- //
	virtual void clear() {}
	
	virtual ~minhash_table_t() {}

	virtual int get_n_streams() { 
		return n_streams;
	}

	virtual void insert(const std::string& seq, const int stream_id)  { 
		minhash_fp_t fp;
		if(!compute_fp(seq, fp)) return;
	
		// compute the projections
		#pragma omp parallel for
		for(int t = 0; t < n_tables; t++) { // for each hash table
	    		proj_hash_t key = compute_fp_proj(fp, t);
			tables[t].insert(key, stream_id);
			//std::cout << key << " ";
		} //std::cout << "\n";
	}
    
	// lookup the kmer count in the sketch
	virtual counter_t lookup(const minhash_fp_t& fp, const int stream_id) const {
		std::vector<counter_t> proj_counts(n_tables);
		for(int t = 0; t < n_tables; t++) { // for each hash table
	    		proj_hash_t key = compute_fp_proj(fp, t);
			proj_counts[t] = static_tables[t].lookup(key, stream_id);		
		}
		// return median
		std::sort(proj_counts.begin(), proj_counts.end());
		return proj_counts[proj_counts.size()/2];
	}

	virtual void done() {
		#pragma omp parallel for
		for(int t = 0; t < n_tables; t++) { // for each hash table
			std::vector<kmer_2bit_t> keys;
			std::vector<counter_t> key_counts;
			tables[t].get_key_values(keys, key_counts);
			tables[t].clear();
			static_tables[t].init(keys, key_counts);
		}
		std::vector<basic_table_t>().swap(tables);
	}
	
	// write the table to file
	virtual void save_to_file(const std::string& fname, const int n_count_bits) {
		std::ofstream file;
		file.open(fname.c_str(), std::ios::out | std::ios::binary | std::ios::app);
		file.write(reinterpret_cast<char*>(&n_streams), sizeof(n_streams));
		file.write(reinterpret_cast<char*>(&k), sizeof(k));
		file.write(reinterpret_cast<char*>(&n_tables), sizeof(n_tables));
		file.write(reinterpret_cast<char*>(&FP_len), sizeof(FP_len));
		file.write(reinterpret_cast<char*>(&FP_proj_len), sizeof(FP_proj_len));
		file.close();	
	
		for(int t = 0; t < n_tables; t++) { // for each hash table
			static_tables[t].save_to_file(fname, n_count_bits);
		}
	}
	
	// load the table from file
	virtual long int load_from_file(const std::string& fname, long int file_offset) {
		std::ifstream file;
		file.open(fname.c_str(), std::ios::in | std::ios::binary);
		file.seekg(file_offset, file.beg);
		file.read(reinterpret_cast<char*>(&n_streams), sizeof(n_streams));
		file.read(reinterpret_cast<char*>(&k), sizeof(k));
		file.read(reinterpret_cast<char*>(&n_tables), sizeof(n_tables));
		file.read(reinterpret_cast<char*>(&FP_len), sizeof(FP_len));
		file.read(reinterpret_cast<char*>(&FP_proj_len), sizeof(FP_proj_len));
		long int s = file.tellg();
		file.close();
		
		init_hash_funcs();
		tables.resize(n_tables);
		static_tables.resize(n_tables);
		std::cout << "Minhash config: " << FP_len << " " << n_tables << "\n"; 
		for(int t = 0; t < n_tables; t++) { // for each hash table
			s = static_tables[t].load_from_file(fname, s); 
		}
		//init_hash_funcs();
		std::cout << "Loaded index \n";
		return s;
	}
	
	virtual void print_stats() {
		for(int t = 0; t < n_tables; t++) { 
			//tables[t].print_stats();
		}
	}

	virtual counter_t lookup(const kmer_2bit_t& kmer, const int stream_id) const {
		return 0;
	}
	
	virtual void insert(const kmer_2bit_t& kmer, const int stream_id) {}

};

#endif
