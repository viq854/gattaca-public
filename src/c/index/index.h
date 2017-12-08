#ifndef INDEX_H_
#define INDEX_H_

#pragma once

#include "../seq/types.h"
#include "../seq/io.h"
#include "counts_table.h"
#include "count_min_sketch.h"
#include "minhash_sketch.h"
#include "mph.h"
#include "basic_map.h"

typedef enum {CMS = 0, EXACT_BASIC = 1, MINHASH = 2,  EXACT_MPHF = 3} index_alg;

struct cms_params_t {
	int n_rows;
	int n_cols;
	float epsilon;
	float delta;
	
	void set_defaults() {
		epsilon = 0;
		delta = 0;
		n_rows = 8;
		n_cols = (1 << 17);
	}
};

// program parameters
struct index_params_t {    
	index_alg alg;		// type of table used for counting
    int k;				// length of the sequence kmers
    int min_kmer_freq;	// kmer filters
    int n_threads;		// multi-threading

	counter_t count_cutoff_max;
	counter_t count_cutoff_min;
	cms_params_t cms_config;

    void set_default_index_params() {
		k = 31;
		min_kmer_freq = 1;
		n_threads = 1;
		cms_config.set_defaults();
		count_cutoff_max = std::numeric_limits<counter_t>::max();
		count_cutoff_min = 0;
		alg = EXACT_MPHF;;
    }
};

// multi-sample counts index
struct index_t {
	// construction
	static void build_and_save(const std::vector<std::vector<std::string>>& files_to_index, index_params_t& params);
	static void build_and_save(const std::vector<std::vector<std::string>>& files_to_index, const int n_count_bits, index_params_t& params);
	static void discretize(const std::string& in_index_fname, const std::string& out_index_fname, const int nbits, const int min_r, const int max_r, const int max_v, index_params_t& params);
	void load(const std::vector<std::string>& index_files, index_params_t& params);

	// lookup per file
	void get_aggregate_counts(const std::string& fasta_fname, const std::string& output_fname, const bool is_median, const index_params_t& params);
	void get_counts_file(const std::string& fasta_fname, const int sample_id,  const std::string& output_fname,  const index_params_t& params);

	// lookup per sequence
	counter_t get_median_count(const std::string& seq, const int sample_id, const index_params_t& params) const;
	counter_t get_avg_count(const std::string& seq, const int sample_id, const index_params_t& params) const;
	void get_counts_seq(const std::string& seq, const int sample_id, const index_params_t& params) const;

	void print_basic_stats() {
		if(!count_tables) return;
		for(int i = 0; i < n_samples; i++) {
			count_tables[i]->print_stats();
		}
	}

	~index_t() {
		if(!count_tables) return;
		for(int i = 0; i < n_samples; i++) {
			delete count_tables[i];
		}
		delete count_tables;
	}

protected:
	static counts_table_t* init_counts_table(index_params_t& params);
	static counts_table_t* count(const std::string& seq_file, const index_params_t& params);
	static counts_table_t* count_mphf(const std::vector<std::string>& seq_file, const index_params_t& params);
	static counts_table_t* count_mphf(const std::string& seq_file, const index_params_t& params);
	static counts_table_t* count_minhash(const std::string& seq_file, const index_params_t& params);

private:
	counts_table_t** count_tables;
	int n_samples; 	// number of samples in ref panel
};

#endif
