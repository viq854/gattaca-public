#ifndef PERFECTHASH_H_
#define PERFECTHASH_H_

#include <iostream>
#include <cstdint>
#include <cstring>
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <ctype.h>
#include <stdlib.h>
#include <bf.h>
#if defined(_OPENMP)
#include <omp.h>
#endif
extern "C" {
#include <cmph.h>
}
#include "../seq/types.h"
#include "counts_table.h"

// exact counts table using a minimal perfect hash function
class mphf_table_t : public counts_table_t {
public:
	static const int MIN_KMER_COUNT = 1;
	typedef struct {
		char val[8];
	} key_t; // converted struct key to match the cmph API
	
	cmph_t* hash; // minimal perfect hash function
	uint64 n_keys; // number of distinct keys
	counter_t* counts; // array of counts
	bf::basic_bloom_filter* bf; // bloom filter of distinct keys used during mphf construction 

	mphf_table_t() {} // todo: clean up empty constructor, refactor
	mphf_table_t(std::vector<kmer_2bit_t>& keys) {
		std::cout << "Sorting and counting... " << keys.size() << " kmers\n";
		std::sort(keys.begin(), keys.end());

		std::vector<kmer_2bit_t> keys_distinct;
		std::vector<counter_t> key_counts;
		keys_distinct.reserve(keys.size());
		key_counts.reserve(keys.size());
		counter_t c = 1;
		for(uint64 i = 1; i < keys.size(); i++) {
			if(keys[i] == keys[i-1]) {
				if(c < std::numeric_limits<counter_t>::max()) {
					c++;
				}
			} else {
				if(c > MIN_KMER_COUNT) {
					keys_distinct.push_back(keys[i-1]);
					key_counts.push_back(c);
				}
				c = 1;
			}
		}
		if(c > MIN_KMER_COUNT) {
			keys_distinct.push_back(keys[keys.size()-1]);
			key_counts.push_back(c);
		}
		n_keys = keys_distinct.size();
		//std::cout << "Number of distinct kmers: " << n_keys << "\n";
		keys.clear();
		init(keys_distinct, key_counts);
	}

	void init(std::vector<kmer_2bit_t>& keys, const std::vector<counter_t>& key_counts) {
		// construct the mphf
		std::cout << "Building the mphf...\n";
		key_t* key_structs = new key_t[n_keys];
		#pragma omp parallel for        
		for(uint64 i = 0; i < n_keys; i++) {
				get_key(keys[i], key_structs[i]);
		}
		build_mpfh(key_structs);
		free(key_structs);

		// store the keys in a bloom filter
		std::cout << "Building the bloom filter...\n";
		bf = new bf::basic_bloom_filter(0.05, n_keys, 0, false, false);
		for(uint64 i = 0; i < n_keys; i++) {
				bf->add(keys[i]);
		}

		// count all input keys using the mphf
		std::cout << "Scattering counts...\n";
		counts = new counter_t[n_keys];
		#pragma omp parallel for
		for(uint64 i = 0; i < n_keys; i++) {
				const unsigned int id = get_id(keys[i]);
				counts[id] = key_counts[i];
		}
		keys.clear();
		std::cout << "Index construction done!\n";
	} 
	
	static inline void get_key(const kmer_2bit_t& key_in, key_t& key_out) {
		std::memcpy(key_out.val, &key_in, sizeof(key_out.val));
	}
	
	inline unsigned int get_id(const kmer_2bit_t& key) const {
		key_t key_struct;
		get_key(key, key_struct);
		return cmph_search(hash, key_struct.val, sizeof(key_struct.val));
	}	
	
	void build_mpfh(key_t* keys) {
		cmph_io_adapter_t* source = cmph_io_struct_vector_adapter(keys, (cmph_uint32) sizeof(key_t), 0, sizeof(key_t), n_keys);
		cmph_config_t* config = cmph_config_new(source);
		cmph_config_set_algo(config, CMPH_BDZ);
		hash = cmph_new(config);
		if(hash == NULL) {
			std::cout << "ERROR: null hash \n";
			exit(-1);
		}
		cmph_config_destroy(config);
	}
	
	virtual ~mphf_table_t() {
		cmph_destroy(hash);
		bf->clear();
		free(bf);
		free(counts);
	}
	virtual void clear() {}
	
	// insert a kmer into the sketch (the kmer must be one of the distinct keys)
	virtual void insert(const kmer_2bit_t& key, const int stream_id)  { 
		const unsigned int id = get_id(key);
		counter_t newval, curr;
		do {
			curr = counts[id];
			if(curr == std::numeric_limits<counter_t>::max()) {
				return;
			}
			newval = curr + 1;
		} while (!__sync_bool_compare_and_swap(&counts[id], curr, newval));
	}

	// lookup the kmer count in the sketch
	virtual counter_t lookup(const kmer_2bit_t& key, const int stream_id) const {
		// check the bloom filter to ensure the key is part of the sketch
		if(!bf->lookup(key)) return 0;
		return counts[get_id(key)];
	}

	// write the table to file
	virtual void save_to_file(const std::string& fname, const int n_count_bits) {
		save_mphf_aux(fname);
		std::ofstream file;
		file.open(fname.c_str(), std::ios::out | std::ios::binary | std::ios::app);
		file.write(reinterpret_cast<char*>(&n_keys), sizeof(n_keys));

		// bloom // TEMP getting around the lack of io in the bf lib
		const bf::bitvector& bits = bf->storage();
		long int bits_size = bits.bits_.size();
		file.write(reinterpret_cast<const char*>(&bits.num_bits_), sizeof(bits.num_bits_));
		file.write(reinterpret_cast<char*>(&bits_size), sizeof(bits_size));		
		file.write(reinterpret_cast<const char*>(&bits.bits_[0]), bits_size*sizeof(bits.bits_[0]));
		
		if(n_count_bits == std::numeric_limits<counter_t>::digits) {
			file.write(reinterpret_cast<char*>(&counts[0]), n_keys*sizeof(counter_t));
		} else {
			int n_clipped = 0;
			uint8 maxv = 255; // reduce to char temp
			for(uint64 i = 0; i < n_keys; i++) {
				if(counts[i] >= maxv) {
					file.write(reinterpret_cast<char*>(&maxv), sizeof(uint8));
					n_clipped++;
				}
				else {
					file.write(reinterpret_cast<char*>(&counts[i]), sizeof(uint8));
				}
			}
			std::cout << "Number of clipped counts: " << n_clipped << "\n";
		}
		file.close();
	}
	
	void save_mphf_aux(const std::string& fname) {
		FILE* mphf_fd = fopen(fname.c_str(), "a");
		cmph_dump(hash, mphf_fd); 
		fclose(mphf_fd);
	}
	
	// load the table from file
	virtual long int load_from_file(const std::string& fname, long int file_offset) {
		long int file_pos = load_mpfh_aux(fname, file_offset);
		std::ifstream file;
		file.open(fname.c_str(), std::ios::in | std::ios::binary);
		file.seekg(file_pos, file.beg);
		file.read(reinterpret_cast<char*>(&n_keys), sizeof(n_keys));

		// bloom
		bf = new bf::basic_bloom_filter(0.05, n_keys, 0, false, false);
		size_t nbits;
		file.read(reinterpret_cast<char*>(&nbits), sizeof(nbits));
		long int bits_size;
        file.read(reinterpret_cast<char*>(&bits_size), sizeof(bits_size));
		file.read(reinterpret_cast<char*>(&bf->bits_.bits_[0]), bits_size*sizeof(bf->bits_.bits_[0]));

		counts = new counter_t[n_keys];
		file.read(reinterpret_cast<char*>(&counts[0]), n_keys*sizeof(counter_t));
		long int s = file.tellg();
		file.close();
		return s;
	}
	
	long int load_mpfh_aux(const std::string& fname, long int file_offset) {
		FILE* mphf_fd = fopen(fname.c_str(), "r");
		fseek(mphf_fd, file_offset, SEEK_SET);
		hash = cmph_load(mphf_fd);
		long int s = ftell(mphf_fd);
		fclose(mphf_fd);
		return s;
	}
	
	// create the mphf for a given set of keys and save to disk
	// + bloom filter
	// note: all input keys are expected to be distinct
	void build_and_save_mphf_aux(const std::string& distinct_keys_fname) {
		std::ifstream keys_file;
		keys_file.open(distinct_keys_fname.c_str(), std::ios::in | std::ios::binary);
		if (!keys_file.is_open()) {
			std::cerr << "ERROR: Could not open keys file: " << distinct_keys_fname << "\n";
			exit(1);
		}
		std::cout << "Loading..." << distinct_keys_fname << "\n";
		keys_file.seekg(0, std::ios::end);
		long long int size = keys_file.tellg();
		keys_file.seekg(0, std::ios::beg);
		n_keys = size/sizeof(kmer_2bit_t);
		key_t* keys = new key_t[n_keys];
		std::cout << "Expected number of keys: " << n_keys << "\n";
		keys_file.read(reinterpret_cast<char*>(&keys[0]), n_keys*sizeof(key_t));
		std::cout << "Finished loading the keys \n";			

		// function
		std::string mphf_fname = std::string(distinct_keys_fname);
		mphf_fname += std::string(".mph");
		build_mpfh(keys);
		save_mphf_aux(mphf_fname);
		free(keys);
	}
	
	virtual void print_stats() {
		std::cout << "Number of distinct keys: " << n_keys << "\n";
	}
	
	virtual int get_n_streams() {
		return 1;
	}
};

#endif
