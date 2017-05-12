#ifndef CMINSKETCH_H_
#define CMINSKETCH_H_

#pragma once
#include <stdlib.h>
#include <math.h>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include "../seq/types.h"
#include "counts_table.h"

// multi-stream count-min sketch frequency table
template <int NBITS>
class cms_table_t : public counts_table_t {
public:
	typedef uint64 counter_block_t;
	static const int BITS_PER_COUNTER_BLOCK = 64;
	static const int COUNTERS_PER_BLOCK = BITS_PER_COUNTER_BLOCK / NBITS;
	static const uint64 MAX_COUNT_VAL = ((1ULL << NBITS) - 1);
	static const uint64  COUNTER_MASK = MAX_COUNT_VAL;
	
	// table configuration
	typedef uint64 table_dim_t;
	table_dim_t n_rows; // number of hash functions
	table_dim_t n_cols; // hash funtion range
	int n_streams;
	table_dim_t n_blocks;
	counter_block_t* counters; // flat compressed counter array
	table_dim_t n_counter_blocks_per_row;
	table_dim_t n_counter_blocks_per_stream;

	// internal coordinates in the compressed representation
	struct cms_dim2_t {
		table_dim_t block_id;
		table_dim_t counter_id;
		cms_dim2_t() : block_id(0), counter_id(0) {};
		cms_dim2_t(const table_dim_t _block_id, const table_dim_t _counter_id) : block_id(_block_id), counter_id(_counter_id) {};
		cms_dim2_t(const int stream_id, const table_dim_t row_id, const table_dim_t col_id, const table_dim_t n_counter_blocks_per_stream, const table_dim_t n_counter_blocks_per_row) {
			block_id = stream_id*n_counter_blocks_per_stream + row_id*(n_counter_blocks_per_row) + col_id/COUNTERS_PER_BLOCK;
			counter_id = col_id  % COUNTERS_PER_BLOCK;
		}
	} ;
	
	// logical coordinates
	struct cms_dim3_t {
		int stream_id;
		table_dim_t row_id;
		table_dim_t col_id;
		cms_dim3_t() : stream_id(0), row_id(0), col_id(0) {};
		cms_dim3_t(const int _stream_id, const table_dim_t _row_id, const table_dim_t _col_id) : stream_id(_stream_id), row_id(_row_id), col_id(_col_id) {};
	} ;

	// counter masking/manipulation
	counter_block_t counter_block_masks[COUNTERS_PER_BLOCK];

	// kmer hashing, multiply-shift hash functions (Dietzfelbinger et al)
	typedef uint64 hash_size_t;
	static const int BITS_IN_KMER_HASH = 64;
	std::vector<hash_size_t> hash_funcs;
	int n_col_log;

	// counts stats
	counter_t max_count;
	counter_t avg_count;
	uint64 n_kmers;
	uint64 n_clipped;

	// constructors
	cms_table_t(): cms_table_t(0,0,0) {}
	cms_table_t(const table_dim_t nrows, const table_dim_t ncols, const uint32 n_input_streams)  {
		n_rows = nrows;
		n_cols = ncols;
		n_streams = n_input_streams;
		n_counter_blocks_per_row = ceil((1.0 * n_cols) / COUNTERS_PER_BLOCK);
		n_counter_blocks_per_stream = n_rows * n_counter_blocks_per_row;
		n_blocks = n_streams * n_counter_blocks_per_stream;
		counters = new counter_block_t[n_blocks];
		init_masks();
		init_hash_funcs();
		std::cout << "Resulting CMS table size: n_rows = " << n_rows << " n_cols = " << n_cols << " n_streams = " << n_streams << " n_blocks = " << n_blocks << "\n";
	}
	
	static cms_table_t* size_params(const table_dim_t nrows, const table_dim_t ncols, const uint32 n_input_streams)  {
		return new cms_table_t(nrows, ncols, n_input_streams);
	}
	static cms_table_t* err_params(const float delta, const float epsilon, const uint32 n_input_streams)  {
		std::cout << "Initializing CMS table with epsilon = " << epsilon << " delta = " << delta << "\n";
		return new cms_table_t(ceil(log(1 / delta)), ceil(exp(1) / epsilon), n_input_streams);
	}
	void init_masks() {
		for(int i = 0; i < COUNTERS_PER_BLOCK; i++) {
			counter_block_masks[i]  = 0;
			counter_block_masks[i] |= COUNTER_MASK;
			counter_block_masks[i] <<= (i * BITS_PER_COUNTER);
		}
	}

	bool is_prime(table_dim_t n) {
		if(n <= 1) return false;
		if(n <= 3) return true;
		if(n % 2 == 0 || n % 3 == 0) return false;
		for(table_dim_t i = 5; i*i <= n; i += 6) {
			if(n % i == 0 || n % (i+2) == 0) return false;
		}
		return true;
	}

	hash_size_t get_closest_prime(table_dim_t n_buckets) {
		if(n_buckets % 2 == 0) n_buckets--;
		while(true) {
			if(is_prime(n_buckets)) return n_buckets;
			n_buckets -= 2;
		}
		return 0;
	}
	
	void init_hash_funcs() {
		table_dim_t bound = n_cols;
		for(table_dim_t i = 0; i < n_rows; i++) {
			table_dim_t prime = get_closest_prime(bound); //2 * rand() + 1; // odd multipliers
			if(prime == 0) break;
			hash_funcs.push_back(prime);
			bound = prime - 2;
		}
		if(hash_funcs.size() != n_rows) {
			std::cerr << "ERROR: Could not generate hash functions for the index \n";
			exit(1);
		}
		//n_col_log = log2(n_cols);
	}
	
	// compute the kmer hash value for a given row
	inline  table_dim_t hash(const kmer_2bit_t& kmer, const table_dim_t row_id) const {
		return kmer % hash_funcs[row_id]; // (hash_funcs[row_id]* kmer.packed_rep) % n_cols; //>> (BITS_IN_KMER_HASH - n_col_log);
	}

	// unpack the counter value of a given stream
	inline counter_t get_counter_val(const cms_dim2_t& c) const {
		return (counters[c.block_id] & counter_block_masks[c.counter_id]) >> (c.counter_id * BITS_PER_COUNTER);
	}

	// increment the counter value of a given stream by the given value
	// increment the count until max is reached (avoid overflow)
	void inc_counter_val(const cms_dim2_t& c, const counter_t i) {
		const counter_block_t current_count = get_counter_val(c);
		if(current_count < MAX_COUNT_VAL) { 
			counters[c.block_id] &= ~counter_block_masks[c.counter_id];
			const counter_block_t new_val = current_count + i;
			counters[c.block_id] |= ((new_val ) << (c.counter_id * BITS_PER_COUNTER));
		} else {
			n_clipped++;
		}
	}

	// ---- interface ---- //

	virtual ~cms_table_t() {
		delete counters;
	}

	virtual void clear() {}
	
	virtual int get_n_streams() { 
		return n_streams;
	}

	// insert a kmer into the sketch
	virtual void insert(const kmer_2bit_t& kmer, const int stream_id)  { 
		for(table_dim_t row = 0; row < n_rows; row++) {
			const cms_dim2_t flat_coords(stream_id, row, hash(kmer, row), n_counter_blocks_per_stream,  n_counter_blocks_per_row);
			inc_counter_val(flat_coords, 1);
		}
	}
    
	// lookup the kmer count in the sketch
	virtual counter_t lookup(const kmer_2bit_t& kmer, const int stream_id) const {
		counter_t min_count = MAX_COUNT_VAL;
		for(table_dim_t row = 0; row < n_rows; row++) {
			const cms_dim2_t flat_coords(stream_id, row, hash(kmer, row), n_counter_blocks_per_stream,  n_counter_blocks_per_row);
			const counter_t c = get_counter_val(flat_coords);
			if(c < min_count) {
				min_count = c;
			}
		}
		return min_count;
    	}

	// write the table to file
	virtual void save_to_file(const std::string& fname, const int n_count_bits) {
		std::ofstream file;
		file.open(fname.c_str(), std::ios::out | std::ios::binary | std::ios::app);
		int nbits = NBITS;
		file.write(reinterpret_cast<char*>(&nbits), sizeof(int));
		file.write(reinterpret_cast<char*>(&n_rows), sizeof(n_rows));
		file.write(reinterpret_cast<char*>(&n_cols), sizeof(n_cols));
		file.write(reinterpret_cast<char*>(&n_streams), sizeof(n_streams));
		file.write(reinterpret_cast<char*>(&n_blocks), sizeof(n_blocks));
		file.write(reinterpret_cast<const char*>(counters), n_blocks*sizeof(counter_block_t));
	}
	
	// load the table from file
	virtual long int load_from_file(const std::string& fname, long int file_offset) {
		std::ifstream file;
		file.open(fname.c_str(), std::ios::in | std::ios::binary);
		int nbits_used;
		file.read(reinterpret_cast<char*>(&nbits_used), sizeof(nbits_used));
		//if(nbits_used != NBITS) {
		//	std::cerr << "ERROR: Mismatch with the number of bits used during index construction: expected " << NBITS << ", built with " << nbits_used << ". [Recompile accordingly.]\n";
		//	exit(1);
		//}
		file.read(reinterpret_cast<char*>(&n_rows), sizeof(n_rows));
		file.read(reinterpret_cast<char*>(&n_cols), sizeof(n_cols));
		file.read(reinterpret_cast<char*>(&n_streams), sizeof(n_streams));
		file.read(reinterpret_cast<char*>(&n_blocks), sizeof(n_blocks));
		n_counter_blocks_per_row = ceil((1.0 * n_cols) / COUNTERS_PER_BLOCK);
		n_counter_blocks_per_stream = n_rows * n_counter_blocks_per_row;
		counters = new counter_block_t[n_blocks];
		file.read(reinterpret_cast<char*>(counters), n_blocks*sizeof(counter_block_t));
		init_hash_funcs();
		init_masks();
		std::cout << "Loaded CMS table size: n_rows = " << n_rows << " n_cols = " << n_cols << " n_streams = " << n_streams << " n_blocks = " << n_blocks << "\n";
		return 0;
	}
	
	virtual void print_stats() {
		int max_count = 0;
		uint64 ncells = n_streams*n_rows*n_cols;
		uint64 nz = 0;
		uint64 sum = 0;
		for(int stream = 0; stream < n_streams; stream++) {
			for(table_dim_t row = 0; row < n_rows; row++) {
				for(table_dim_t col = 0; col < n_cols; col++) {
					const counter_t count = get_counter_val(cms_dim2_t(stream, row, col, n_counter_blocks_per_stream,  n_counter_blocks_per_row));
					sum += count;
					if(count > max_count) max_count = count;
					if(count == 0) nz++;
				}
			}
		}
		
		float avg = ((float) sum)/(ncells-nz);
		std::cout << "Index stats:" << "\n n_entries: " << ncells  << "\n n_zeroes: " << nz  << "\nmax value: " << max_count  << "\navg value: " << avg << "\n size(GB): " << ((float)n_blocks)*sizeof(counter_block_t)/1024/1024/1024 << "\n";
	}
};

template <int NBITS>
class discretized_cms_table_t : public cms_table_t<NBITS> {
public:
	typedef typename cms_table_t<NBITS>::table_dim_t table_dim_t;
	counter_t bucket_range;
	counter_t  min_range;
	counter_t  max_range;
	counter_t max_counter_value;
	
	//discretized_cms_table_t() {}

	// construct a discretized CMS table from a pre-computed initial CMS counts table
	// values under min_range are treated as zero (first bucket)
	// values over max_range are capped at max_range (last bucket)
	template <int BITS_PER_COUNTER_INIT>
	discretized_cms_table_t(const cms_table_t<BITS_PER_COUNTER_INIT>* init_table, const counter_t  _min_range, const counter_t  _max_range, const counter_t _max_val)
	: cms_table_t<NBITS>(init_table->n_rows, init_table->n_cols, init_table->n_streams) {
		min_range = _min_range;
		max_range = _max_range;
		max_counter_value = _max_val;
		bucket_range = (max_range-min_range)/(this->MAX_COUNT_VAL-2);
		
		for(int stream = 0; stream < this->n_streams; stream++) {
			for(table_dim_t row = 0; row < this->n_rows; row++) {
				for(table_dim_t col = 0; col < this->n_cols; col++) {
					typename cms_table_t<BITS_PER_COUNTER_INIT>::cms_dim2_t coords(stream, row, col, init_table->n_counter_blocks_per_stream,  init_table->n_counter_blocks_per_row);
					const counter_t count = init_table->get_counter_val(coords);
					const counter_t bucket_id = discretize(count);
					cms_table_t<NBITS>::inc_counter_val(typename discretized_cms_table_t<NBITS>::cms_dim2_t(stream, row, col, this->n_counter_blocks_per_stream,  this->n_counter_blocks_per_row), bucket_id);
				}
			}
		}
	}
	 
	// returns the bucket id into which the counter value falls
	inline counter_t discretize(const counter_t counter) const {
		if(counter < min_range) return 0;	
		if(counter >= max_range) return this->MAX_COUNT_VAL;	
		return counter/bucket_range + 1;
    }
	
	// returns the representitive counter value for the given bucket
	// mid value for internal buckets
	inline counter_t revert(const counter_t bucket_id) const {
		if(bucket_id == 0) return 0;	
		if(bucket_id == this->MAX_COUNT_VAL) return max_counter_value;	
		return min_range + (bucket_id-1)*bucket_range+bucket_range/2;
	}
	
	virtual void insert(const kmer_t& kmer, const int stream_id) { /* todo: unsupported */ }
	
	virtual ~discretized_cms_table_t() {}

	 // return the mid value in the bucket
	 virtual counter_t lookup(const kmer_t& kmer, const int stream_id) const {
		 counter_t bucket_id = cms_table_t<NBITS>::lookup(kmer, stream_id); 
		 counter_t count =  revert(bucket_id);
		 return count;
	 }
	 
	 virtual void save_to_file(std::ofstream& file) {
		cms_table_t<NBITS>::save_to_file(file);
		file.write(reinterpret_cast<char*>(&min_range), sizeof(min_range));
		file.write(reinterpret_cast<char*>(&max_range), sizeof(max_range));
		file.write(reinterpret_cast<char*>(&max_counter_value), sizeof(max_counter_value));
		file.write(reinterpret_cast<char*>(&bucket_range), sizeof(bucket_range));
	 }
	 
	 virtual void load_from_file(std::ifstream& file) {
		cms_table_t<NBITS>::load_from_file(file);
		file.read(reinterpret_cast<char*>(&min_range), sizeof(min_range));
		file.read(reinterpret_cast<char*>(&max_range), sizeof(max_range));
		file.read(reinterpret_cast<char*>(&max_counter_value), sizeof(max_counter_value));
		file.read(reinterpret_cast<char*>(&bucket_range), sizeof(bucket_range));
		
		std::cout << "range " << (int) bucket_range << " min " << (int) min_range << " max " << (int) max_range << "\n";
	 }
};
	
#endif
