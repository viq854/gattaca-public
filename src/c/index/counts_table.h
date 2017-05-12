#ifndef COUNTS_H_
#define COUNTS_H_

#include "../seq/types.h"
#include "../seq/io.h"

#if BITS_PER_COUNTER <= 8
	typedef uint8 counter_t;
#elif BITS_PER_COUNTER <= 16
	typedef uint16 counter_t;
#elif BITS_PER_COUNTER <= 32
	typedef uint32 counter_t;
#else
	typedef uint64 counter_t;
#endif

//#define MAX_COUNT_VAL ((1ULL << BITS_PER_COUNTER) - 1)
//#define COUNTER_MASK MAX_COUNT_VAL

//typedef uint64 counter_block_t;
//#define BITS_PER_COUNTER_BLOCK 64
//#define COUNTERS_PER_BLOCK (BITS_PER_COUNTER_BLOCK / BITS_PER_COUNTER)

// discretization
//typedef uint8 bucket_t;
//#define BITS_PER_BUCKET 4
//#define BUCKETS_PER_BLOCK (BITS_PER_COUNTER_BLOCK / BITS_PER_COUNTER)

// counts table base class (supports multi-stream tables)
class counts_table_t {
public:
	virtual counter_t lookup(const kmer_2bit_t& kmer, const int stream_id) const = 0;
	virtual void insert(const kmer_2bit_t& kmer, const int stream_id) = 0;
	virtual void insert(const std::string& seq, const int stream_id) {}
	virtual void save_to_file(const std::string& file, const int n_count_bits) = 0;
	virtual long int load_from_file(const std::string& file,  long int file_offset) = 0;
	virtual int get_n_streams() = 0;
	virtual void print_stats() = 0;
	virtual void clear() = 0;
	virtual ~counts_table_t() {};
};

#endif
