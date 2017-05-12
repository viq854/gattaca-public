#ifndef IO_H_
#define IO_H_

#include <vector>
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/basic.h>
#include "types.h"

#define BASE_IGNORE		4
static const unsigned char dna5_char[5] =  {'A', 'G', 'C', 'T', 'N'};
static const unsigned char dna5_comp[5] = {3/*A*/, 2/*G*/, 1/*C*/, 0/*T*/, 4/*N*/};
static const unsigned char dna5_table[256] = { // encoding: A=0, G=1, C=2, T=3, N=4
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4,0/*A*/,4,2/*C*/,4,4,4,1/*G*/,4,4,4,4,4,4,4/*N*/,4,
	4, 4, 4, 4,  3/*T*/,4,4,4,4, 4, 4, 4,  4, 4, 4, 4,
	4,0/*a*/,4,2/*c*/,4,4,4,1/*g*/,4,4,4,4,4,4,4/*n*/,4,
	4, 4, 4, 4,  3/*t*/,4,4,4,4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

struct contig_t {
	std::string name; 			// seq name
	std::string seq;				// seq sequence
	
	void print_contig() const {
		std::cout << name.c_str() << "\n";
		std::cout << seq.c_str() << "\n";
	}
};

struct read_t {
	std::string name; 			// read name
	std::string seq;				// read sequence
	std::string qual;				// quality scores

	void print_read() const {
		std::cout << name.c_str() << "\n";
		std::cout << seq.c_str() << "\n";
	}
};

typedef enum {FASTA, FASTQ} seq_file_reader_mode_t;
struct seq_file_reader_t {
	seqan::SeqFileIn file_handle;
	seq_file_reader_mode_t mode;
	int n_records;

	void open_file(const std::string& fname, seq_file_reader_mode_t file_mode = FASTA) {
		if (!seqan::open(file_handle, seqan::toCString(fname))) {
			std::cerr << "ERROR: Could not open FASTA/FASTQ file: " << fname << "\n";
			exit(1);
		}
		n_records = 0;
		mode = file_mode;
	}
	
	// load FASTA/FASTQ read records
	bool load_next_read(read_t& r) {
		//if(n_records != 0 && n_records % 1000000 == 0) {
			//std::cout << "Processed " << n_records << " reads\n";
		//}
		if(!seqan::atEnd(file_handle)) {
			if(mode == FASTQ) {
				seqan::readRecord(r.name, r.seq, r.qual, file_handle);
			}  else {
				seqan::readRecord(r.name, r.seq, file_handle);
			}
			n_records++;
			return true;
		} else {
			return false;
		}
	}
	
	// load FASTA contig records
	bool load_next_contig(contig_t & contig) {
		if(!seqan::atEnd(file_handle)) {
			seqan::readRecord(contig.name, contig.seq, file_handle);
			n_records++;
			return true;
		} else {
			return false;
		}
	}
	
	void close_file() {
		seqan::close(file_handle);
	}
};

struct kmer_t {
	kmer_2bit_t packed_f;
	kmer_2bit_t packed_rc;
	kmer_2bit_t packed_rep;
	bool valid; // does not contain any ambiguous bases
	
	void set_rep() {
		packed_rep = packed_f < packed_rc ? packed_f :packed_rc;
	}
	
	void print() {
			std::cout << "kmer f: " << packed_f << " rc: " << packed_rc <<  " valid " << valid << "\n"; 
	}
};

struct kmer_parser_t {
	std::string s; // sequence to parse
	int kmer_len;
	seq_t pos; // position in the sequence
	kmer_t kmer;
	bool first;

	void init(const std::string& seq, int k) {
		s = seq;
		kmer_len = k;
		pos = 0;
		first = true;
	}

	// compression
	inline void pack_init() {
 		kmer.valid = true;
		kmer.packed_f = 0;
		kmer.packed_rc = 0;
		for (int i = 0; i < kmer_len; i++) {
			uint8 c = dna5_table[(int) s[pos+i]];
			uint8 c_comp = dna5_comp[(int) c];
			if(c == BASE_IGNORE) {
				kmer.valid = false;
				break;
			}
			kmer.packed_f |= ((c & CHAR_MASK) << (BITS_IN_PACKED_WORD - (i+1) * BITS_PER_CHAR));
			kmer.packed_rc |= ((c_comp & CHAR_MASK) << (BITS_IN_PACKED_WORD - (kmer_len -  i) * BITS_PER_CHAR));
		}
		if(kmer.valid) {
			kmer.set_rep();
			pos += kmer_len;
		} else {
			pos++;
		}
	}

	void pack_roll() {
		uint8 c = dna5_table[(int) s[pos]];
		uint8 c_comp = dna5_comp[(int) c];
		if(c == BASE_IGNORE) {
			kmer.valid = false;
			pos++; // skip all the kmers containing this char
			return;
		}
		kmer.packed_f <<= BITS_PER_CHAR;
		kmer.packed_f |= ((c & CHAR_MASK) << (BITS_IN_PACKED_WORD - kmer_len * BITS_PER_CHAR));
		kmer.packed_rc >>= BITS_PER_CHAR;
		kmer.packed_rc |= ((c_comp & CHAR_MASK) << (BITS_IN_PACKED_WORD - BITS_PER_CHAR));
		kmer.packed_rc &= ~(CHAR_MASK << (BITS_IN_PACKED_WORD - (kmer_len + 1) * BITS_PER_CHAR));
		kmer.set_rep();
		pos++;
	}

	bool get_next_kmer(kmer_t& new_kmer) {
		if(pos >= s.size()) return false;
		// process the next char in the sequence
		if (first || !kmer.valid) {
			pack_init();
			if(first) {
				first = false;
			}
		} else {
			pack_roll();
		}
		new_kmer = kmer;
		return true;
	}

	 bool get_next_kmer(kmer_t& new_kmer, const int step_size) { // temp
                if(pos >= s.size()) return false;
               	seq_t init_pos = pos;
		pack_init();
                pos = init_pos + step_size;
		new_kmer = kmer;
                return true;
        }
};

// load the index file containing a list of FASTQ files (one per sample)
void load_panel_file(const std::string& fastq_index_fname, std::vector<std::string>& fastq_files) ;
void get_seq_kmers_packed(const std::string& seq_file, const int kmer_len, std::vector<kmer_2bit_t>& keys);
void get_seq_kmers_hashed(const std::string& seq_file, const int kmer_len, std::vector<kmer_2bit_t>& keys);
void store_kmers(const std::string& fname, const std::vector<kmer_2bit_t>& keys);
void load_kmers(const std::string& fname, std::vector<kmer_2bit_t>& keys);
void get_unique_kmers(std::vector<kmer_2bit_t>& keys);

#endif
