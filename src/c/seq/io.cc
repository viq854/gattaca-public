#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <time.h>
#include <limits.h>
#include <fstream>
#include <sstream>
#include "types.h"
#include "io.h"
#include "city.h"

// load the panel file containing a list of sample sequence/index files (one per sample)
void load_panel_file(const std::string& panel_fname, std::vector<std::string>& files) {
	std::ifstream file;
	file.open(panel_fname.c_str(), std::ios::in);
	if (!file.is_open()) {
		std::cerr << "load_panel_file: Cannot open the panel file: " << panel_fname << "!\n";
		exit(1);
	}
	std::string fname;
	while(std::getline(file, fname)) {
		files.push_back(fname);
	}
	file.close();
	std::cout<< "Panel file contains " << files.size() << " entries. \n";
}

void load_panel_file(const std::string& panel_fname, std::vector<std::vector<std::string>>& files) {
	std::ifstream file;
	file.open(panel_fname.c_str(), std::ios::in);
	if (!file.is_open()) {
		std::cerr << "load_panel_file: Cannot open the panel file: " << panel_fname << "!\n";
		exit(1);
	}
	std::string line;
	while(std::getline(file, line)) {
	    std::stringstream file_names(line);
	    std::vector<std::string> per_sample_files;
	    std::string fname;
	    while(std::getline(file_names, fname, '\t')) {
		    per_sample_files.push_back(fname);
	    }
	    files.push_back(per_sample_files);
	}
	file.close();
	std::cout<< "Panel file contains " << files.size() << " entries. \n";
}

void get_seq_kmers_packed(const std::string& seq_file, const int kmer_len, std::vector<kmer_2bit_t>& keys) {
	seq_file_reader_t reader;
	reader.open_file(seq_file, FASTA);
	read_t r;
	while(reader.load_next_read(r)) {
		kmer_parser_t seq_parser;
		seq_parser.init(r.seq, kmer_len);
		kmer_t kmer;
		while(seq_parser.get_next_kmer(kmer)) {
			if(!kmer.valid) continue;
			keys.push_back(kmer.packed_rep);
		}
	}
	reader.close_file();
}

void get_seq_kmers_hashed(const std::string& seq_file, const int kmer_len, std::vector<kmer_2bit_t>& keys) {
	keys.reserve(100000000);
	seq_file_reader_t reader;
	std::cout << "Opening file " << seq_file << "\n";
	reader.open_file(seq_file, FASTQ);
	read_t r;
	int i = 0;
	while(reader.load_next_read(r)) {
		i++;
		for(int j = 0; j < (int) r.seq.size() - kmer_len + 1; j++) {
			kmer_2bit_t kmer_hash = CityHash64(&r.seq[j], kmer_len);
			if(kmer_hash != 0) {
				keys.push_back(kmer_hash);
			}
		}
	}
	std::cout << "Loaded " << i << " reads\n";
	reader.close_file();
}

void get_unique_kmers(std::vector<kmer_2bit_t>& keys) {
	std::sort(keys.begin(), keys.end());
	auto last = std::unique(keys.begin(), keys.end());
	keys.erase(last, keys.end());
}

void load_kmers(const std::string& fname, std::vector<kmer_2bit_t>& keys) {
	std::ifstream file;
	file.open(fname.c_str(), std::ios::in | std::ios::binary);
	uint64 n_keys;
	file.read(reinterpret_cast<char*>(&n_keys), sizeof(n_keys));
	keys.resize(n_keys);
	file.read(reinterpret_cast<char*>(&(keys[0])), n_keys*sizeof(kmer_2bit_t));
	file.close();
}

void store_kmers(const std::string& fname, std::vector<kmer_2bit_t>& keys) {
	std::ofstream file;
	file.open(fname.c_str(), std::ios::out | std::ios::binary);
	uint64 n_keys = keys.size();
	file.write(reinterpret_cast<char*>(&n_keys), sizeof(n_keys));
	for(uint64 j = 0; j < n_keys; j++) {
		file.write(reinterpret_cast<const char*>(&keys[j]), sizeof(kmer_2bit_t));
	}
	file.close();
}
