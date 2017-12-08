#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <limits.h>
#include <fstream>
#include "../seq/io.h"
#include "index.h"

void index_t::build_and_save(const std::vector<std::vector<std::string>>& files_to_index, index_params_t& params) {
	build_and_save(files_to_index, std::numeric_limits<counter_t>::digits, params);
}

void index_t::build_and_save(const std::vector<std::vector<std::string>>& files_to_index, const int n_count_bits, index_params_t& params) {
	std::string index_fname_suffix(".k");
	index_fname_suffix += std::to_string(params.k);
	index_fname_suffix += std::string(".gatc"); // TODO: enforce identical params upon loading

	//double t = omp_get_wtime();
	#if defined(_OPENMP)
	#pragma omp parallel for
	#endif
	for(unsigned int i = 0; i < files_to_index.size(); i++) { // over samples
	    std::cout << "Processing sample " << i << "\n";
		counts_table_t* table = count_mphf(files_to_index[i], params);
		table->print_stats();
		std::string index_fname(files_to_index[i][0].c_str());
		index_fname += index_fname_suffix;
		table->save_to_file(index_fname, n_count_bits);
		delete table;
	}
	//printf("Index construction time (total): %.2f sec\n", omp_get_wtime() - t);
}

void index_t::load(const std::vector<std::string>& index_files, index_params_t& params) {
	n_samples = index_files.size();
	count_tables = new counts_table_t*[n_samples];
	for(int i = 0; i < n_samples; i++) {
		std::cout << "Loading sample " << i << " from file " << index_files[i] << "\n";
		count_tables[i] = init_counts_table(params);
		count_tables[i]->load_from_file(index_files[i], 0);
		//count_tables[i]->print_stats();
	}
}

counts_table_t* index_t::init_counts_table(index_params_t& params) {
	switch(params.alg) { // TODO: store index type and info in file
		case EXACT_MPHF:
			return new mphf_table_t();
		case CMS: // TODO: separate discrete loading
			return new cms_table_t<BITS_PER_COUNTER>();
		case EXACT_BASIC: 
			return new basic_table_t(1);
		case MINHASH: 
			return new minhash_table_t();
		default: 
			std::cout << "Invalid index type option \n"; exit(-1);
	}
}

counts_table_t* index_t::count_mphf(const std::string& seq_file, const index_params_t& params) {
    std::vector<std::string> per_sample_files;
    per_sample_files.push_back(seq_file);
    return count_mphf(per_sample_files, params);
}

counts_table_t* index_t::count_mphf(const std::vector<std::string>& seq_files, const index_params_t& params) {
	std::vector<kmer_2bit_t> keys;
	keys.reserve(100000000);

    for(unsigned int i = 0; i < seq_files.size(); i++) {
        std::cout << "... using file: " << seq_files[i] << "\n";
	    get_seq_kmers_packed(seq_files[i], params.k, keys);
	}
	return new mphf_table_t(keys);
}

counts_table_t* index_t::count_minhash(const std::string& seq_file, const index_params_t& params) {
	minhash_table_t* counts_table = new minhash_table_t(params.k, 128, 48, 2, 1);
	seq_file_reader_t reader;
	reader.open_file(seq_file);
	read_t r;
	while(reader.load_next_read(r)) {
		counts_table->insert(r.seq, 0);
	}
	counts_table->done();
	reader.close_file();
	return counts_table;
}

// count the kmers in the given sample file and return the counts table
counts_table_t* index_t::count(const std::string& seq_file, const index_params_t& params)  {
	if(params.alg == EXACT_MPHF) {
		return count_mphf(seq_file, params);
	} else if(params.alg == MINHASH) {
		return count_minhash(seq_file, params);
	}
	counts_table_t* counts_table;
	if(params.alg == CMS) { 
			if(params.cms_config.epsilon != 0 && params.cms_config.delta != 0) { // initialize the CMS index based on desired error parameters
				counts_table = cms_table_t<BITS_PER_COUNTER>::err_params(params.cms_config.delta, params.cms_config.epsilon, 1);
			} else {
				counts_table = cms_table_t<BITS_PER_COUNTER>::size_params(params.cms_config.n_rows, params.cms_config.n_cols, 1);
			} 
	} else {
		counts_table = new basic_table_t(1);
	}

	seq_file_reader_t reader;
	reader.open_file(seq_file);
	read_t r;
	while(reader.load_next_read(r)) {
		kmer_parser_t seq_parser;
		seq_parser.init(r.seq, params.k);
		kmer_t kmer;
		while(seq_parser.get_next_kmer(kmer)) {
			if(!kmer.valid) continue;
			counts_table->insert(kmer.packed_rep, 0);
		}
	}
	reader.close_file();
	return counts_table;
}

// outputs kmer frequencies for each input contig in each indexed sample (1 line / sample)
void index_t::get_counts_file(const std::string& fasta_fname, const int sample_id, const std::string& output_fname, const index_params_t& params) {
	std::ofstream file;
	file.open(output_fname.c_str());
	if(!file.is_open()) {
		std::cerr << "ERROR: Could not open file: " << output_fname << "\n";
		exit(-1);
	}
	seq_file_reader_t reader;
	reader.open_file(fasta_fname);
	contig_t c;
	while(reader.load_next_contig(c)) {
		if(params.alg == MINHASH) {
			std::cerr << "Operation not supported yet. \n";
		} else {
			kmer_parser_t seq_parser;
			seq_parser.init(c.seq, params.k);
			kmer_t kmer;
			while(seq_parser.get_next_kmer(kmer)) {
				if(!kmer.valid) continue;
				counter_t count = count_tables[sample_id]->lookup(kmer.packed_rep, 0);
				file << (int) count << " ";
			}
			file << "\n";
		}
	}
	reader.close_file();
	file.close();
}

counter_t index_t::get_avg_count(const std::string& seq, const int sample_id, const index_params_t& params) const {
	uint64 sum = 0;
	uint64 n_counts = 0;
	kmer_parser_t seq_parser;
	seq_parser.init(seq, params.k);
	kmer_t kmer;
	while(seq_parser.get_next_kmer(kmer)) {
		if(!kmer.valid) continue;
		counter_t count = count_tables[sample_id]->lookup(kmer.packed_rep, sample_id);
		if(count <= params.count_cutoff_max  && count >= params.count_cutoff_min) {
			sum += count;
			n_counts += 1;
		}
	}
	if(n_counts == 0) return 0;
	return 1.0f*sum/n_counts;
}

counter_t index_t::get_median_count(const std::string& seq, const int sample_id, const index_params_t& params) const {
	std::vector<counter_t> kmer_counts;
	uint64 n_zeros = 0;
	
	kmer_parser_t seq_parser;
	seq_parser.init(seq, params.k);
	kmer_t kmer;
	while(seq_parser.get_next_kmer(kmer)) {
		if(!kmer.valid) continue;
		counter_t count = count_tables[sample_id]->lookup(kmer.packed_rep, sample_id);
		if(count == 0) {
			n_zeros++;
		} else if(count <= params.count_cutoff_max  && count >= params.count_cutoff_min) {
			kmer_counts.push_back(count);
		}
	}
	if(kmer_counts.size() == 0) return 0;
	//std::sort(kmer_counts.begin(), kmer_counts.end());
	const uint64 n_counts = n_zeros + kmer_counts.size();
	const uint64 median_idx = n_counts/2;
	if(median_idx < n_zeros) return 0;
	const uint64 median_arr_idx = median_idx - n_zeros;
	std::nth_element(kmer_counts.begin(), kmer_counts.begin() + median_arr_idx, kmer_counts.end());
	return kmer_counts[median_arr_idx];
}

void index_t::get_aggregate_counts(const std::string& fasta_fname, const std::string& output_fname, const bool is_median,
		const index_params_t& params) {
	//double t = omp_get_wtime();
	std::cout << "Generating aggregate coverage for " << fasta_fname << " , # samples: " << n_samples << "\n";
	std::ofstream file;
	file.open(output_fname.c_str());
	if(!file.is_open()) {
		std::cerr << "ERROR: Could not open file: " << output_fname << "\n";
		exit(-1);
	}
	seq_file_reader_t reader;
	reader.open_file(fasta_fname);
	int idx = 0;
	contig_t c;
	while(reader.load_next_contig(c)) {
		std::vector<counter_t> sample_counts(n_samples);
		#if defined(_OPENMP)
        	#pragma omp parallel for
        	#endif
		for(int sample_id = 0; sample_id < n_samples; sample_id++) {
			if(is_median) {
				sample_counts[sample_id] = get_median_count(c.seq, sample_id, params);
			} else {
				sample_counts[sample_id] = get_avg_count(c.seq, sample_id, params);
			}
		} 
		for(int sample_id = 0; sample_id < n_samples; sample_id++) {
			file << (int) sample_counts[sample_id] << "\t";
		}
		file << "\n";
		idx++;
		if(idx % 10000 == 0) {
			std::cout << "Processed " << idx << " contigs\n";
		}
	}
	if(idx % 10000 != 0) {
		std::cout << "Processed " << idx << " contigs total\n";
	}
	reader.close_file();
	file.close();
	//printf("Lookup time (total): %.2f sec\n", omp_get_wtime() - t);
}

void index_t::discretize(const std::string& in_index_fname, const std::string& out_index_fname, const int nbits, const int min_r, const int max_r, const int max_v, index_params_t& params) {
//	load_index(in_index_fname, 0, params);
//	std::ofstream file;
//	file.open(out_index_fname.c_str(), std::ios::out | std::ios::binary);
//	if (!file.is_open()) {
//		std::cerr << "ERROR: Could not open file: " << out_index_fname << "\n";
//		exit(1);
//	}
//	counts_table_t* d_counts_table;
//	switch(nbits){
//			case 4: d_counts_table = new discretized_cms_table_t<4>((cms_table_t<BITS_PER_COUNTER>*) counts_table, min_r, max_r, max_v); break;
//			case 5: d_counts_table = new discretized_cms_table_t<5>((cms_table_t<BITS_PER_COUNTER>*) counts_table, min_r, max_r, max_v); break;
//			case 6: d_counts_table = new discretized_cms_table_t<6>((cms_table_t<BITS_PER_COUNTER>*) counts_table, min_r, max_r, max_v); break;
//			case 7: d_counts_table = new discretized_cms_table_t<7>((cms_table_t<BITS_PER_COUNTER>*) counts_table, min_r, max_r, max_v); break;
//			default: exit(1);
//	}
//	d_counts_table->save_to_file(file);
//	file.close();
}
