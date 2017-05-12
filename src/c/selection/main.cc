#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <assert.h>
#include "../seq/types.h"
#include "../seq/io.h"
#include "minhash.h"
#include "minhash_params.h"

void print_usage(const minhash_params_t& params) {
	std::cout<<"Usage: ./gattaca [command] [options] \n";
	std::cout<<"\nCommand: fp \n\n";
	std::cout<<"\t Options:\n\n";
	std::cout<<"\t		-i <arg>	sample FASTA/FASTQ file to fingerprint with MinHash \n";
	std::cout<<"\t		-f <arg>	text file containing a list of sample FASTA/FASTQ files to fingerprint (one line per sample)\n";
	std::cout<<"\t		-k <arg>	length of the sequence kmers (default: " << params.k << ")\n";
	std::cout<<"\t		-L <arg>	length of the fingerprint (default: " << params.L << ")\n";
	std::cout<<"\nCommand: relevance \n\n";
	std::cout<<"\t Options:\n\n";
	std::cout<<"\t		-q <arg>	MinHash fingerprint of the query sample [required]\n";
	std::cout<<"\t		-i <arg>	MinHash fingperint index of the reference panel \n";
	std::cout<<"\t		-f <arg>	text file containing a list of reference MinHash fingerprints\n";
	std::cout<<"\t		-L <arg>	length of the fingerprint (default: " << params.L << ")\n";
	std::cout<<"\nCommand: diversity \n\n";
	std::cout<<"\t Options:\n\n";
	std::cout<<"\t		-f <arg>	text file containing a list of reference MinHash fingerprints\n";
	std::cout<<"\t		-L <arg>	length of the fingerprint (default: " << params.L << ")\n";
	std::cout<<"\t		-N <arg>	number of samples to select, must be smaller than the number of samples in the panel [required]\n";
	std::cout<<"\nCommand: fp_index \n\n";
	std::cout<<"\t Options:\n\n";
	std::cout<<"\t		-f <arg>	text file containing a list of reference MinHash fingerprints\n";
	std::cout<<"\t		-o <arg>	output filename \n";
	std::cout<<"\t		-L <arg>	length of the fingerprint (default: " << params.L << ")\n";
	std::cout<<"\nCommand: compare_pairwise_fp \n\n";
	std::cout<<"\t Options:\n\n";
	std::cout<<"\t		-f <arg>	text file containing a list of reference MinHash fingerprints\n";
	std::cout<<"\t		-L <arg>	length of the fingerprint (default: " << params.L << ")\n";
	std::cout<<"\nGeneral options:\n\n";
	std::cout<<"		-t <arg>	number of threads (default: " << params.n_threads << ")\n";
}
	
int main(int argc, char *argv[]) {
	minhash_params_t params;
	params.set_default_params();
	if (argc < 2) {
		print_usage(params);
		exit(1);
	}
	
	srand(1);
	if(strcmp(argv[1], "fp") == 0) {
		std::string single_sample_fname;
		std::string panel_fname;
		int c;
		while ((c = getopt(argc-1, argv+1, "i:f:k:t:L:")) >= 0) {
			switch (c) {
				case 'i': single_sample_fname = std::string(optarg); break;
				case 'f': panel_fname = std::string(optarg); break;
				case 'k': params.k = atoi(optarg); break;
				case 't': params.n_threads = atoi(optarg); break;
				case 'L': params.L = atoi(optarg); break;
				default: print_usage(params); return 0;
			}
		}
		#if defined(_OPENMP)
		omp_set_num_threads(params.n_threads);
		#endif

		std::vector<std::string> files_to_fp;
		if(!single_sample_fname.empty()) {
			files_to_fp.push_back(single_sample_fname);
		}
		if(!panel_fname.empty()) {
			load_panel_file(panel_fname, files_to_fp);
		}
		if(files_to_fp.empty()) {
			std::cout<<"[!]No samples specified. Please provide either the -i or the -f parameter or both.\n";
			print_usage(params);
			return 0;
		}
		minhash_t minhash(params.k, params.L);
		#if defined(_OPENMP)
		#pragma omp parallel for
		#endif
		for(unsigned int i = 0; i < files_to_fp.size(); i++) {
			std::cout << "Processing sample " << i << " from file " << files_to_fp[i] << "\n";
			std::vector<kmer_2bit_t> keys;
			get_seq_kmers_hashed(files_to_fp[i], params.k, keys);
			get_unique_kmers(keys);
			minhash_fp_t fingerprint;
			minhash.compute_fp(keys, fingerprint);
			minhash_t::save_fp_to_file(files_to_fp[i], params.k, fingerprint);
		}
		printf("Minhash fingerprint construction done!\n");
	} else if(strcmp(argv[1], "relevance") == 0) {
		std::string query_fname;
		std::string minhash_index_fname;
		std::string panel_fname;
		int c;
		while ((c = getopt(argc-1, argv+1, "q:f:i:t:L:")) >= 0) {
			switch (c) {
				case 'q': query_fname = std::string(optarg); break;
				case 'f': panel_fname = std::string(optarg); break;
				case 'i': minhash_index_fname = std::string(optarg); break;
				case 't': params.n_threads = atoi(optarg); break;
				case 'L': params.L = atoi(optarg); break;
				default: print_usage(params); return 0;
			}
		}
		if(query_fname.empty() || (minhash_index_fname.empty() && panel_fname.empty())) {
			print_usage(params);
			return 0;
		}
		#if defined(_OPENMP)
		omp_set_num_threads(params.n_threads);
		#endif

		minhash_fp_t query_fp;
		minhash_t::load_fp_from_file(query_fname, params.L, query_fp);
		if(!minhash_index_fname.empty()) {
			minhash_t minhash_index(params.k, params.L);
			minhash_index.load_index_from_file(minhash_index_fname);
			std::vector<sample_support_t> id_counts;
			minhash_index.lookup(query_fp, id_counts);
			for(unsigned int i = 0; i < id_counts.size(); i++) {
				std::cout << "Sample " << id_counts[i].id+1 << ":\t" << id_counts[i].support << "\n";
			}
		} else {
			std::vector<std::string> files_to_compare;
			load_panel_file(panel_fname, files_to_compare);
			for(unsigned int i = 0; i < files_to_compare[i].size(); i++) {
				minhash_fp_t fp;
				minhash_t::load_fp_from_file(files_to_compare[i], params.L, fp);
				if(fp.v.size() != query_fp.v.size()) {
					std::cout<<"[!]The fingerprint lengths do not match.\n";
					return 0;
				}
				std::cout << "Sample " << i+1 << ":\t" << minhash_t::compare_fps(query_fp, fp) << "\n";
			}
		}
	} else if(strcmp(argv[1], "fp_index") == 0) {
		std::string panel_fname;
		std::string output_fname;
		int c;
		while ((c = getopt(argc-1, argv+1, "f:o:t:L:")) >= 0) {
			switch (c) {
				case 'f': panel_fname = std::string(optarg); break;
				case 'o': output_fname = std::string(optarg); break;
				case 't': params.n_threads = atoi(optarg); break;
				case 'L': params.L = atoi(optarg); break;
				default: print_usage(params); return 0;
			}
		}
		if(panel_fname.empty() || output_fname.empty()) {
			print_usage(params);
			return 0;
		}
		std::vector<std::string> files_to_index;
		load_panel_file(panel_fname, files_to_index);

		std::cout << "Indexing the MinHash FPs....\n";
		minhash_t minhash_index(params.k, params.L);
		minhash_index.init_index_tables();
		for(unsigned int i = 0; i < files_to_index.size(); i++) {
			minhash_fp_t fp;
			minhash_t::load_fp_from_file(files_to_index[i], params.L, fp);
			minhash_index.insert(fp, i);
		}
		minhash_index.save_index_to_file(output_fname);
		printf("Minhash fingerprint indexing done!\n");

	} else if(strcmp(argv[1], "diversity") == 0) {
		std::string panel_fname;
		int num_samples = 0; // number of samples to select
		int c;
		while ((c = getopt(argc-1, argv+1, "f:t:L:N:")) >= 0) {
			switch (c) {
				case 'f': panel_fname = std::string(optarg); break;
				case 'N': num_samples = atoi(optarg); break;
				case 't': params.n_threads = atoi(optarg); break;
				case 'L': params.L = atoi(optarg); break;
				default: print_usage(params); return 0;
			}
		}
		if(panel_fname.empty()) {
			print_usage(params);
			return 0;
		}

		// load the fingerprints
		std::vector<std::string> files_to_compare;
		load_panel_file(panel_fname, files_to_compare);
		std::vector<minhash_fp_t> fingerprints(files_to_compare.size());
		for(unsigned int i = 0; i < files_to_compare.size(); i++) {
			minhash_t::load_fp_from_file(files_to_compare[i], params.L, fingerprints[i]);
		}

		// compute the distances
		std::vector<std::vector<int>> distances(fingerprints.size());
		for(unsigned int i = 0; i < fingerprints.size(); i++) {
			distances[i].resize(fingerprints.size());
		}
		int max_i = 0;
		int max_j = 0;
		int max_d = 0;
		std::vector<int> ids;
		for(unsigned int i = 0; i < fingerprints.size(); i++) {
			ids.push_back(i);
			for(unsigned int j = 0; j < fingerprints.size(); j++) {
				if(i == j) continue;
				int distance = 0;
				for(int k = 0; k < params.L; k++) {
					if(fingerprints[i].v[k] == fingerprints[j].v[k]) {
							distance++;
					}
				}
				distances[i][j] = params.L - distance;
				if(distances[i][j] > max_d) {
					max_d = distances[i][j];
					max_i = i;
					max_j = j;
				}
			}
		}

		// MaxMin algorithm
		std::vector<int> selection;
		selection.push_back(max_i);
		selection.push_back(max_j);
		ids.erase(std::remove(ids.begin(), ids.end(), max_i), ids.end());
		ids.erase(std::remove(ids.begin(), ids.end(), max_j), ids.end());
		while((int) selection.size() < num_samples) {
			int max_min = 0;
			int max_min_id = 0;
			for(unsigned int i = 0; i < ids.size(); i++) {
				const int id = ids[i];
				int min_dist = distances[id][selection[0]];
				for(unsigned int j = 1; j < selection.size(); j++) {
					const int dist = distances[id][selection[j]];
					if(dist < min_dist) {
						min_dist = dist;
					}
				}
				if(min_dist > max_min) {
					max_min = min_dist;
					max_min_id = id;
				}
			}
			selection.push_back(max_min_id);
			ids.erase(std::remove(ids.begin(), ids.end(), max_min_id), ids.end());
		}
		for(unsigned int i = 0; i < selection.size(); i++) {
			std::cout << i << ": " << selection[i] << "\n";
		}
	} else if(strcmp(argv[1], "compare_fp") == 0) {
		std::string fp1_fname;
		std::string fp2_fname;
		int c;
		while ((c = getopt(argc-1, argv+1, "1:2:L:")) >= 0) {
			switch (c) {
				case '1': fp1_fname = std::string(optarg); break;
				case '2': fp2_fname = std::string(optarg); break;
				case 'L': params.L = atoi(optarg); break;
				default: print_usage(params); return 0;
			}
		}
		// load the fingerprints
		minhash_fp_t fp1;
		minhash_t::load_fp_from_file(fp1_fname, params.L, fp1);
		minhash_fp_t fp2;
		minhash_t::load_fp_from_file(fp2_fname, params.L, fp2);
		std::cout << "Number of shared entries: " << minhash_t::compare_fps(fp1, fp2) << "\n";

		// load kmers from file....
		/*uint64 idx1 = 0;
		uint64 idx2 = 0;
		while(idx1 < keys1.size() && idx2 < keys2.size()) {
			if(keys1[idx1] < keys2[idx2]) {
				idx1++;
			} else if (keys1[idx1] > keys2[idx2]) {
				idx2++;
			} else {
				count++;
				idx1++;
				idx2++;
			}
		}	
		std::cout << "" << count << "\n";
		std::cout << "" << n_keys << " " << n_keys2 << " " << (float) count*1.0f/(n_keys + n_keys2 - count) << "\n";*/
	} else if(strcmp(argv[1], "compare_pairwise_fp") == 0) {
		std::string panel_fname;
		int c;
		while ((c = getopt(argc-1, argv+1, "f:L:")) >= 0) {
				switch (c) {
						case 'f': panel_fname = std::string(optarg); break;
						case 'L': params.L = atoi(optarg); break;
						default: print_usage(params); return 0;
				}
		}
		std::vector<std::string> fp_files;
		load_panel_file(panel_fname, fp_files);
		const int n_samples = fp_files.size();
		std::vector<minhash_fp_t> fingerprints(n_samples);
        for(int i = 0; i < n_samples; i++) {
        	minhash_t::load_fp_from_file(fp_files[i], params.L, fingerprints[i]);
		}
		for(int s1 = 0; s1 < n_samples; s1++) {
			for(int s2 = 0; s2 < n_samples; s2++) {
				int count = 0;
				for(int i = 0; i < params.L; i++) {
					if(fingerprints[s1].v[i] == fingerprints[s2].v[i]) {
							count++;
					}
				}
				if(s2 == 0) {
					std::cout << count;
				} else {
					std::cout << "\t" << count;
				}
			}
			std::cout << "\n";
		}	
	} else {
		print_usage(params);
		exit(1);
	}
	return 0;
}
