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
#include "index.h"

void print_usage(const index_params_t& params) {
	std::cout<<"Usage: ./gattaca [command] [options] \n";
	std::cout<<"\nCommand: index \n\n";
	std::cout<<"\t Options:\n\n";
	std::cout<<"\t		-i <arg>	sample FASTA/FASTQ file to index\n";
	std::cout<<"\t		-f <arg>	text file containing a list of sample FASTA/FASTQ files to index (one line per sample)\n";
	std::cout<<"\t		-k <arg>	length of the sequence kmers (default: " << params.k << ")\n";
	std::cout<<"\nCommand: lookup \n\n";
	std::cout<<"\t Options:\n\n";
	std::cout<<"\t		-c <arg>	input FASTA contig file [required]\n";
	std::cout<<"\t		-i <arg>	GAC index \n";
	std::cout<<"\t		-f <arg>	text file containing a list of GAC indices (one index file per line) \n";
	std::cout<<"\t		-o <arg>	output file name for the results [required]\n";
	std::cout<<"\t		-k <arg>	length of the sequence kmers (default: " << params.k << ")\n";
	std::cout<<"\t 		-m		output the median coverage of each contig in each sample (recommended) \n";
	std::cout<<"\t 		-a		output the average coverage of each contig in each sample \n";
	std::cout<<"\t 		-p		output the kmer counts in each contig by position (single sample index only, otherwise will use the first sample listed)  \n";
	std::cout<<"\nGeneral options:\n\n";
	std::cout<<"		-t <arg>	number of threads (default: " << params.n_threads << ")\n";
}
	
int main(int argc, char *argv[]) {
	index_params_t params;
	params.set_default_index_params();
	if (argc < 2) {
		print_usage(params);
		exit(1);
	}
	
	srand(1);
	if (strcmp(argv[1], "index") == 0) {
		std::string single_sample_fname;
		std::string panel_fname;
		int c;
		while ((c = getopt(argc-1, argv+1, "i:f:k:t:")) >= 0) {
			switch (c) {
				case 'i': single_sample_fname = std::string(optarg); break;
				case 'f': panel_fname = std::string(optarg); break;
				case 'k': params.k = atoi(optarg); break;
				case 't': params.n_threads = atoi(optarg); break;
				default: print_usage(params); return 0;
			}
		}
		#if defined(_OPENMP)
		omp_set_num_threads(params.n_threads);
		#endif

		std::vector<std::string> files_to_index;
		if(!single_sample_fname.empty()) {
			files_to_index.push_back(single_sample_fname);
		}
		if (!panel_fname.empty()) {
			load_panel_file(panel_fname, files_to_index);
		}
		if(files_to_index.empty()) {
			std::cout<<"[!]Please provide either the -i or the -f parameter.\n";
			print_usage(params);
			return 0;
		}
		index_t::build_and_save(files_to_index, params);

	} else if(strcmp(argv[1], "lookup") == 0) {
			std::string input_fasta_fname;
			std::string single_sample_fname;
			std::string panel_fname;
			std::string output_fname;
			bool median = false;
			bool avg = false;
			bool counts_by_pos = false;
			int c;
			while ((c = getopt(argc-1, argv+1, "i:c:f:o:t:k:mpa")) >= 0) {
				switch (c) {
					case 'c': input_fasta_fname = std::string(optarg); break;
					case 'i': single_sample_fname = std::string(optarg); break;
					case 'f': panel_fname = std::string(optarg); break;
					case 'o': output_fname = std::string(optarg); break;
					case 'k': params.k = atoi(optarg); break;
					case 'm': median = true; break;
					case 'a': avg = true; break;
					case 'p': counts_by_pos = true; break;
					case 't': params.n_threads = atoi(optarg); break;
					default: print_usage(params); return 0;
				}
			}
			if(input_fasta_fname.empty() || output_fname.empty()) {
				std::cout<<"[!]Parameters -c and -o are required.\n";
				print_usage(params);
				return 0;
			}

			if(!single_sample_fname.empty() && !panel_fname.empty()) {
				std::cout<<"[!]Please provide either the -r or the -f parameter, not both.\n";
				print_usage(params);
				return 0;
			}

			if(!median && !counts_by_pos && !avg) {
				std::cout << "Indicate which lookup operation you'd like to perform (median/avg/counts by position) \n";
				print_usage(params);
				return 0;
			}

			#if defined(_OPENMP)
			omp_set_num_threads(params.n_threads);
			#endif

			// load the index/indices
			std::vector<std::string> index_files;
			if(!single_sample_fname.empty()) {
				index_files.push_back(single_sample_fname);
			} else if (!panel_fname.empty()) {
				load_panel_file(panel_fname, index_files);
			} else {
				std::cout<<"[!]Please provide either the -i or the -f parameter.\n";
				print_usage(params);
				return 0;
			}
			index_t index;
			index.load(index_files, params);
			if(median || avg) index.get_aggregate_counts(input_fasta_fname, output_fname, median, params);
			if(counts_by_pos) index.get_counts_file(input_fasta_fname, 0, output_fname, params);
	} else if(strcmp(argv[1], "count2char") == 0) {
		// only mphf suported now
		//params.alg = index_alg::EXACT_MPHF;
		//index_t idx;
		//idx.load(argv[optind+1], params);
		//idx.save_index(argv[optind+2], 8, params);
	} else if(strcmp(argv[1], "discretize") == 0) {
		std::string in_index_fname;
		std::string out_index_fname; 
		int nbits = 0;
		int min = 0; 
		int max = 255;
		int max_v = 255;
		int c;
		while ((c = getopt(argc-1, argv+1, "i:o:n:t:m:M:V:")) >= 0) {
			switch (c) {
				case 'i': in_index_fname = std::string(optarg); break;
				case 'o': out_index_fname = std::string(optarg); break;
				case 'n': nbits = atoi(optarg); break;
				case 'm': min = atoi(optarg); break;
				case 'M': max = atoi(optarg); break;
				case 'V': max_v = atoi(optarg); break;
				case 't': params.n_threads = atoi(optarg); break;
				default: print_usage(params); return 0;
			}
		}
		if(in_index_fname.empty() || out_index_fname.empty()) {
			print_usage(params); 
			return 0;
		}
		index_t idx;
		idx.discretize(in_index_fname, out_index_fname, nbits, min, max, max_v, params);
	} else {
		print_usage(params);
		exit(1);
	}
	return 0;
}
