GATTACA: Framework for lightweight metagenomic binning
============

### About

GATTACA is a framework for rapid and accurate binning of metagenomic contigs, which (1) does not require read mapping and (2) enables efficient stand-alone analysis of single metagenomic samples. 

At a high level, GATTACA quickly estimates co-abundance profiles within a panel of reference metagnomes (used as features during clustering) from kmer counts stored in a compact index. This results in a significant speedup in coverage estimation. It also provides a way to index metagenomic samples (e.g. from public repositories) and easily reuse them across experiments. Once the kmer count index (GAC) of a sample is computed, this sample can be efficiently re-used in other studies: the small size of the index (~150MB) allows it to be easily downloaded and shared; and the abundances can be computed via fast lookups.

Leveraging the MinHash technique to quickly compare metagenomic samples, GATTACA also provides an efficient way to identify metagenomic samples (e.g. from publicly repositories) that can be incorporated into the panel of reference metagenomes to further improve binning accuracy. More specifically, GATTACA computes small MinHash fingerprints (GAFs) for each sample and offers functionality to quickly compare the fingerprints and select appropriate samples based on relevance and diversity criteria. 

For more information about the algorithm please see the following publication:  
*GATTACA: Lightweight Metagenomic Binning With Compact Indexing Of Kmer Counts And MinHash-based Panel Selection*. Popic V, Kuleshov V, Snyder M, Batzoglou S.(RECOMB 2017)  
bioRxiv 130997; doi: https://doi.org/10.1101/130997   

### Installation

```
git clone --recursive https://github.com/viq854/gattaca-public.git
cd src/c
make all
```

Internally, the tool uses the following three libraries: ```seqan``` (as submodule), ```libbf``` (as source), and ```cmph``` (as submodule). These are automatically downloaded and installed. Please make sure to specify the ```--recursive ``` flag when cloning.

#### System Requirements

To build the third-party libraries and the GATTACA C++ source code:  
```GCC >= 4.7```  
```CMake >= 2.8```  
```OpenMP``` (optional)

To run the clustering code (Python):  
```Python 2.7```  
```NumPy, SciPy, Scikit-Learn```

*In case there are any issues when installing on your system, please let us know and we'll address them quickly. This is our first version of the tool and it hasn't been tested yet across different platforms. Thanks!*

### Usage

Two executables are provided: ```gattaca``` to build the GAC indices and estimate coverage and ```gattaca-minh``` to compare samples using MinHash.

#### GAC index construction and coverage estimation
```
./gattaca [command] [options] 

Command: index 
	 Options:
			-i <arg>	sample FASTA/FASTQ file to index
			-f <arg>	text file containing a list of sample FASTA/FASTQ files to index (one line per sample)
			-k <arg>	length of the sequence kmers (default: 31)

Command: lookup 
	 Options:
			-c <arg>	input FASTA contig file [required]
			-i <arg>	GAC index 
			-f <arg>	text file containing a list of GAC indices (one index file per line) 
			-o <arg>	output file name for the results [required]
			-k <arg>	length of the sequence kmers (default: 31)
	 		-m		output the median coverage of each contig in each sample (recommended) 
	 		-a		output the average coverage of each contig in each sample 
	 		-p		output the kmer counts in each contig by position (single sample index only, otherwise will use the first sample listed)  

General options:
		-t <arg>	number of threads (default: 1)
```

Basic workflow: (1) build the kmer count indices using the ```index``` command for the samples in the reference panel (will create one ```.gatc``` file for each sample)  and (2) estimate the contig coverages in the panel using the ```lookup``` command. The output file for the aggregate (median and average flags) is a tab-delimited text file, where row *i* corresponds to contig *i* in the input FASTA file and  column *j* is the coverage in sample *j* of the GAC index list.


#### Clustering

To cluster the contigs, please execute the following script. 

```
python $(GATTACA_HOME)/src/python/gattaca.py cluster \
		--contigs contigs.fasta \
          	--coverage coverage_inputtable.tsv \
	  	--algorithm dirichlet \
          	--clusters output.clusters
          
```

The required ```--coverage``` file is a tab-delimited text file summarizing the co-abundance profiles across the sample panel. It can be created from the output file of the ```lookup``` command above, by adding the header consisting of the sample names and adding the contig name and length to each row. It should have the following format:

```
contig	length	Sample1	... SampleN
contig-name1	1000	5 ...	12	
...
```

The clustering results are written to ```output.clusters```.

#### Sample Comparisongs with MinHash

The ```gattaca-minh``` executable is used to compare samples using MinHash.
```
./gattaca-minh [command] [options] 

Command: fp 
	 Options:
			-i <arg>	sample FASTA/FASTQ file to fingerprint with MinHash 
			-f <arg>	text file containing a list of sample FASTA/FASTQ files to fingerprint (one line per sample)
			-k <arg>	length of the sequence kmers (default: 31)
			-L <arg>	length of the fingerprint (default: 1024)

Command: relevance 
	 Options:
			-q <arg>	MinHash fingerprint of the query sample [required]
			-i <arg>	MinHash fingperint index of the reference panel 
			-f <arg>	text file containing a list of reference MinHash fingerprints
			-L <arg>	length of the fingerprint (default: 1024)

Command: diversity 
	 Options:
			-f <arg>	text file containing a list of reference MinHash fingerprints
			-L <arg>	length of the fingerprint (default: 1024)
			-N <arg>	number of samples to select, must be smaller than the number of samples in the panel [required]

Command: fp_index 
	 Options:
			-f <arg>	text file containing a list of reference MinHash fingerprints
			-o <arg>	output filename 
			-L <arg>	length of the fingerprint (default: 1024)

Command: compare_pairwise_fp 
	 Options:
			-f <arg>	text file containing a list of reference MinHash fingerprints
			-L <arg>	length of the fingerprint (default: 1024)

General options:
		-t <arg>	number of threads (default: 1)
```

### Citation
*GATTACA: Lightweight Metagenomic Binning With Compact Indexing Of Kmer Counts And MinHash-based Panel Selection*. Popic V, Kuleshov V, Snyder M, Batzoglou S.(RECOMB 2017)  
bioRxiv 130997; doi: https://doi.org/10.1101/130997   

### License
MIT License

### Support
For help running the program or any questions/suggestions/bug reports, please contact viq@stanford.edu.
