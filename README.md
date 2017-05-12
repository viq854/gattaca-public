GATTACA: Frameworkd for lightweight metagenomic binning
============

### About

GATTACA is a framework for rapid and accurate binning of metagenomic contigs, which (1) does not require read mapping and (2) enables efficient stand-alone analysis of single metagenomic samples. 

At a high level, GATTACA quickly estimates co-abundance profiles within a panel of reference metagnomes (used as features during clustering) from kmer counts stored in a compact index. This results in a significant speedup in coverage estimation. It also provides a way to index metagenomic samples (e.g. from public repositories) and easily reuse them across experiments. Once the kmer count index (GAC) of a sample is computed, this sample can be efficiently re-used in other studies: the small size of the index (~150MB) allows it to be easily downloaded and shared; and the abundances can be computed via fast lookups.

Leveraging the MinHash technique to quickly compare metagenomic samples, GATTACA also provides an efficient way to identify metagenomic samples (e.g. from publicly repositories) that can be incorporated into the panel of reference metagenomes to further improve binning accuracy. More specifically, GATTACA computes small MinHash fingerprints (GAFs) for each sample and offers functionality to quickly compare the fingerprints and select appropriate samples based on relevance and diversity criteria. 

### Installation

git clone --recursive 

Internally, the tool uses the following three libraries: seqan (as submodule), libbf (as source), and cmph (as submodule). These should be automatically downloaded and installed.

[In case there are any issues when installing on your system, please let us know and we'll address them quickly. This is our first version of the tool and it hasn't been tested yet across different platforms. Thanks!]

### License

MIT License

### Support

For help running the program or any questions/suggestions/bug reports, please contact viq@stanford.edu.
