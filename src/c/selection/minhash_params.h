#ifndef INDEX_H_
#define INDEX_H_

#pragma once

// program parameters
struct minhash_params_t {
    int k;				// length of the sequence kmers
    int n_threads;		// multi-threading
    int L; // length of sample fingerprint
	int N; // number of samples to select

    void set_default_params() {
		k = 31;
		n_threads = 1;
		L = 1024;
    }
};

#endif
