#ifndef TYPES_H_
#define TYPES_H_

#pragma once
#include <string>
#include <set>
#include <map>
#include <queue>
#include <vector>

// general
typedef unsigned int uint32;
typedef unsigned long long int uint64;
typedef unsigned short uint16;
typedef unsigned char uint8;
#define BITS_IN_WORD 	32
#define BITS_IN_LWORD 	64

// sequences
typedef uint32 seq_t;

// kmer packing
typedef uint64 kmer_2bit_t;
typedef uint64 kmer_hash_t;
#define BITS_PER_CHAR 2
#define BITS_IN_PACKED_WORD 64
#define CHAR_MASK 3ULL

#endif
