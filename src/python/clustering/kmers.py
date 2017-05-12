from itertools import product, tee, izip
from collections import Counter

import numpy as np

# ----------------------------------------------------------------------------

def kmer_profile(seq, k, feature_mapping, nr_features):
  seq_len = len(seq)
  
  # create a list containing all kmers, translated to integers
  kmers = [
          feature_mapping[kmer_tuple]
          for kmer_tuple 
          in _window(str(seq).upper(), k)
          if kmer_tuple in feature_mapping
          ]

  # numpy.bincount returns an array of size = max + 1
  # so we add the max value and remove it afterwards
  # numpy.bincount was found to be much more efficient than
  # counting manually or using collections.Counter
  kmers.append(nr_features - 1)
  composition_v = np.bincount(np.array(kmers))
  composition_v[-1] -= 1

  return composition_v

def init_feature_mapping(kmer_len):
  BASE_COMPLEMENT = {"A":"T","T":"A","G":"C","C":"G"}
  kmer_hash = {}
  counter = 0
  for kmer in product("ATGC",repeat=kmer_len):
    if kmer not in kmer_hash:
      kmer_hash[kmer] = counter
      rev_compl = tuple([BASE_COMPLEMENT[x] for x in reversed(kmer)])
      kmer_hash[rev_compl] = counter
      counter += 1

  return kmer_hash, counter

# ----------------------------------------------------------------------------
# helpers

def _window(seq,n):
  els = tee(seq,n)
  for i,el in enumerate(els):
    for _ in xrange(i):
      next(el, None)

  return izip(*els)