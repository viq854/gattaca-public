import logging
import numpy as np

from sklearn.decomposition import PCA

# ----------------------------------------------------------------------------

def create_matrix(lengths, covs, kmers):
  """Create coverage matrix from coverages and kmers"""

  C = _create_and_normalize_covs(covs, lengths)
  K = _create_and_normalize_kmers(kmers)

  return np.hstack((K,C))

def reduce_dimensionality(X, target_var=0.9):
  """Project to lower dimensional subspace using PCA"""
  
  pca = PCA(n_components=target_var).fit(X)
  return pca.transform(X)

# ----------------------------------------------------------------------------
# helpers

def _create_and_normalize_covs(covs, lengths, rlen=100):
  """Normalize coverage matrix"""

  X = covs.astype(np.dtype('float64'))
  n, f = X.shape

  # add pseudocounts
  X += float(rlen) / np.array(lengths).reshape((n,1))

  # normalize read count across samples (contigs/rows sum to one)
  # TODO: try removing this, or replace by whitening transform
  X /= np.sum(X, axis=0).reshape((1,f))
  tot_cov = np.sum(X, axis=1)

  # normalize read count across contigs (samples/columns sum to one)
  X /= np.sum(X, axis=1).reshape((n,1))

  # add total coverage as a feature
  X = np.hstack((X,tot_cov.reshape(n,1)))

  # take logs
  X = np.log(X)

  return X

def _create_and_normalize_kmers(kmers):
  """Create and normalize kmer matrix"""

  # create matrix
  names = kmers.keys()
  n, f = len(names), len(kmers[names[0]])
  X = np.empty((n,f))

  for i, ctg in enumerate(names):
    X[i] = kmers[ctg].astype(np.dtype('float64'))

  # add pseudocounts to matrix
  X += 1.

  # normalize matrix
  X /= np.sum(X, axis=1).reshape((n,1))

  # take logs
  X = np.log(X)

  return X
