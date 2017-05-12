import logging
import pickle
import numpy as np
from sklearn.mixture import DPGMM, GMM

from preprocess import create_matrix, reduce_dimensionality
from algorithms.dirichlet import dirichlet_em
from algorithms.ard import variational_em
from algorithms.agglomerative import agglomerative

# ----------------------------------------------------------------------------

def clustering_algorithm(lengths, covs, kmers, algorithm='dirichlet', K=300,
                         max_epoch=25, t=0, seed=None, mu_pkl=None):
  """Clusters using given algorithm

  Takes as argument cluster names, lengths, and coverage/kmer matrices.
  """

  # create matrix for clustering
  logging.info('Creating data matrix')
  X = create_matrix(lengths, covs, kmers)

  # project down dimension
  logging.info('Performing dimensionality reduction')
  X = reduce_dimensionality(X)

  # do the clustering
  logging.info('Starting clustering algorithm')
  if algorithm == 'sk-gmm':
    gmm = GMM(n_components=K, covariance_type='full', n_iter=500)
    gmm.fit(X)
    z = gmm.predict(X)
    return z
  elif algorithm == 'sk-dpgmm':
    gmm = DPGMM(n_components=K, covariance_type='full', n_iter=500)
    gmm.fit(X)
    z = gmm.predict(X)
    return z
  elif algorithm == 'dirichlet':
    n_data = X.shape[0]
    mu_pred_dem, Sigma_pred_dem, asgn_dem, llik = dirichlet_em(X.T, K=K, 
      n_minibatch=n_data, max_epoch=max_epoch, seed=seed)

    if mu_pkl:
      with open(mu_pkl, 'wb') as f:
        pickle.dump((mu_pred_dem, asgn_dem), f)

    # further compress stuff 
    compressed_clusters = agglomerative(mu_pred_dem.T, t=t)
    transl_dict = { i : c for i, c in enumerate(compressed_clusters)}
    asign_dem_agg = np.array([transl_dict[i] for i in asgn_dem])

    return asign_dem_agg
  elif algorithm == 'ard':
    n_data = X.shape[0]
    _, _, asgn, _ = variational_em(X.T, K=K, n_minibatch=n_data, max_epoch=max_epoch)
    return asgn
  else:
    raise ValueError("Invalid algorithm name")
