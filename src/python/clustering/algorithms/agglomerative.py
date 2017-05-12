from scipy.cluster.hierarchy import fcluster, linkage
import numpy as np

# ----------------------------------------------------------------------------

def agglomerative(X, t=0, linkage_type='ward'):
  Z = linkage(X, linkage_type)
  clusters = fcluster(Z, t, criterion='distance')
  return clusters
