import logging
import pickle
import numpy as np

from parse import parse_table, parse_kmers, save_clusters, parse_input_table
from cluster import clustering_algorithm
from algorithms.dirichlet import dirichlet_em
from algorithms.agglomerative import agglomerative

# ----------------------------------------------------------------------------

def perform_clustering(args):
  """Starts clustering process.

  Takes as argument the arg dict from the cluster subparser in main script.
  """

  # for debugging
  if args.algorithm == 'preload':
    with open(args.cluster_means, 'rb') as f:
      (mu_pred_dem, asgn_dem) = pickle.load(f)

    # further compress stuff 
    compressed_clusters = agglomerative(mu_pred_dem.T, t=args.t)
    transl_dict = { i : c for i, c in enumerate(compressed_clusters)}
    asign_dem_agg = np.array([transl_dict[i] for i in asgn_dem])

    names, lengths, coverages = parse_table(args.coverage)
    save_clusters(names, asign_dem_agg, args.clusters)
    return

  if args.input_table:
    # debug mode
    names, X = parse_input_table(args.input_table)
    mu_pred_dem, Sigma_pred_dem, asgn_dem, llik \
      = dirichlet_em(X, K=300, n_minibatch=40000, max_epoch=args.iter)

    save_clusters(names, asign_dem_agg, args.clusters)
  else:
    # the real thing
    logging.info("Loading coverage data")

    # # load coverage data
    if args.coverage:
      names, lengths, coverages = parse_table(args.coverage)
    elif args.contigs and args.index:
      logging.error("Handling contig/indices is not yet implemented")
    else:
      raise ValueError("Invalid arguments")

    logging.info("Loading composition data")

    # # load kmer data
    kmers = parse_kmers(args.contigs)
    
    logging.info("Starting clustering using algorithm '%s'" % args.algorithm)

    clusters = clustering_algorithm(lengths, coverages, kmers, args.algorithm, 
                                    args.K, args.iter, args.t, args.seed, 
                                    args.cluster_means)

    logging.info("Saving clusters in %s" % args.clusters)

    save_clusters(names, clusters, args.clusters)

    logging.info("wohoo!")
