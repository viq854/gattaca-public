#!/usr/bin/env python
import argparse
import logging

from clustering import perform_clustering

# ----------------------------------------------------------------------------

def main():
  parser = argparse.ArgumentParser()
  subparsers = parser.add_subparsers(title='Commands')

  # indexing arguments

  index_parser = subparsers.add_parser('index')
  index_parser.set_defaults(func=index)

  index_parser.add_argument('sequences', 
    help='List of files to index')
  index_parser.add_argument('-k', type=int, default=16,
    help='k-mer size used for indexing')
  index_parser.add_argument('--log',
    help='Logfile path (default is stdout)')

  # clustering arguments

  cluster_parser = subparsers.add_parser('cluster')
  cluster_parser.set_defaults(func=cluster)

  cluster_parser.add_argument('--input_table',
    help='Input table from Concoct (for testing)')
  cluster_parser.add_argument('--contigs',
    help='Contigs to cluster in FASTA format')
  cluster_parser.add_argument('--coverage',
    help='Coverage table for contigs')
  cluster_parser.add_argument('--clusters',
    help='Output file')
  cluster_parser.add_argument('--K', type=int, default=300,
    help='Maximum number of clusters to consider')
  cluster_parser.add_argument('--t', type=int, default=0,
    help='Post-processing mergin threshold')
  cluster_parser.add_argument('--iter', type=int, default=25,
    help='Number of iterations of EM')
  cluster_parser.add_argument('--cluster-means',
    help='Load precomputed cluster means (for internal dev purposes)')
  cluster_parser.add_argument('--seed', type=int,
    help='Random seed for initializing kmeans')
  cluster_parser.add_argument('--algorithm', action='store',
                              choices=['sk-gmm', 'sk-dpgmm', 
                                       'ard', 'dirichlet', 'preload'], 
                              default='dirichlet',
                              help='Clustering algorithm')
  cluster_parser.add_argument('--log', 
    help='Logfile path (default is stdout)')
  
  # parse arguments
  args = parser.parse_args()

  # initialize logging module
  logging.basicConfig(filename=args.log, filemode='a', level=logging.INFO)

  # run handler function
  args.func(args)

# ----------------------------------------------------------------------------
# entry point functions

def index(args):
  """Launches the indexing of a set of reference metagenomes.

  Takes as argument a textfile with one fasta/fastq file per line.
  """

  logging.info('Starting indexing sequences in %s' % args.sequences)
  logging.error('TODO: Implement indexing!')

def cluster(args):
  """Launches the clustering algorithm

  Takes either contigs together with an index (in which case it first computes
  coverages within reference contigs) or a coverage table in CONCOCT format.
  """

  # if not (args.coverage or args.index):
  #   logging.error('Must specify a coverage file or contigs + reference index.')

  logging.info('Starting clustering process')
  perform_clustering(args)

if __name__ == '__main__':
  main()