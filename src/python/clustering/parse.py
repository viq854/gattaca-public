import logging
import numpy as np
import pysam
from collections import OrderedDict

from kmers import kmer_profile, init_feature_mapping

# ----------------------------------------------------------------------------
# load data

def parse_table(coverage_table_path, min_len=1000):
  """Loads coverage table into a dictionary

  Takes as argument paths; outputs contig names, lengths, and coverage matrix
  """

  table = np.genfromtxt(coverage_table_path, dtype=None, delimiter='\t', 
                        skip_header=1)

  names = [row[0] for row in table]
  lengths = [row[1] for row in table]

  long_ctgs = [(i, n, l) for i, (n,l) in enumerate(zip(names, lengths)) if l >= min_len]
  long_i, long_names, long_lengths = zip(*long_ctgs)
  long_covs = np.array([list(table[i])[2:] for i in long_i])

  # TODO: check if the matrix is okay in terms of data types

  # for debugging:
  # long_names = [row[0] for row in table]
  # long_lengths = [row[1] for row in table]
  # long_covs = np.zeros((len(long_names),1))

  return long_names, long_lengths, long_covs

def parse_kmers(contigs_path, k=4, min_len=1000):
  """Extracts kmer profiles from contig FASTA

  Returns dict of the form kmer[ctg] = profile
  """

  kmers = OrderedDict()
  contigs = pysam.FastaFile(contigs_path)
  feature_mapping, feature_num = init_feature_mapping(k)

  for i, ctg in enumerate(contigs.references):
    if i % 5000 == 0:
      logging.info('%d/%d contigs parsed' % (i, len(contigs)))
    seq = contigs.fetch(ctg)
    if len(seq) < min_len:
      continue
    kmers[ctg] = kmer_profile(seq, k, feature_mapping, feature_num)

  return kmers

def parse_input_table(table_path):
  """Parse pre-processed CONCOCT input table (for debugging)"""
  names = list()
  data = list()
  with open('./khmer-full.csv') as f:
    f.readline()
    for i, line in enumerate(f):
      fields = line.split(',')
      names.append(fields[0])
      data.append(np.array([float(x) for x in fields[1:]]))

  return names, np.array(data).T

# ----------------------------------------------------------------------------
# save data

def save_clusters(names, clusters, cluster_file):
  with open(cluster_file, 'w') as out:
    for name, cluster in zip(names, clusters):
      out.write('%s,%d\n' % (name, int(cluster)))

