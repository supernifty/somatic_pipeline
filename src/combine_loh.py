#!/usr/bin/env python
'''
  given tumour and normal vcf pairs, explore msi status
'''

import argparse
import csv
import gzip
import logging
import sys

import intervaltree

def main(lohs, transcripts, min_accept):
  logging.info('processing %s', transcripts)
  chroms = {}
  rows = 0
  for rows, row in enumerate(csv.DictReader(gzip.open(transcripts, 'rt'), delimiter='\t')):
    chrom = row['chrom'].replace('chr', '')
    if chrom not in chroms:
      chroms[chrom] = intervaltree.IntervalTree()
      logging.info('added %s', chrom)
    chroms[chrom][int(row['txStart']):int(row['txEnd'])] = row['name2']
  logging.info('%i records processed', rows)

  written = set()
  sys.stdout.write('sample\tgene\taccept\n')
  for loh in lohs:
    logging.info('processing %s', loh)
    sample = loh.split('/')[-1].split('.')[0]
    for line in open(loh, 'r'):
      #12      6709614 6711147 50.0    1       0       1       1533
      fields = line.strip('\n').split('\t')
      chrom, start, finish, accept = fields[0], int(fields[1]), int(fields[2]), int(fields[4])
      if accept < min_accept:
        continue
      chrom = chrom.replace('chr', '')
      if chrom in chroms:
        overlaps = chroms[chrom][start:finish]
        for overlap in overlaps:
          if (sample, overlap.data) not in written:
            sys.stdout.write('{}\t{}\t{}\n'.format(sample, overlap.data, accept))
            written.add((sample, overlap.data))
      else:
        logging.warn('chrom %s in %s not found in %s', chrom, loh, transcripts)

  logging.info('done. wrote %i records for %i samples', len(written), len(lohs))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Combine LOH')
  parser.add_argument('--lohs', required=True, nargs='+', help='loh files')
  parser.add_argument('--transcripts', required=True, help='transcripts')
  parser.add_argument('--min_accept', required=False, default=1, type=int, help='transcripts')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.lohs, args.transcripts, args.min_accept)
