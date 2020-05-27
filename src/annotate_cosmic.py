#!/usr/bin/env python
'''
  given tumour and normal vcf pairs, explore msi status
'''

import argparse
import csv
import logging
import sys

import cyvcf2

def main(cosmic):
  logging.info('reading cosmic file...')
  #cds = collections.defaultdict(int)
  #aa = collections.defaultdict(int)
  #g = collections.defaultdict(int)
  #for row_count, row in enumerate(csv.reader(cosmic, delimiter='\t')):
  #  if row_count % 10000 == 0:
  #    logging.info('%i records read...', row_count)
  #
  #    cds['{}:{}'.format(row['Gene name'], row['Mutation CDS'])] += 1
  #    aa['{}:{}'.format(row['Gene name'], row['Mutation AA'])] += 1
  #    if '-' in row['Mutation genome position']:
  #      g[row['Mutation genome position'].split('-')[0]] += 1


  counts = {}
  total = 0
  for total, variant in enumerate(cyvcf2.VCF(cosmic)):
    position = '{}:{} {}/{}'.format(variant.CHROM, variant.POS, variant.REF, variant.ALT[0])
    counts[position] = variant.INFO['CNT'] # just overwrite repeats
    if total % 100000 == 0:
      logging.debug('read %i lines from COSMIC, last was %s: %i', total, position, counts[position])

  logging.info('annotating vcf...')

  vcf_in = cyvcf2.VCF('-')
  vcf_in.add_info_to_header({'ID': 'cosmic', 'Description': 'Number of times position seen in COSMIC', 'Type':'Character', 'Number': '1'})
  sys.stdout.write(vcf_in.raw_header)
  total = seen = 0
  summary = {'max': 0, 'sum': 0, 'maxpos': None}
  for total, variant in enumerate(vcf_in):
    position = '{}:{} {}/{}'.format(variant.CHROM, variant.POS, variant.REF, variant.ALT[0])
    if position in counts:
      variant.INFO['cosmic'] = counts[position]
      seen += 1
      if counts[position] > summary['max']:
        summary['max'] = counts[position]
        summary['maxpos'] = position
      summary['sum'] += counts[position]
    else:
      variant.INFO['cosmic'] = 0 # not seen
    if total % 1000 == 0:
      logging.debug('read %i lines from vcf, last was %s', total, position)

    sys.stdout.write(str(variant))

  logging.info('done updating %i records. saw %i cosmic variants. max count %i at %s. total count %i', total, seen, summary['max'], summary['maxpos'], summary['sum'])

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Annotate VCF with COSMIC data')
  parser.add_argument('--cosmic', required=True, help='cosmic file')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.cosmic)


