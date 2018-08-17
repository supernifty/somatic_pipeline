#!/usr/bin/env python
'''
  write out details of a variant if found
'''

import argparse
import logging
import sys

import cyvcf2

def main(sample, chrom, pos, nofilter):
  logging.info('reading from stdin...')

  vcf_in = cyvcf2.VCF('-')  
  sample_id = vcf_in.samples.index(sample)

  for variant in vcf_in:
    if not nofilter and variant.FILTER is not None:
      continue

    if variant.POS == pos and variant.CHROM == chrom:
      # check gt 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT
      gt = variant.gt_types[sample_id]
      if gt == 1 or gt == 3:
        ad = variant.format('AD')[sample_id]
        gt_str = ['0/0', '0/1', './.', '1/1'][gt]
        sys.stdout.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(sample, chrom, pos, ', '.join([str(x) for x in ad]), gt_str, '1'))
        logging.info('done')
        sys.exit(0)

  # not found
  sys.stdout.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(sample, chrom, pos, 'NA', 'NA', '0'))
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Max coverage')
  parser.add_argument('--sample', required=True, help='sample name')
  parser.add_argument('--chrom', required=True, help='location of variant')
  parser.add_argument('--pos', required=True, type=int, help='location of variant')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--nofilter', action='store_true', help='allow non-pass')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.sample, args.chrom, args.pos, args.nofilter)
