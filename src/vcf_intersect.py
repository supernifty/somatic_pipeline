#!/usr/bin/env python
'''
  intersect vcf files
'''

import argparse
import logging
import sys

import numpy

import cyvcf2

def main(vcfs, rejected, pass_only):
  '''
  '''
  vcf_in = cyvcf2.VCF(vcfs[0])  
  variant_base_count = 0
  variant_base_pass = 0
  base = {}
  logging.info('reading %s...', vcfs[0])
  for variant_base_count, variant in enumerate(vcf_in):
    if (variant_base_count + 1 ) % 100000 == 0:
      logging.info('reading %s: %i variants processed...', vcfs[0], variant_base_count + 1)
    if pass_only and variant.FILTER is not None:
      continue
    variant_base_pass += 1
    if variant.CHROM not in base:
      base[variant.CHROM] = set()
    base[variant.CHROM].add(variant.POS)

  logging.info('done reading %s: %i variants processed...', vcfs[0], variant_base_count + 1)

  logging.info('reading %s...', vcfs[1])
  vcf_cand = cyvcf2.VCF(vcfs[1])  
  variant_cand_count = 0
  variant_cand_pass = 0
  included = 0
  reject = 0
  sys.stdout.write(vcf_cand.raw_header)

  if rejected is not None:
    rejected_fh = open(rejected, 'w')
    rejected_fh.write(vcf_cand.raw_header)

  for variant_cand_count, variant in enumerate(vcf_cand):
    if (variant_cand_count + 1 ) % 100000 == 0:
      logging.info('reading %s: processed %i variants. wrote %i...', vcfs[1], variant_cand_count + 1, included)
    if pass_only and variant.FILTER is not None:
      continue
    variant_cand_pass += 1
    if variant.CHROM in base and variant.POS in base[variant.CHROM]:
      sys.stdout.write(str(variant))
      included += 1
    else:
      reject += 1
      if rejected is not None:
        rejected_fh.write(str(variant))

  logging.info('done. %s: %i passed variants. %s: %i passed variants. wrote %i variants. rejected %i variants', vcfs[1], variant_cand_pass + 1, vcfs[0], variant_base_pass + 1, included, reject)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Intersect vcfs')
  parser.add_argument('--inputs', required=True, nargs='+', help='input vcf files')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--pass_only', action='store_true', help='remove non-pass')
  parser.add_argument('--rejected', required=False, help='file to write rejected to')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
  main(args.inputs, args.rejected, args.pass_only)

