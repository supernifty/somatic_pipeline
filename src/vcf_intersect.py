#!/usr/bin/env python
'''
  intersect vcf files
'''

import argparse
import logging
import sys

import numpy

import cyvcf2

def is_pass(variant, allowed_filters):
  if variant.FILTER is None:
    return True
  elif allowed_filters is None:
    return variant.FILTER is None
  else:
    for f in variant.FILTER:
      if f not in allowed_filters:
        return False
    return True

def main(vcfs, rejected, pass_only, pass_one, allowed_filters):
  '''
  '''
  vcf_in = cyvcf2.VCF(vcfs[0])  
  variant_base_count = 0
  variant_base_pass = 0
  base = {}
  already_passed = {}
  logging.info('reading %s...', vcfs[0])
  for variant_base_count, variant in enumerate(vcf_in):
    if (variant_base_count + 1 ) % 100000 == 0:
      logging.info('reading %s: %i variants processed...', vcfs[0], variant_base_count + 1)
    if pass_only and not is_pass(variant, allowed_filters):
      continue
    # variant is a pass or pass_only is false
    variant_base_pass += 1
    if variant.CHROM not in base:
      base[variant.CHROM] = set()
      already_passed[variant.CHROM] = set()

    if pass_one and is_pass(variant, allowed_filters):
      already_passed[variant.CHROM].add(variant.POS) # wait for others
    else:
      base[variant.CHROM].add(variant.POS) # wait and see the others

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
    if pass_only and not is_pass(variant, allowed_filters): # pass is required
      continue
    # it's a pass or pass_only is false
    variant_cand_pass += 1
    if variant.CHROM in base and variant.POS in base[variant.CHROM]: # seen previously
      if pass_one:
        if is_pass(variant, allowed_filters): # 2nd is a pass
          sys.stdout.write(str(variant))
          included += 1
        elif variant.CHROM in already_passed and variant.POS in already_passed[variant.CHROM]:
          sys.stdout.write(str(variant))
          included += 1
        else:
          reject += 1
          if rejected is not None:
            rejected_fh.write(str(variant))
      else: # print regardless
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
  parser.add_argument('--pass_one', action='store_true', help='just one pass is required')
  parser.add_argument('--allowed_filters', required=False, nargs='*', help='input vcf files')
  parser.add_argument('--rejected', required=False, help='file to write rejected to')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
  main(args.inputs, args.rejected, args.pass_only, args.pass_one, args.allowed_filters)

