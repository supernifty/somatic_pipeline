#!/usr/bin/env python
'''
  adds AF and DP to strelka output
  - assume that this is a single sample
'''

import argparse
import logging
import sys

import cyvcf2

def main(sample_name):
  logging.info('reading vcf from stdin...')
  vcf_in = cyvcf2.VCF('-')

  vcf_in.add_info_to_header({'ID': 'AF', 'Description': 'Calculated allele frequency', 'Type':'Float', 'Number': '1'})
  vcf_in.add_info_to_header({'ID': 'DP', 'Description': 'Calculated depth', 'Type':'Float', 'Number': '1'})
  sample_id = vcf_in.samples.index(sample_name)

  sys.stdout.write(vcf_in.raw_header)

  stats = { 'min_af': 1e6, 'max_af': -1, 'min_dp': 1e6, 'max_dp': -1, 'allowed': 0, 'denied': 0}
  for variant in vcf_in:
    ok = True

    # pindel
    if all([x in variant.FORMAT for x in ('PP', 'NP', 'PR', 'NR')]):
      dp = variant.format('PR')[0][0] + variant.format('NR')[0][0]
      af = (variant.format('PP')[0][0] + variant.format('NP')[0][0]) / dp
      variant.INFO["AF"] = af
      variant.INFO["DP"] = float(dp)

    # strelka
    elif all([x in variant.FORMAT for x in ('TAR', 'TIR')]):
      tir = sum(variant.format('TIR')[sample_id])
      tar = sum(variant.format('TAR')[sample_id])
      if tir + tar == 0:
        logging.warn('No reads found for TIR or TAR')
        af = 0.0
      else:
        af = tir / (tir + tar)

      dp = int(sum([sum(x) for x in variant.format('DP')]))
      variant.INFO["AF"] = af
      variant.INFO["DP"] = dp

    else: # not strelka or pindel
      logging.warn('Failed to add AF to variant')
      af = 0.0
      dp = 0.0

    stats['min_af'] = min(stats['min_af'], af)
    stats['max_af'] = max(stats['max_af'], af)
    stats['min_dp'] = min(stats['min_dp'], dp)
    stats['max_dp'] = max(stats['max_dp'], dp)

    if ok:
      sys.stdout.write(str(variant))
      stats['allowed'] += 1
    else:
      stats['denied'] += 1

    if (stats['allowed']) % 100000 == 0:
      logging.debug('stats: %s.', ', '.join(['{}: {}'.format(x, stats[x]) for x in stats]))
    
  logging.info('done. stats: %s.', ', '.join(['{}: {}'.format(x, stats[x]) for x in stats]))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Extract samples from VCF')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--sample', required=True, help='sample name')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
  main(args.sample)


