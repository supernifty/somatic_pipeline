#!/usr/bin/env python
'''
  very simple filter
'''

import argparse
import logging
import sys

import numpy

import cyvcf2

def main(sample, af_threshold, dp_threshold, info_af, pass_only, dp_field, min_sample_ad):

  logging.info('reading from stdin...')

  vcf_in = cyvcf2.VCF("-")  
  vcf_in.add_info_to_header({'ID': 'AF', 'Description': 'Calculated allele frequency', 'Type':'Float', 'Number': '1'})

  sys.stdout.write(vcf_in.raw_header)

  variant_count = 0
  skipped_pass = skipped_dp = skipped_af = skipped_ad = allowed = 0
  sample_id = vcf_in.samples.index(sample)
  for variant_count, variant in enumerate(vcf_in):
    # GL000220.1      135366  .       T       C       .       LowEVS;LowDepth SOMATIC;QSS=1;TQSS=1;NT=ref;QSS_NT=1;TQSS_NT=1;SGT=TT->TT;DP=2;MQ=60.00;MQ0=0;ReadPosRankSum=0.00;SNVSB=0.00;SomaticEVS=0.71    DP:FDP:SDP:SUBDP:AU:CU:GU:TU    1:0:0:0:0,0:0,0:0,0:1,1 1:0:0:0:0,0:1,1:0,0:0,0
    if (variant_count + 1 ) % 100000 == 0:
      logging.info('%i variants processed...', variant_count + 1)

    if len(variant.ALT) > 1:
      logging.warn('variant %i is multi-allelic', variant_count + 1)

    if pass_only and variant.FILTER is not None and variant.FILTER != 'alleleBias': # PASS only, or alleleBias for platypus
      skipped_pass += 1
      continue

    if dp_threshold > 0:
      try:
        if variant.INFO[dp_field] < dp_threshold: # somatic + germline
          logging.debug('variant at %i skipped with info depth %i', variant.POS, variant.INFO[dp_field])
          skipped_dp += 1
          continue
      except:
        if variant.format(dp_field)[sample_id][0] < dp_threshold: # somatic
          logging.debug('variant at %i skipped with sample depth %i', variant.POS, variant.format(dp_field)[sample_id][0])
          skipped_dp += 1
          continue

    logging.debug('variant at %i has OK depth', variant.POS)

    if info_af:
      try:
        af = variant.INFO["AF"]
      except:
        skipped_af += 1
        logging.debug('variant at %i has no info AF', variant.POS)
        continue
    else:
      try:
        af = variant.format("AF")[sample_id][0] 
      except:
        skipped_af += 1
        logging.debug('variant at %i has no sample AF', variant.POS)
        continue

    if af < af_threshold:
      skipped_af += 1
      continue

    if min_sample_ad is not None:
      low_ad = False
      for msad in min_sample_ad:
        soi, min_ad = msad.split('=')
        min_ad = int(min_ad)
        soi_idx = vcf_in.samples.index(soi)
        if sum(variant.format("AD")[soi_idx]) < min_ad:
          skipped_ad += 1
          logging.debug('variant at %i has low AD', variant.POS)
          low_ad = True
          break
      if low_ad:
        continue

    logging.debug('variant at %i is OK', variant.POS)
    allowed += 1
    sys.stdout.write(str(variant))

  logging.info('processed %i variants. no pass %i. low af %i. low dp %i. low sample ad %i. allowed %i', variant_count + 1, skipped_pass, skipped_af, skipped_dp, skipped_ad, allowed)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Filter VCF')
  parser.add_argument('--sample', required=True,  help='sample name for af filter')
  parser.add_argument('--af', type=float, required=True,  help='minimum af')
  parser.add_argument('--dp', type=int, required=True,  help='minimum dp')
  parser.add_argument('--dp_field', required=False, default='DP',  help='name of DP field')
  parser.add_argument('--info_af', action='store_true', help='use af from info')
  parser.add_argument('--pass_only', action='store_true', help='reject non-pass')
  parser.add_argument('--min_sample_ad', nargs='*', required=False,  help='minimum ad sum for sample in the form sample=ad...')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  # sample af vcf
  main(args.sample, args.af, args.dp, args.info_af, args.pass_only, args.dp_field, args.min_sample_ad)
