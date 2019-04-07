#!/usr/bin/env python
'''
  calculate and plot af from ad
'''

import argparse
import logging
import sys

import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import cyvcf2

def main(sample, dp_threshold, target, info_af, log):

  logging.info('reading from stdin...')

  vcf_in = cyvcf2.VCF("-")  
  ads = []

  variant_count = 0
  skipped_pass = skipped_dp = skipped_af = allowed = 0
  sample_id = vcf_in.samples.index(sample)
  for variant_count, variant in enumerate(vcf_in):
    # GL000220.1      135366  .       T       C       .       LowEVS;LowDepth SOMATIC;QSS=1;TQSS=1;NT=ref;QSS_NT=1;TQSS_NT=1;SGT=TT->TT;DP=2;MQ=60.00;MQ0=0;ReadPosRankSum=0.00;SNVSB=0.00;SomaticEVS=0.71    DP:FDP:SDP:SUBDP:AU:CU:GU:TU    1:0:0:0:0,0:0,0:0,0:1,1 1:0:0:0:0,0:1,1:0,0:0,0
    if (variant_count + 1 ) % 100000 == 0:
      logging.info('%i variants processed...', variant_count + 1)

    if len(variant.ALT) > 1:
      logging.warn('variant %i is multi-allelic', variant_count + 1)

    if variant.FILTER is not None and variant.FILTER != 'alleleBias': # PASS only, or alleleBias for platypus
      skipped_pass += 1
      continue

    if variant.INFO["DP"] < dp_threshold: # somatic + germline
      skipped_dp += 1
      continue

    if info_af:
      ads.append(variant.INFO["AF"])
    else:
      ad = variant.format("AD")[sample_id]
      ref = ad[0]
      alt = ad[1]
      if ref + alt > 0:
        ads.append(alt / (ref + alt))
      else:
        ads.append(0)
    allowed += 1

  logging.info('processed %i variants. no pass %i. low af %i. low dp %i. allowed %i AF range %.2f to %.2f', variant_count + 1, skipped_pass, skipped_af, skipped_af, allowed, min(ads), max(ads))

  # now plot histogram
  plt.hist(ads, bins=int(max(ads) * 100))
  if log:
    plt.yscale('log')
  plt.grid(which='both')
  plt.savefig(target)
  
  

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Filter VCF')
  parser.add_argument('--sample', required=True,  help='sample name')
  parser.add_argument('--target', required=True,  help='image output')
  parser.add_argument('--dp', type=int, required=False, default=0, help='minimum dp')
  parser.add_argument('--info_af', action='store_true', help='info af')
  parser.add_argument('--log', action='store_true', help='log on y scale')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  # sample af vcf
  main(args.sample, args.dp, args.target, args.info_af, args.log)
