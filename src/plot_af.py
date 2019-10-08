#!/usr/bin/env python
'''
  calculate and plot af from ad
'''

import argparse
import collections
import logging
import random
import sys

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import rcParams
rcParams['figure.figsize'] = 16, 10

import cyvcf2

colors = {"SBS1": "#de3860", "SBS2": "#41ac2f", "SBS3": "#7951d0", "SBS4": "#73d053", "SBS5": "#b969e9", "SBS6": "#91ba2c", "SBS7a": "#b4b42f", "SBS7b": "#5276ec", "SBS7c": "#daae36", "SBS7d": "#9e40b5", "SBS8": "#43c673", "SBS9": "#dd4cb0", "SBS10a": "#3d9332", "SBS10b": "#de77dd", "SBS11": "#7bad47", "SBS12": "#9479e8", "SBS13": "#487b21", "SBS14": "#a83292", "SBS15": "#83c67d", "SBS16": "#664db1", "SBS17a": "#e18d28", "SBS17b": "#588de5", "SBS18": "#e2672a", "SBS19": "#34c7dd", "SBS20": "#cf402b", "SBS21": "#5acdaf", "SBS22": "#d74587", "SBS23": "#428647", "SBS24": "#7b51a7", "SBS25": "#b4ba64", "SBS26": "#646cc1", "SBS27": "#a27f1f", "SBS28": "#3b63ac", "SBS29": "#dca653", "SBS30": "#505099", "SBS31": "#7d8529", "SBS32": "#bf8ade", "SBS33": "#516615", "SBS34": "#b65da7", "SBS35": "#57a87a", "SBS36": "#c84249", "SBS37": "#37b5b1", "SBS38": "#a14622", "SBS39": "#58b5e1", "SBS40": "#ba6e2f", "SBS41": "#589ed8", "SBS42": "#e98261", "SBS43": "#3176ae", "SBS44": "#656413", "SBS45": "#a19fe2", "SBS46": "#756121", "SBS47": "#7e4a8d", "SBS48": "#326a38", "SBS49": "#dd8abf", "SBS50": "#1a6447", "SBS51": "#e78492", "SBS52": "#30876c", "SBS53": "#9d4d7c", "SBS54": "#919d5b", "SBS55": "#9d70ac", "SBS56": "#5b6f34", "SBS57": "#65659c", "SBS58": "#c9a865", "SBS59": "#a1455d", "SBS60": "#5e622c", "SBS84": "#b66057", "SBS85": "#dca173", "DBS1": "#855524", "DBS2": "#9f7846", "DBS3": "#7951d0", "DBS4": "#73d053", "DBS5": "#b969e9", "DBS6": "#91ba2c", "DBS7": "#3656ca", "DBS8": "#b4b42f", "DBS9": "#5276ec", "DBS10": "#daae36", "DBS11": "#9e40b5", "ID1": "#de3860", "ID2": "#41ac2f", "ID3": "#7951d0", "ID4": "#73d053", "ID5": "#b969e9", "ID6": "#91ba2c", "ID7": "#9e40b5", "ID8": "#43c673", "ID9": "#dd4cb0", "ID10": "#3d9332", "ID11": "#de77dd", "ID12": "#7bad47", "ID13": "#9479e8", "ID14": "#487b21", "ID15": "#a83292", "ID16": "#83c67d", "ID17": "#664db1"}

def main(sample, dp_threshold, target, info_af, log, filter, just_pass, use_likelihoods, percent, title):

  logging.info('reading from stdin...')

  vcf_in = cyvcf2.VCF("-")  
  ads = []
  ads_nopass = []
  sig_ads = collections.defaultdict(list)

  if filter is not None:
    if ':' in filter:
      chromosome, rest = filter.split(':')
      start, finish = [int(x) for x in rest.split('-')]
    else:
      chromosome = filter.split(':')[0]
      start = finish = None
  else:
    chromosome = None

  variant_count = 0
  skipped_pass = skipped_dp = skipped_af = allowed = 0
  sample_id = vcf_in.samples.index(sample)
  for variant_count, variant in enumerate(vcf_in):
    # GL000220.1      135366  .       T       C       .       LowEVS;LowDepth SOMATIC;QSS=1;TQSS=1;NT=ref;QSS_NT=1;TQSS_NT=1;SGT=TT->TT;DP=2;MQ=60.00;MQ0=0;ReadPosRankSum=0.00;SNVSB=0.00;SomaticEVS=0.71    DP:FDP:SDP:SUBDP:AU:CU:GU:TU    1:0:0:0:0,0:0,0:0,0:1,1 1:0:0:0:0,0:1,1:0,0:0,0
    if (variant_count + 1 ) % 100000 == 0:
      logging.info('%i variants processed...', variant_count + 1)

    if len(variant.ALT) > 1:
      logging.warn('variant %i is multi-allelic', variant_count + 1)

    if chromosome is not None:
      if chromosome != variant.CHROM:
        continue
      if start is not None and not (start <= variant.POS < finish):
        continue

    is_pass = variant.FILTER is None or variant.FILTER == 'alleleBias'
    if is_pass:
      target_ad = ads
    else:
      target_ad = ads_nopass
      skipped_pass += 1 # but we'll still record

    if variant.INFO["DP"] < dp_threshold: # somatic + germline
      skipped_dp += 1
      continue

    if info_af:
      value = variant.INFO["AF"]
    else:
      ad = variant.format("AD")[sample_id]
      ref = ad[0]
      alt = ad[1]
      if ref + alt > 0:
        value = alt / (ref + alt)
      else:
        value = 0
    target_ad.append(value) # overall set of afs
    allowed += 1
    if use_likelihoods:
      try:
        likelihoods = variant.INFO["signature_likelihood"]
      except:
        logging.warn('Signature likelihood not found. Is the VCF annotated?')
        continue
      items = [sigvalue.split('/') for sigvalue in likelihoods.split(',')]
      guess = random.random()
      start = 0.0
      sig = None
      for item in items:
        if start < guess < start + float(item[1]):
          sig = item[0]
          break
        start += float(item[1])
      if is_pass or not just_pass:
        if sig is not None:
          sig_ads[sig].append(value)

  min_ads = min(ads + ads_nopass)
  max_ads = max(ads + ads_nopass)
  logging.info('processed %i variants. no pass %i. low af %i. low dp %i. allowed %i AF range %.2f to %.2f', variant_count + 1, skipped_pass, skipped_af, skipped_af, allowed, min_ads, max_ads)

  # now plot histogram
  if use_likelihoods:
    sig_names = sorted(sig_ads.keys())
    if just_pass:
      xmax = max(ads + [0.01])
    else:
      xmax = max(ads + ads_nopass + [0.01])
    if len(sig_ads) == 0: # no variants do a dummy
      plt.hist([0.0])
    else:
      if percent:
        plt.ylabel('Proportion of variants')
        ax = plt.axes()
        totals = [0] * int(xmax * 100)
        sig_totals = {}
        values = [sig_ads[n] for n in sig_names]
        for value_idx, value in enumerate(values):
          sig_totals[sig_names[value_idx]] = [0] * int(xmax * 100)
          for bin in range(0, int(xmax * 100)):
            totals[bin] += len([x for x in value if bin <= int(x * 100) < bin + 1])
            sig_totals[sig_names[value_idx]][bin] = len([x for x in value if bin <= int(x * 100) < bin + 1])
        #plt.hist(values, label=sig_names, color=[colors[sig] for sig in sig_names], bins=int(xmax * 100), stacked=True, weights=weights)
        bottom = [0] * int(xmax * 100)
        tick_labels = ['{:.2f}'.format(i / 100) for i in range(len(totals))]
        ax.set_xticklabels([tick_label for idx, tick_label in enumerate(tick_labels) if idx % 10 == 0])
        ax.set_xticks([pos for pos in range(len(totals)) if pos % 10 == 0])

        for value_idx, value in enumerate(values): # each signature
          sig = sig_names[value_idx]
          new_ys = [sig_totals[sig][bin] / max(1, totals[bin]) for bin in range(len(totals))]
          #bars = ax.bar([round(bin / 100.0, 2) for bin in range(len(totals))], new_ys, label=sig, bottom=bottom, color=colors[sig])
          #tick_label = ['{:.2f}'.format(i / 100) if i % 10 == 0 else '' for i in range(len(totals))]
          bars = ax.bar(x=range(len(totals)), height=new_ys, width=0.99, label=sig, bottom=bottom, color=colors[sig]) #, tick_label=tick_label)
          bottom = [x + y for x,y in zip(bottom, new_ys)]
          #for h in range(len(x)):
          #  bars[h].set_color(colors[sig])
        #ax.xaxis.set_minor_formatter(plt.NullFormatter())
        #ax.grid(which='y')
        #logging.info(np.arange(0, len(totals), 10))
        #ax.xaxis.set_ticks(np.arange(0, len(totals), 10))
      else:
        plt.ylabel('Number of variants')
        plt.hist([sig_ads[n] for n in sig_names], label=sig_names, color=[colors[sig] for sig in sig_names], bins=int(xmax * 100), stacked=True)
      plt.legend()
  else: # just the totals
    plt.ylabel('Number of variants')
    if just_pass:
      plt.hist(ads, bins=int(max(ads) * 100))
    else:
      plt.hist([ads, ads_nopass], label=('PASS', 'No PASS'), bins=int(max(ads + ads_nopass) * 100), stacked=True)
      plt.legend()
  if log:
    plt.yscale('log')
  plt.xlabel('Allele Fraction')
  plt.title(title)
  plt.grid(which='both')
  plt.savefig(target)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Filter VCF')
  parser.add_argument('--sample', required=True,  help='sample name')
  parser.add_argument('--title', required=False, default='Variants seen as a function of allele fraction', help='sample name')
  parser.add_argument('--target', required=True,  help='image output')
  parser.add_argument('--filter', required=False,  help='chromosome or chromosome:range')
  parser.add_argument('--dp', type=int, required=False, default=0, help='minimum dp')
  parser.add_argument('--info_af', action='store_true', help='info af')
  parser.add_argument('--just_pass', action='store_true', help='only plot passes')
  parser.add_argument('--log', action='store_true', help='log on y scale')
  parser.add_argument('--percent', action='store_true', help='show as percentage')
  parser.add_argument('--signature_likelihoods', action='store_true', help='use signature likelihood annotations')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  # sample af vcf
  main(args.sample, args.dp, args.target, args.info_af, args.log, args.filter, args.just_pass, args.signature_likelihoods, args.percent, args.title)
