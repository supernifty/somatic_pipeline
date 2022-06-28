#!/usr/bin/env python
'''
  calculate and plot af from ad
'''

import argparse
import collections
import logging
import math
import random
import sys

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import rcParams

from matplotlib.offsetbox import AnchoredText

import cyvcf2

SQUIGGEM=0.013
YSQUIGGEM=0.01
DPI=300

# mutect pass, mutect no pass, strelka pass, strelka no pass
AF_COLORS=['#92b3eb', '#8293ab', '#83e063', '#63a033']

colors = {"SBS1": "#de3860", "SBS2": "#41ac2f", "SBS3": "#7951d0", "SBS4": "#73d053", "SBS5": "#b969e9", "SBS6": "#91ba2c", "SBS7a": "#b4b42f", "SBS7b": "#5276ec", "SBS7c": "#daae36", "SBS7d": "#9e40b5", "SBS8": "#43c673", "SBS9": "#dd4cb0", "SBS10a": "#3d9332", "SBS10b": "#de77dd", "SBS11": "#7bad47", "SBS12": "#9479e8", "SBS13": "#487b21", "SBS14": "#a83292", "SBS15": "#83c67d", "SBS16": "#664db1", "SBS17a": "#e18d28", "SBS17b": "#588de5", "SBS18": "#e2672a", "SBS19": "#34c7dd", "SBS20": "#cf402b", "SBS21": "#5acdaf", "SBS22": "#d74587", "SBS23": "#428647", "SBS24": "#7b51a7", "SBS25": "#b4ba64", "SBS26": "#646cc1", "SBS27": "#a27f1f", "SBS28": "#3b63ac", "SBS29": "#dca653", "SBS30": "#505099", "SBS31": "#7d8529", "SBS32": "#bf8ade", "SBS33": "#516615", "SBS34": "#b65da7", "SBS35": "#57a87a", "SBS36": "#c84249", "SBS37": "#37b5b1", "SBS38": "#a14622", "SBS39": "#58b5e1", "SBS40": "#ba6e2f", "SBS41": "#589ed8", "SBS42": "#e98261", "SBS43": "#3176ae", "SBS44": "#656413", "SBS45": "#a19fe2", "SBS46": "#756121", "SBS47": "#7e4a8d", "SBS48": "#326a38", "SBS49": "#dd8abf", "SBS50": "#1a6447", "SBS51": "#e78492", "SBS52": "#30876c", "SBS53": "#9d4d7c", "SBS54": "#919d5b", "SBS55": "#9d70ac", "SBS56": "#5b6f34", "SBS57": "#65659c", "SBS58": "#c9a865", "SBS59": "#a1455d", "SBS60": "#5e622c", "SBS84": "#b66057", "SBS85": "#dca173", "DBS1": "#855524", "DBS2": "#9f7846", "DBS3": "#7951d0", "DBS4": "#73d053", "DBS5": "#b969e9", "DBS6": "#91ba2c", "DBS7": "#3656ca", "DBS8": "#b4b42f", "DBS9": "#5276ec", "DBS10": "#daae36", "DBS11": "#9e40b5", "ID1": "#de3860", "ID2": "#41ac2f", "ID3": "#7951d0", "ID4": "#73d053", "ID5": "#b969e9", "ID6": "#91ba2c", "ID7": "#9e40b5", "ID8": "#43c673", "ID9": "#dd4cb0", "ID10": "#3d9332", "ID11": "#de77dd", "ID12": "#7bad47", "ID13": "#9479e8", "ID14": "#487b21", "ID15": "#a83292", "ID16": "#83c67d", "ID17": "#664db1"}

H = '0123456789ABCDEF'
def random_color():
  return '#' + ''.join([H[random.randint(0, 15)] for _ in range(6)])

def main(samples, dp_threshold, target, info_af, log, filter, just_pass, use_likelihoods, percent, title, genes, consequences, vep_format, impacts, gene_colors, annotate, vcfs, vcf_names, width, height, annotate_graph):
  logging.info('width, height = %i, %i', width, height)
  rcParams['figure.figsize'] = width, height

  if vcfs is None:
    logging.info('reading from stdin...')
    vcfs = ["-"]
    vcf_names = [""]

  vcf_ads = {}

  if genes is not None:
    genes = set(genes)
  
  if consequences is not None:
    consequences = set(consequences)
  
  if impacts is not None:
    impacts = set(impacts)
  
  if vep_format is not None:
    vep_format = vep_format.split('|')
  
  vafs = []
  
  for vcf, name, sample in zip(vcfs, vcf_names, samples):
    logging.info('reading %s...', vcf)
    vcf_in = cyvcf2.VCF(vcf)  
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
  
      try:
        value = variant.INFO["AF"]
      except:
        try:
          value = variant.INFO["VAF"]
        except:
          ad_samples = variant.format("AD")
          if ad_samples is None:
            value = 0
          else:
            ad = ad_samples[sample_id]
            ref = ad[0]
            alt = ad[1]
            if ref + alt > 0:
              value = alt / (ref + alt)
            else:
              value = 0
      target_ad.append(value) # overall set of afs
  
      # if this a variant of interest?
      if genes is not None and len(genes) > 0:
        # parse vep
        veps = variant.INFO['CSQ'].split(',')
        for vep in veps:
          values = {x[0]: x[1] for x in zip(vep_format, vep.split('|'))}
          #logging.debug('vep: %s', values)
          if values['PICK'] == '1' and values['SYMBOL'] in genes and (consequences is None or values['Consequence'] in consequences) and (impacts is None or values['IMPACT'] in impacts):
            logging.debug('interesting variant %s', variant)
            vafs.append({'gene': values['SYMBOL'], 'HGVSc': values['HGVSc'], 'HGVSp': values['HGVSp'], 'vaf': value})
            break
  
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
    vcf_ads[name] = (ads, ads_nopass)
    logging.info("finished reading %s", name)

  # done with vcf
  all_ads = []
  pass_ads = []
  for name in vcf_ads.keys():
    all_ads.extend(vcf_ads[name][0] + vcf_ads[name][1])
    pass_ads.extend(vcf_ads[name][0])

  min_ads = min(all_ads)
  max_ads = max(all_ads)
  logging.info('processed %i variants. no pass %i. low af %i. low dp %i. allowed %i AF range %.2f to %.2f', variant_count + 1, skipped_pass, skipped_af, skipped_af, allowed, min_ads, max_ads)

  if just_pass:
    xmax = max(pass_ads + [0.01])
  else:
    xmax = max(all_ads + [0.01])

  # now plot histogram
  if use_likelihoods:
    sig_names = sorted(sig_ads.keys())
    if len(sig_ads) == 0: # no variants do a dummy
      yh, xh, _ = plt.hist([0.0])
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

          if sig not in colors:
            colors[sig] = random_color()
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
        for sig in sig_names:
          if sig not in colors:
            colors[sig] = random_color()
        yh, xh, _ = plt.hist([sig_ads[n] for n in sig_names], label=sig_names, color=[colors[sig] for sig in sig_names], bins=int(xmax * 100), stacked=True)
  else: # just the totals
    plt.ylabel('Number of variants')
    if just_pass and len(pass_ads) > 0:
      yh, xh, _ = plt.hist(pass_ads, bins=int(max(ads) * 100)) # todo this won't work with multiple vcfs
    else:
      labels = []
      ys = []
      for key in vcf_ads:
        labels.append('{} pass'.format(key))
        labels.append('{} non-pass'.format(key))
        ys.append(vcf_ads[key][0])
        ys.append(vcf_ads[key][1])
      yh, xh, _ = plt.hist(ys, label=labels, bins=int(max(all_ads) * 100), stacked=True, color=AF_COLORS[:len(ys)])
  if log:
    plt.yscale('log')

  # vafs of interest
  if annotate is not None and len(annotate) > 0:
    for a in annotate:
      if '@' not in a:
        logging.warn('annotation %s did not contain @', a)
        continue
      message, vaf = a.split('@')
      # vaf can have a color too
      if ':' in vaf:
        vaf, color = vaf.split(':')
      else:
        color = None

      if float(vaf) > 1:
        vafs.append({'custom': message, 'vaf': float(vaf)/100, 'color': color})
      else:
        vafs.append({'custom': message, 'vaf': float(vaf), 'color': color})

  # vafs.append({'gene': values['SYMBOL'], 'HGVSp': values['HGVSp'], 'vaf': value})
  # sorted(list_to_be_sorted, key=lambda k: k['name']) 
  labelled = set()
  for i, vaf in enumerate(sorted(vafs, key=lambda k: k['vaf'])):
    logging.debug('adding %s', vaf)
    if 'color' in vaf and vaf['color'] is not None:
      color=vaf['color']
    elif gene_colors is None or 'gene' not in vaf:
      color='#800000'
    else: # look up color based on gene
      for gc in gene_colors:
        g, c = gc.split('=')
        if g == vaf['gene']:
          color=c
          break
    if 'custom' in vaf:
      plt.axvline(vaf['vaf'], 0, 1, color=color)
    elif vaf['gene'] not in labelled:
      plt.axvline(vaf['vaf'], 0, 1, color=color, label=vaf['gene'])
      labelled.add(vaf['gene'])
    else:
      plt.axvline(vaf['vaf'], 0, 1, color=color)

    if log:
      ypos = math.pow(i % 3 + 1, 10) * yh.max() / 1000000 + 2
    else:
      ypos = i % 4 / 4 * yh.max() + yh.max() * YSQUIGGEM 

    logging.debug('ypos %s', ypos)
    if 'custom' in vaf:
      annot = vaf['custom']
      plt.text(float(vaf['vaf'])-SQUIGGEM * xmax, ypos, '{}'.format(annot), rotation=90, verticalalignment='baseline', fontsize='smaller', color='black')
      continue
    elif ':' in vaf['HGVSp']:
      annot = vaf['HGVSp'].split(':')[1]
    elif ':' in vaf['HGVSc']:
      annot = vaf['HGVSc'].split(':')[1]
    else:
      annot = ''
    plt.text(float(vaf['vaf'])-SQUIGGEM * xmax, ypos, '{}:{} {}%'.format(vaf['gene'], annot, '{:.0f}'.format(vaf['vaf'] * 100)), rotation=90, verticalalignment='baseline', fontsize='smaller', color='black')

  if annotate_graph is not None:
    text_box = AnchoredText(annotate_graph, frameon=True, loc=4, pad=0.5)
    plt.setp(text_box.patch, facecolor='white', alpha=0.5)
    ax = plt.axes()
    ax.add_artist(text_box)

  plt.legend()
  plt.xlabel('Allele Fraction')
  plt.title(title)
  plt.grid(which='both')
  plt.tight_layout()
  plt.savefig(target, dpi=DPI)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Filter VCF')
  parser.add_argument('--sample', nargs='+', required=True,  help='sample names')
  parser.add_argument('--title', required=False, default='Variants seen as a function of allele fraction', help='sample name')
  parser.add_argument('--target', required=True,  help='image output')
  parser.add_argument('--filter', required=False,  help='chromosome or chromosome:range')
  parser.add_argument('--dp', type=int, required=False, default=0, help='minimum dp')
  parser.add_argument('--info_af', action='store_true', help='info af')
  parser.add_argument('--just_pass', action='store_true', help='only plot passes')
  parser.add_argument('--log', action='store_true', help='log on y scale')
  parser.add_argument('--percent', action='store_true', help='show as percentage')
  parser.add_argument('--signature_likelihoods', action='store_true', help='use signature likelihood annotations')
  parser.add_argument('--genes', nargs='*', required=False, help='genes to annotate with vaf lines')
  parser.add_argument('--gene_colors', nargs='*', required=False, help='genes to annotate with vaf lines')
  parser.add_argument('--consequences', nargs='*', required=False, help='variants with these consequence will be included')
  parser.add_argument('--impacts', nargs='*', required=False, help='variants with these impacts will be included')
  parser.add_argument('--vep_format', required=False, help='vep format')
  parser.add_argument('--annotate', nargs='*', required=False, help='custom annotation of the form annotation@af[:color]')
  parser.add_argument('--annotate_graph', required=False, help='text to place at bottom right')
  parser.add_argument('--vcfs', required=False, nargs='+', help='vcfs to read')
  parser.add_argument('--vcf_names', required=False, nargs='+', help='vcfs names')
  parser.add_argument('--width', required=False, default=16, type=int, help='width')
  parser.add_argument('--height', required=False, default=12, type=int, help='width')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  # sample af vcf
  main(args.sample, args.dp, args.target, args.info_af, args.log, args.filter, args.just_pass, args.signature_likelihoods, args.percent, args.title, args.genes, args.consequences, args.vep_format, args.impacts, args.gene_colors, args.annotate, args.vcfs, args.vcf_names, args.width, args.height, args.annotate_graph)
