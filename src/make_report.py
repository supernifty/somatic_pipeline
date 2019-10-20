#!/usr/bin/env python
'''
  summarize the analysis
'''

import argparse
import collections
import csv
import logging
import sys

def highlight_if(condition, msg):
  if condition:
    return '*{}*'.format(msg)
  else:
    return msg

def message_if(condition, msg_yes, msg_no):
  if condition:
    return msg_yes
  else:
    return msg_no

def convert_range(v, ranges, texts, otherwise):
  for r, t in zip(ranges, texts):
    if v < r:
      return t
  return otherwise

def tsv_to_md(in_fh, out_fh):
  for line_num, line in enumerate(in_fh):
    line = line.strip('\n').replace('\t', ' | ')
    out_fh.write('{}\n'.format(line))
    if line_num == 0:
      out_fh.write('|'.join('-' * (line.count('|') + 1)))
      out_fh.write('\n')
  out_fh.write('\n')

def main(versions, signatures, burden, msisensor, qc, selected_variants, all_variants, signature_detail):
  logging.info('starting...')

  samples = collections.defaultdict(dict)
  # get samples from qc
  for row in csv.DictReader(open(qc, 'r'), delimiter='\t'):
    sample, category = row['Sample'].split('_', 1)
    if category not in samples[sample]:
      samples[sample][category] = {}
    samples[sample][category]['qc'] = row['Assessment']

  if burden is not None:
    for row in csv.DictReader(open(burden, 'r'), delimiter='\t'):
      sample, category = row['Filename'].split('/')[-1].split('.')[0].split('_', 1)
      if category not in samples[sample]:
        samples[sample][category] = {}
      if 'genome' not in samples[sample][category]:
        samples[sample][category]['genome'] = []
      samples[sample][category]['genome'].append('TMB: {} ({})'.format(row['PerMB'], convert_range(float(row['PerMB']), [10, 100], ['Normal', 'Hypermutated'], 'Ultrahypermutated'))) 
      

  if msisensor is not None:
    for row in csv.DictReader(open(msisensor, 'r'), delimiter='\t'):
      sample, category = row['Sample'].split('_', 1)
      if category not in samples[sample]:
        samples[sample][category] = {}
      if 'genome' not in samples[sample][category]:
        samples[sample][category]['genome'] = []
      samples[sample][category]['genome'].append('MSI Sensor: {}% ({})'.format(row['%'], row['Class']))
 
  # signatures
  sig_map = {}
  if signature_detail is not None:
    for row in csv.DictReader(open(signature_detail, 'r'), delimiter='\t'):
      sig_map[row['Signature']] = row['Summary']

  for row in csv.DictReader(open(signatures, 'r'), delimiter='\t'):
    sample, category = row['Filename'].split('_', 1)
    if category not in samples[sample]:
      samples[sample][category] = {}
    samples[sample][category]['sigs'] = []
    # add number of mutations, error, top 3 sigs
    samples[sample][category]['sigs'].append('Number of mutations: {} ({})'.format(row['Mutations'], message_if(int(row['Mutations']) > 50, 'OK', '*Low*')))
    samples[sample][category]['sigs'].append('Reconstruction error: {} ({})'.format(row['Error'], message_if(float(row['Error']) < 0.1, 'OK', '*High*')))
    sorted_sbs = [(y[0], y[1]) for y in sorted([x for x in row.items() if x[0].startswith('SBS') and x[0][3].isdigit()], key=lambda item: -float(item[1]))]
    for idx in range(min(len(sorted_sbs), 3)):
      samples[sample][category]['sigs'].append('{} {:.0f}% (associated with {})'.format(sorted_sbs[idx][0], round(float(sorted_sbs[idx][1]) * 100, 0), sig_map.get(sorted_sbs[idx][0], 'unspecified')))

  # somatic variants of interest
  for row in csv.DictReader(open(selected_variants, 'r'), delimiter='\t'):
    sample, category = row['Filename'].split('_', 1)
    if category not in samples[sample]:
      samples[sample][category] = {}
    if 'variants' not in samples[sample][category]:
      samples[sample][category]['variants'] = []
    samples[sample][category]['variants'].append({'gene': row['vep_SYMBOL'], 'category': row['vep_Consequence'], 'cchange': row['vep_HGVSc'], 'polyphen': row['vep_PolyPhen'], 'sift': row['vep_SIFT'], 'gnomad': 'https://gnomad.broadinstitute.org/variant/{CHROM}-{POS}-{REF}-{ALT}'.format(**row)})
    
  # write out results for each sample
  for sample in samples:
    for category in samples[sample]:
      sys.stdout.write('\n## {} {}\n'.format(sample, category))
      sys.stdout.write('### QC\n')
      assessment = samples[sample][category]['qc']
      sys.stdout.write('* {}\n'.format(highlight_if(assessment != 'OK', assessment)))
      
      if 'sigs' in samples[sample][category]:
        sys.stdout.write('\n### COSMIC Signatures\n')
        items = samples[sample][category]['sigs']
        for item in items:
          sys.stdout.write('* {}\n'.format(item))
      
      if 'genome' in samples[sample][category]:
        sys.stdout.write('\n### Burden\n')
        items = samples[sample][category]['genome']
        for item in items:
          sys.stdout.write('* {}\n'.format(item))

      if 'variants' in samples[sample][category]:
        sys.stdout.write('\n### Variants of interest ({})\n'.format(len(samples[sample][category]['variants'])))
        sys.stdout.write('|Gene|c change|Category|Polyphen|SIFT|Links|\n')
        sys.stdout.write('|-|-|-|-|-|-|\n')
        items = samples[sample][category]['variants']
        for item in items:
          sys.stdout.write('|{}|{}|{}|{}|{}|[Gnomad]({})|\n'.format(item['gene'], item['cchange'], item['category'], item['polyphen'], item['sift'], item['gnomad']))

  return

  if selected_variants is not None:
    sys.stdout.write('## Selected Variants\n')
    tsv_to_md(open(selected_variants, 'r'), sys.stdout)
    
    sys.stdout.write('\n')

  if all_variants is not None:
    sys.stdout.write('## All Variants\n')
    sys.stdout.write('\n')
  

  if versions is not None:
    sys.stdout.write('## Versions\n')
    sys.stdout.write('Tool|Version\n-|-\n')
    for line in csv.DictReader(open(versions, 'rt'), delimiter='\t'):
      sys.stdout.write('{} | {}\n'.format(line['Tool'], line['Version']))
  
  logging.info('done.')

#    "src/make_report.py --versions {input.versions} --signatures {input.signatures} --burden {input.burden} --msi_burden {input.msi_burden} --qc {input.qc} --selected_variants {input.selected_variants} --all_variants {input.all_variants} > {output.md} 2>{log.stderr} && "
if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Extract deamination and oxidation scores from picard\'s output')
  parser.add_argument('--versions', help='provenance')
  parser.add_argument('--signatures', help='signature results')
  parser.add_argument('--signature_detail', help='signature results')
  parser.add_argument('--burden', help='mutational burden - exonic snvs and indels')
  parser.add_argument('--msisensor', help='mutational burden - microsatellite indels')
  parser.add_argument('--qc', help='qc results')
  parser.add_argument('--selected_variants', help='variants in genes of interest')
  parser.add_argument('--all_variants', help='all found variants')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.versions, args.signatures, args.burden, args.msisensor, args.qc, args.selected_variants, args.all_variants, args.signature_detail)

