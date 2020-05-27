#!/usr/bin/env python
'''
  summarize the analysis
'''

import argparse
import collections
import csv
import logging
import sys

#  173348 .
#    342 Benign
#    226 Benign/Likely_benign
#      1 clinvar_pathogenic
#    262 Conflicting_interpretations_of_pathogenicity
#      2 drug_response
#    702 Likely_benign
#    146 Likely_pathogenic
#     60 not_provided
#    352 Pathogenic
#     56 Pathogenic/Likely_pathogenic
#      2 risk_factor
#   1504 Uncertain_significance

CLINVAR_PASS = set(['Conflicting_interpretations_of_pathogenicity', 'drug_response', 'Likely_pathogenic', 'Pathogenic', 'Pathogenic/Likely_pathogenic', 'risk_factor', 'Uncertain_significance'])

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

def main(versions, signatures, burden, msisensor, qc, selected_somatic_variants, all_somatic_variants, all_germline_variants, signature_detail):
  logging.info('starting...')

  samples = collections.defaultdict(dict)
  # get samples from qc
  if qc is not None:
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

  if signatures is not None:
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

  if selected_somatic_variants is not None:
    # somatic variants of interest
    for row in csv.DictReader(open(selected_somatic_variants, 'r'), delimiter='\t'):
      sample, category = row['Filename'].split('_', 1)
      if category not in samples[sample]:
        samples[sample][category] = {}
      if 'selected_somatic_variants' not in samples[sample][category]:
        samples[sample][category]['selected_somatic_variants'] = []
      samples[sample][category]['selected_somatic_variants'].append({'gene': row['vep_SYMBOL'], 'category': row['vep_Consequence'], 'cchange': row['vep_HGVSc'], 'polyphen': row['vep_PolyPhen'], 'sift': row['vep_SIFT'], 'clinvar': row.get('clinvar_pathogenic', 'unavailable'), 'gnomad': 'https://gnomad.broadinstitute.org/variant/{CHROM}-{POS}-{REF}-{ALT}'.format(**row)})

  if all_somatic_variants is not None:
    logging.info('processing %s...', all_somatic_variants)
    # germline variants of interest
    added = 0
    for idx, row in enumerate(csv.DictReader(open(all_somatic_variants, 'r'), delimiter='\t')):
      if idx % 100000 == 0 or added > 0 and added % 100 == 0:
        logging.info('%i variants processed, added %i...', idx, added)
      if row['VCF_SAMPLE_ID'] is None:
        continue # skip missing sample
      if row['GT'] == '0/0' or row['GT'] == './.':
        continue # skip no genotype
      if row.get('clinvar_pathogenic', 'unavailable') not in CLINVAR_PASS:
        continue # skip not annotated as possibly pathogenic
      sample, category = row['VCF_SAMPLE_ID'].split('_', 1)
      if category not in samples[sample]:
        samples[sample][category] = {}
      if 'somatic_variants' not in samples[sample][category]:
        samples[sample][category]['somatic_variants'] = []
      samples[sample][category]['somatic_variants'].append({'gene': row['vep_SYMBOL'], 'category': row['vep_Consequence'], 'cchange': row['vep_HGVSc'], 'polyphen': row['vep_PolyPhen'], 'sift': row['vep_SIFT'], 'clinvar': row.get('clinvar_pathogenic', 'unavailable'), 'cosmic': row['cosmic'], 'gnomad': 'https://gnomad.broadinstitute.org/variant/{CHROM}-{POS}-{REF}-{ALT}'.format(**row)})
      added += 1


  if all_germline_variants is not None:
    logging.info('processing %s...', all_germline_variants)
    # germline variants of interest
    added = 0
    for idx, row in enumerate(csv.DictReader(open(all_germline_variants, 'r'), delimiter='\t')):
      if idx % 1000000 == 0 or added > 0 and added % 100 == 0:
        logging.info('%i variants processed, added %i...', idx, added)
      if row['VCF_SAMPLE_ID'] is None:
        continue # skip missing sample
      if row['GT'] == '0/0' or row['GT'] == './.':
        continue # skip no genotype
      if row.get('clinvar_pathogenic', 'unavailable') not in CLINVAR_PASS:
        continue # skip not annotated as possibly pathogenic
      sample, category = row['VCF_SAMPLE_ID'].split('_', 1)
      if category not in samples[sample]:
        samples[sample][category] = {}
      if 'germline_variants' not in samples[sample][category]:
        samples[sample][category]['germline_variants'] = []
      samples[sample][category]['germline_variants'].append({'gene': row['vep_SYMBOL'], 'category': row['vep_Consequence'], 'cchange': row['vep_HGVSc'], 'polyphen': row['vep_PolyPhen'], 'sift': row['vep_SIFT'], 'clinvar': row.get('clinvar_pathogenic', 'unavailable'), 'gnomad': 'https://gnomad.broadinstitute.org/variant/{CHROM}-{POS}-{REF}-{ALT}'.format(**row)})
      added += 1


  # write out results for each sample
  for sample in samples:
    for category in samples[sample]:
      sys.stdout.write('\n## {} {}\n'.format(sample, category))
      if 'qc' in samples[sample][category]:
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

      if 'somatic_variants' in samples[sample][category]:
        sys.stdout.write('\n### Somatic variants of interest ({})\n'.format(len(samples[sample][category]['somatic_variants'])))
        sys.stdout.write('|Gene|c change|Category|Polyphen|SIFT|Clinvar|COSMIC|Links|\n')
        sys.stdout.write('|-|-|-|-|-|-|-|-|\n')
        items = samples[sample][category]['somatic_variants']
        for item in items:
          sys.stdout.write('|{}|{}|{}|{}|{}|{}|{}|[Gnomad]({})|\n'.format(item['gene'], item['cchange'], item['category'], item['polyphen'], item['sift'], item['clinvar'], item['cosmic'], item['gnomad']))

      if 'germline_variants' in samples[sample][category]:
        sys.stdout.write('\n### Germline variants of interest ({})\n'.format(len(samples[sample][category]['germline_variants'])))
        sys.stdout.write('|Gene|c change|Category|Polyphen|SIFT|Clinvar|Links|\n')
        sys.stdout.write('|-|-|-|-|-|-|-|\n')
        items = samples[sample][category]['germline_variants']
        for item in items:
          sys.stdout.write('|{}|{}|{}|{}|{}|{}|[Gnomad]({})|\n'.format(item['gene'], item['cchange'], item['category'], item['polyphen'], item['sift'], item['clinvar'], item['gnomad']))

      if 'selected_somatic_variants' in samples[sample][category]:
        sys.stdout.write('\n### Somatic variants of interest ({})\n'.format(len(samples[sample][category]['selected_somatic_variants'])))
        sys.stdout.write('|Gene|c change|Category|Polyphen|SIFT|Clinvar|Links|\n')
        sys.stdout.write('|-|-|-|-|-|-|-|\n')
        items = samples[sample][category]['selected_somatic_variants']
        for item in items:
          sys.stdout.write('|{}|{}|{}|{}|{}|{}|[Gnomad]({})|\n'.format(item['gene'], item['cchange'], item['category'], item['polyphen'], item['sift'], item['clinvar'], item['gnomad']))

  if versions is not None:
    sys.stdout.write('## Versions\n')
    sys.stdout.write('Tool|Version\n-|-\n')
    for line in csv.DictReader(open(versions, 'rt'), delimiter='\t'):
      sys.stdout.write('{} | {}\n'.format(line['Tool'], line['Version']))
  
  logging.info('done.')

#    "src/make_report.py --versions {input.versions} --signatures {input.signatures} --burden {input.burden} --msi_burden {input.msi_burden} --qc {input.qc} --selected_somatic_variants {input.selected_somatic_variants} --all_somatic_variants {input.all_somatic_variants} > {output.md} 2>{log.stderr} && "
if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Extract deamination and oxidation scores from picard\'s output')
  parser.add_argument('--versions', help='provenance')
  parser.add_argument('--signatures', help='signature results')
  parser.add_argument('--signature_detail', help='signature results')
  parser.add_argument('--burden', help='mutational burden - exonic snvs and indels')
  parser.add_argument('--msisensor', help='mutational burden - microsatellite indels')
  parser.add_argument('--qc', help='qc results')
  parser.add_argument('--selected_somatic_variants', help='variants in genes of interest')
  parser.add_argument('--all_somatic_variants', help='all found variants')
  parser.add_argument('--all_germline_variants', help='all found variants')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.versions, args.signatures, args.burden, args.msisensor, args.qc, args.selected_somatic_variants, args.all_somatic_variants, args.all_germline_variants, args.signature_detail)

