#!/usr/bin/env python
'''
  summarize the analysis
'''

import argparse
import csv
import logging
import sys

def tsv_to_md(in_fh, out_fh):
  for line_num, line in enumerate(in_fh):
    line = line.strip('\n').replace('\t', ' | ')
    out_fh.write('{}\n'.format(line))
    if line_num == 0:
      out_fh.write('|'.join('-' * (line.count('|') + 1)))
      out_fh.write('\n')
  out_fh.write('\n')

def main(versions, signatures, burden, msi_burden, qc, selected_variants, all_variants):
  logging.info('starting...')

  if qc is not None:
    sys.stdout.write('## QC\n')
    tsv_to_md(open(qc, 'r'), sys.stdout)

  if signatures is not None:
    sys.stdout.write('## COSMIC Signatures\n')
    for line_num, line in enumerate(open(signatures, 'r')):
      line = line.strip('\n').replace('\t', ' | ').replace('Signature.', '')
      sys.stdout.write('{}\n'.format(line))
      if line_num == 0:
        sys.stdout.write('|'.join('-' * (line.count('|') + 1)))
        sys.stdout.write('\n')
    sys.stdout.write('\n')

  if burden is not None and msi_burden is not None:
    sys.stdout.write('## Burden\n')
    sys.stdout.write('\n')

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
  parser.add_argument('--burden', help='mutational burden - exonic snvs and indels')
  parser.add_argument('--msi_burden', help='mutational burden - microsatellite indels')
  parser.add_argument('--qc', help='qc results')
  parser.add_argument('--selected_variants', help='variants in genes of interest')
  parser.add_argument('--all_variants', help='all found variants')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.versions, args.signatures, args.burden, args.msi_burden, args.qc, args.selected_variants, args.all_variants)

