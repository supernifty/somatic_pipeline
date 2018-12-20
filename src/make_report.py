#!/usr/bin/env python
'''
  summarize the analysis
'''

import argparse
import csv
import logging
import sys

def main(versions, signatures):
  logging.info('starting...')

  if signatures is not None:
    sys.stdout.write('## COSMIC Signatures\n')
    for line_num, line in enumerate(open(signatures, 'r')):
      line = line.strip('\n').replace('\t', ' | ').replace('Signature.', '')
      sys.stdout.write('{}\n'.format(line))
      if line_num == 0:
        sys.stdout.write('|'.join('-' * (line.count('|') + 1)))
        sys.stdout.write('\n')

  if versions is not None:
    sys.stdout.write('## Versions\n')
    sys.stdout.write('Tool|Version\n-|-\n')
    for line in csv.DictReader(open(versions, 'rt'), delimiter='\t'):
      sys.stdout.write('{} | {}\n'.format(line['Tool'], line['Version']))
  
  logging.info('done.')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Extract deamination and oxidation scores from picard\'s output')
  parser.add_argument('--versions', help='provenance')
  parser.add_argument('--signatures', help='signature results')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.versions, args.signatures)


