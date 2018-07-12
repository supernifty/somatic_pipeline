#!/usr/bin/env python
'''
  given fastqs what is the maximum possible coverage over a specified region
'''

import argparse
import gzip
import logging
import sys

def main(bed, fastqs):
  logging.info('calculating bed coverage...')
  seen = set()
  current = None
  coords = None
  total = 0
  for line in open(bed, 'r'):
    if line.startswith('#'):
      continue
    fields = line.strip('\n').split('\t')
    if len(fields) < 3:
      continue
    chr = fields[0]
    if chr != current:
      if coords is not None:
        total += len(coords)
        logging.info('%s: %i bases. total: %i', current, len(coords), total)
      current = chr
      logging.info('%s: processing...', current)
      coords = set()
      if current in seen:
        logging.error('bed file appears unsorted')
      else:
        seen.add(current)

    # add range to set
    for x in range(int(fields[1]), int(fields[2])):
      coords.add(x)

  if coords is not None:
    total += len(coords)
    logging.info('%s: %i bases. total: %i', current, len(coords), total)
    coords = None

  logging.info('calculating sequence amount...')
  sequence = 0
  total_min_rl = 1e6
  total_max_rl = -1
  for fastq in fastqs:
    logging.info('%s...', fastq)
    current = 0
    min_rl = 1e6
    max_rl = -1
    for idx, line in enumerate(gzip.open(fastq, 'r')):
      if idx % 4 == 1:
        rl = len(str(line).strip('\n'))
        min_rl = min(min_rl, rl)
        max_rl = max(max_rl, rl)
        current += rl
      if (idx + 1) % 10000000 == 0:
        logging.debug('%s: %i lines. %i sequence. read length: %i to %i', fastq, idx + 1, current, min_rl, max_rl)
    logging.info('%s: %i sequence. read length: %i to %i', fastq, current, min_rl, max_rl)
    sequence += current
    total_min_rl = min(min_rl, total_min_rl)
    total_max_rl = max(max_rl, total_min_rl)

  sys.stdout.write('Total sequence:\t{}\n'.format(sequence))
  sys.stdout.write('Min read length:\t{}\n'.format(total_min_rl))
  sys.stdout.write('Max read length:\t{}\n'.format(total_max_rl))
  sys.stdout.write('Target region:\t{}\n'.format(total))
  sys.stdout.write('Max coverage:\t{:.1f}\n'.format(1.0 * sequence / total))

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Max coverage')
  parser.add_argument('--bed', required=True, help='regions')
  parser.add_argument('--fastqs', required=True, nargs='+', help='fastq files')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.bed, args.fastqs)
