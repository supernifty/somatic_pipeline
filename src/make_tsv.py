#!/usr/bin/env python
'''
  make a tsv from values
'''

import argparse
import logging
import sys

def main(columns, rows):
  logging.info('starting...')

  sys.stdout.write('{}\n'.format('\t'.join(columns)))
  for row in rows:
    fields = row.split(',')
    if len(fields) != len(columns):
      logging.warn('%s has %i fields. Expected %i', row, len(fields), len(columns))
    sys.stdout.write('{}\n'.format('\t'.join(fields)))

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--columns', nargs='+', required=True, help='column names')
  parser.add_argument('--rows', nargs='+', help='rows')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.columns, args.rows)
