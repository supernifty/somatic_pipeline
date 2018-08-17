#!/usr/bin/env python
'''
  filter tsv file
'''

import argparse
import csv
import logging
import sys

def main(column, values, delimiter='\t'):

  logging.info('reading from stdin...')

  first = True
  accepted = row_idx = 0
  writer = csv.writer(sys.stdout, delimiter=delimiter)

  for row_idx, row in enumerate(csv.reader(sys.stdin, delimiter=delimiter)):
    if first:
      first = False
      column_idx = row.index(column)
      writer.writerow(row)
      continue

    if row[column_idx] in values:
      writer.writerow(row)
      accepted += 1

    if (row_idx + 1) % 100000 == 0:
      logging.info('wrote %i of %i', accepted, row_idx + 1)

  logging.info('wrote %i of %i', accepted, row_idx)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--column', required=True, help='column to filter')
  parser.add_argument('--values', required=True, nargs='+', help='values to match')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
  # sample af vcf
  main(args.column, args.values)
