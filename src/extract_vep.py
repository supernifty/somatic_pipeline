#!/usr/bin/env python
'''
  given tumour and normal vcf pairs, explore msi status
'''

import argparse
import csv
import logging
import sys


def main(vep_header):
  logging.info('starting...')
  vep_fields = vep_header.split('|')
  header = None
  writer = csv.writer(sys.stdout, delimiter='\t')
  for row_count, row in enumerate(csv.reader(sys.stdin, delimiter='\t')):
    if header is None:
      header = row
      csq = header.index('CSQ')
      vep_fields = ['vep_{}'.format(x) for x in vep_fields]
      logging.debug('%i vep columns', len(vep_fields))
      new_header = header[:csq-1] + vep_fields + header[csq+1:]
      writer.writerow(new_header)
      logging.debug('new header has %i columns', len(new_header))
      continue

    for tx in row[csq].split(','):
      vep_cols = tx.split('|')
      if vep_cols[-1] == '1':
        break
    new_row = row[:csq-1] + vep_cols + row[csq+1:]
    writer.writerow(new_row)

    if row_count % 10000 == 0:
      logging.info('%i records read...', row_count)

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--header', required=True, help='vep header')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.header)

