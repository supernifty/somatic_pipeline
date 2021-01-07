#!/usr/bin/env python
'''
  given a tsv containing a vep field, explode that out into new columns
'''

import argparse
import csv
import logging
import sys


def main(vep_header, transcript):
  logging.info('starting...')
  vep_fields = vep_header.split('|')
  transcript_field, transcript_value = transcript.split('=')
  transcript_field_idx = vep_fields.index(transcript_field)
  header = None
  writer = csv.writer(sys.stdout, delimiter='\t')
  for row_count, row in enumerate(csv.reader(sys.stdin, delimiter='\t')):
    if header is None:
      header = row
      csq = header.index('CSQ')
      vep_fields = ['vep_{}'.format(x) for x in vep_fields]
      logging.debug('%i vep columns', len(vep_fields))
      new_header = header[:csq] + vep_fields + header[csq+1:]
      writer.writerow(new_header)
      logging.debug('new header has %i columns: %s', len(new_header), new_header)
      continue

    if csq < len(row):
      found = False
      for tx in row[csq].split(','):
        vep_cols = tx.split('|')
        if len(vep_cols) > transcript_field_idx and vep_cols[transcript_field_idx] == transcript_value:
          # report multiple canonicals
          found = True
          new_row = row[:csq] + vep_cols + row[csq+1:]
          writer.writerow(new_row)

      if not found:
        logging.warn('line %i: transcript matching %s not found. writing anyway.', row_count, transcript)
        new_row = row[:csq] + vep_cols + row[csq+1:]
        writer.writerow(new_row)

    else:
      logging.warn('skipping line %i: only %i rows, need %i', row_count, len(row), csq + 1)

    if row_count % 10000 == 0:
      logging.info('%i records read...', row_count)

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Explode VEP fields')
  parser.add_argument('--header', required=True, help='vep header')
  parser.add_argument('--transcript', required=False, default='PICK=1', help='which transcript field=value')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.header, args.transcript)

