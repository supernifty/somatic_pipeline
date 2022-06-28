#!/usr/bin/env python
'''
  given a tsv containing a vep field, explode that out into new columns
'''

import argparse
import csv
import logging
import sys


def main(vep_header, transcript, override):
  logging.info('starting...')
  vep_fields = vep_header.split('|')
  transcript_field, transcript_value = transcript.split('=')
  transcript_field_idx = vep_fields.index(transcript_field)

  overrides = {}
  if override is not None:
    symbol_idx = vep_fields.index('SYMBOL')
    for o in override:
      logging.debug('override: %s', o)
      g, t = o.split('=')
      overrides[g] = t.split('|')

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
      for tx in row[csq].split(','): # each transcript
        vep_cols = tx.split('|')
        gene_override = False
        if override is not None:
          if len(vep_cols) > symbol_idx:
            gene = vep_cols[symbol_idx]
            if gene in overrides:
              gene_override = True
              logging.debug('looking for %s in %s', overrides[gene][0], vep_fields)
              if vep_cols[vep_fields.index('vep_{}'.format(overrides[gene][0]))] == overrides[gene][1]: # feature matches override
                found = True
                new_row = row[:csq] + vep_cols + row[csq+1:]
                writer.writerow(new_row)
          else:
            logging.warn('line %i: not enough columns (%i) to extract SYMBOL (%i): %s', row_count, len(vep_cols), symbol_idx, tx)

        if not gene_override and len(vep_cols) > transcript_field_idx and vep_cols[transcript_field_idx] == transcript_value:
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
  parser.add_argument('--override', required=False, default=None, nargs='+', help='override preferred transcript for specific gene using the form Symbol=Fieldname|Value e.g. POLD1=Feature|NM_002691.4')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.header, args.transcript, args.override)

