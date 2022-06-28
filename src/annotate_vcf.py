#!/usr/bin/env python
'''
  take fields from vcf
  TODO improve memory usage perhaps with bed filter
'''

import argparse
import collections
import csv
import gzip
import logging
import os.path
import sys

import cyvcf2

CANONICAL='CANONICAL'
CANONICAL_VALUE='YES'

def make_key(variant, alt):
  if variant.REF in 'ACGT' and alt in 'ACGT': # shorter key where possible
    return int(variant.POS) * (1 + 'ACGT'.index(alt))
  else:
    return (int(variant.POS), variant.REF, variant.ALT[0])

def annotate_vcf(annotations, fields, vcf_in, vcf_out, new_names):
    for field in fields:
      vcf_in.add_info_to_header({'ID': new_names.get(field, field), 'Description': 'Annotated field {}'.format(new_names.get(field, field)), 'Type':'Character', 'Number': '1'})
    vcf_out.write(vcf_in.raw_header)
  
    annotated = 0
    count = 0
    seen = set()
    for count, variant in enumerate(vcf_in):
      chr = variant.CHROM.replace('chr', '')
      if chr in annotations:
        key = make_key(variant, variant.ALT[0]) #'{}/{}/{}'.format(variant.POS, variant.REF, variant.ALT[0])
        if count < 100:
          logging.debug('looking for %s', key)
        if key in annotations[chr]:
          for i, field in enumerate([f for f in fields if '|' not in f]): # skip | annotations for now
            #logging.debug('annotating %s at %s with %s', field, key, annotations[chr][key])
            if new_names.get(field, field) in variant.INFO:
              variant.INFO[new_names.get(field, field)] = ','.join(variant.INFO[new_names.get(field, field)], annotations[chr][key][i])
            else:
              variant.INFO[new_names.get(field, field)] = annotations[chr][key][i]
          annotated += 1
      else:
        if chr not in seen:
          seen.add(chr)
          logging.warn('chromosome %s not seen in annotations', chr)
      vcf_out.write(str(variant))
      if (count + 1) % 100000 == 0:
        logging.info('reading %s: %i lines processed %i annotated', chr, count + 1, annotated)
    
    logging.info('done. annotated %i of %i variants', annotated, count)


def main(vcf_in, vcfs, fields, definitions, suffix, rename, no_overwrite):
  logging.info('reading source vcf...')

  new_names = {}
  if rename is not None:
    for r in rename:
      src, dest = r.split('=')
      new_names[src] = dest

  annotations = {}
  for count, variant in enumerate(vcf_in):
    chr = variant.CHROM.replace('chr', '')
    if chr not in annotations:
      annotations[chr] = {}
      logging.info('adding %s to annotations', chr)
    # normal fields
    for alt_idx, alt in enumerate(variant.ALT):
      key = make_key(variant, alt) #'{}/{}/{}'.format(variant.POS, variant.REF, alt)
      annotations[chr][key] = []
      for name in fields:
        if '|' in name:
          continue
        try:
          if isinstance(variant.INFO[name], (list, tuple)):
            annotations[chr][key].append(variant.INFO[name][alt_idx])
          else:
            annotations[chr][key].append(variant.INFO[name])
        except KeyError:
          logging.debug('line %i %s:%i: name %s not found', count, chr, variant.POS, name) # this is no big deal
          annotations[chr][key].append('')
      if count < 100:
        logging.debug('line %i %s:%s: %s', count, chr, key, annotations[chr][key])

    # skip for now
    #key = '{}/{}/{}'.format(variant.POS, variant.REF, variant.ALT[0])
    #annotations[chr][key] = []
    #try:
    #  for name in fields:
    #    if '|' in name and definitions is not None: # subfields
    #      head, tail = name.split('|')
    #      for definition in definitions:
    #        top, rest = definition.split('=')
    #        if top == head: # e.g. CSQ CSQ
    #          items = variant.INFO[head].split(',') # each item
    #          logging.debug('checking %i items from %s', len(items), head)
    #          for item in items:
    #            if item.split('|')[rest.split('|').index(CANONICAL)] == CANONICAL_VALUE: # use canonical
    #              value = item.split('|')[rest.split('|').index(tail)]
    #              logging.debug('annotated %s %s %s with %s...', chr, key, name, value)
    #              annotations[chr][key].append(value)
    #              logging.debug('annotated %s %s %s with %s', chr, key, name, value)
    #              break
    #            else:
    #              pass
    #    else:
    #      annotations[chr][key].append(variant.INFO[name])
    #except KeyError:
    #  logging.debug('line %i %s:%i: extended fields not found', count, chr, variant.POS)

    if (count + 1) % 10000000 == 0:
      logging.info('%i lines...', count + 1)
  logging.debug('reading %s: %i lines processed', chr, count + 1)

  logging.info('reading src vcf: done')

  # now annotate input
  if vcfs is None:
    logging.info('reading from stdin...')
    vcf_in = cyvcf2.VCF('-')
    annotate_vcf(annotations, fields, vcf_in, sys.stdout, new_names)

  else:
    logging.info('processing %i vcfs...', len(vcfs))
    for vcf_count, vcf_fn in enumerate(vcfs):
      vcf_out_fn = vcf_fn.replace('.vcf', '.{}.vcf'.format(suffix))
      if no_overwrite and os.path.isfile(vcf_out_fn):
        logging.info('skipping writing to existing file %s', vcf_out_fn)
        continue

      logging.info('reading from {} and writing to {}...'.format(vcf_fn, vcf_out_fn))
      vcf_in = cyvcf2.VCF(vcf_fn)
      if vcf_out_fn.endswith('.gz'):
        vcf_out = gzip.open(vcf_out_fn, 'wt')
      else:
        vcf_out = open(vcf_out_fn, 'w')
      annotate_vcf(annotations, fields, vcf_in, vcf_out, new_names)
      if vcf_count % 100 == 0:
        logging.debug('processed %i of %i vcfs...', vcf_count, len(vcfs))

def open_file(fn, is_gzipped):
  if is_gzipped:
    return gzip.open(fn, 'rt')
  else:
    return open(fn, 'rt')

def tsv_to_vcf(vcf, chrom_col, pos_col, ref_col, alt_col, delimiter, is_zipped, fields):

  Variant = collections.namedtuple('Variant', 'CHROM POS REF ALT INFO')

  # enumeration a maf into a variant
  logging.debug('reading %s as tsv...', vcf)
  for line, row in enumerate(csv.DictReader(open_file(vcf, is_zipped), delimiter=delimiter)):
    chrom = row[chrom_col].replace('chr', '')
    pos = int(row[pos_col])
    ref = row[ref_col]
    alt = row[alt_col]
    info = {}
    for field in fields:
      info[field] = row[field]

    yield Variant(chrom, pos, ref, (alt,), info)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Annotate VCF with another VCF or TSV')
  parser.add_argument('--vcf', required=True, help='vcf or tsv to annotate')
  parser.add_argument('--is_tsv', action='store_true', help='vcf is actually a tsv')
  parser.add_argument('--tsv_chrom_column', required=False, default='Chromosome', help='tsv chrom column name')
  parser.add_argument('--tsv_pos_column', required=False, default='Start_Position', help='tsv pos column name')
  parser.add_argument('--tsv_ref_column', required=False, default='Reference_Allele', help='tsv ref column name')
  parser.add_argument('--tsv_alt_column', required=False, default='Tumor_Seq_Allele2', help='tsv alt column name')
  parser.add_argument('--tsv_delimiter', required=False, default=',', help='tsv delimiter')
  parser.add_argument('--tsv_zipped', action='store_true', help='is tsv zipped')

  parser.add_argument('--vcfs', required=False, nargs='*', help='vcfs to annotate. stdin if not specified')
  parser.add_argument('--fields', required=True, nargs='+', help='info fields')
  parser.add_argument('--rename', required=False, nargs='*', help='rename field srcname=destname')
  parser.add_argument('--suffix', required=False, default='annot', help='new filename')
  parser.add_argument('--definitions', required=False, nargs='*', help='definitions of fields e.g. CSQ=a|b|...')
  parser.add_argument('--no_overwrite', action='store_true', help='do not overwrite existing vcf')

  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  if args.is_tsv:
    vcf_in = tsv_to_vcf(args.vcf, args.tsv_chrom_column, args.tsv_pos_column, args.tsv_ref_column, args.tsv_alt_column, args.tsv_delimiter, args.tsv_zipped, args.fields)
  else:
    vcf_in = cyvcf2.VCF(args.vcf)
  main(vcf_in, args.vcfs, args.fields, args.definitions, args.suffix, args.rename, args.no_overwrite)

