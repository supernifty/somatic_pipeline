#!/usr/bin/env python
'''
  combine batch results into single array
'''

import argparse
import collections
import logging
import os
import re
import sys

import csv

V2_SBS=False
V3_SBS=True
V3_DBS=True
V31_SBS=False

V2_ID=False
V3_ID=True
V31_ID=False

TMB_CLEANED=False

GENES_OF_INTEREST = set([
  'MUTYH',
  'NRAS',
  'MAP3K21',
  'MSH2',
  'MSH6',
  'ACVR2A',
  'CASP8',
  'TGFBR2',
  'MLH1',
  'CTNNB1',
  'PIK3CA',
  'FBXW7',
  'MSH3',
  'APC',
  'PMS2',
  'BRAF',
  'WRN',
  'DKK1',
  'PTEN',
  'HRAS',
  'IGF2',
  'ATM',
  'KRAS',
  'MUC19',
  'IGF1',
  'POLE',
  'BRCA2',
  'FAN1',
  'NTHL1',
  'TP53',
  'BRCA1',
  'RNF43',
  'AXIN2',
  'SMAD2',
  'SMAD4',
  'POLD1'])

def add_signature(fn, samples, header, source, prefix):
  for row in csv.DictReader(open(fn, 'r'), delimiter='\t'):
    sample = row['Filename']
    del row['Filename'] # include everything but
    if 'Error' in row:
      row['SignatureError'] = row['Error']
      del row['Error']

    # figure the top sigs
    signature_names = [x for x in row.keys() if x not in ['Mutations', 'SignatureError']]
    signature_values = [float(row[signature_name]) for signature_name in signature_names]
    top = sorted(zip(signature_values, signature_names), reverse=True)

    for idx, x in enumerate(top[:5]):
      #row['sigs'] = ' '.join(['{} ({})'.format(x[1].replace('Signature.', '').replace('SBS', '').replace('ID', ''), x[0]) for x in top[:5]])
      row['R{}'.format(idx + 1)] = '{} {}%'.format(x[1], int(x[0] * 100)) # as %

    # also remove original signatures
    for key in signature_names:
      del row[key]

    # add prefix to all row names
    for key in row.copy():
      row['{}{}'.format(prefix, key)] = row[key]
      del row[key]

    logging.debug('add_signature: adding %s to %s', source, sample)
    samples['{}/{}'.format(sample, source)].update(row) # add results
    samples['{}/{}'.format(sample, source)]['source'] = source
    header.update(row.keys())
    logging.debug('add_signature: header is now %s', header)

def display_value(samples, key, name):
  if name == 'LOH':
    return ' '.join(sorted(list(samples[key].get(name, set()))))
  else:
    return samples[key].get(name, 'NA') 

def main(directories, phenotype, require):
  logging.info('starting...')

  samples = collections.defaultdict(dict)
  header = set()

  for directory in directories:
    logging.info('parsing %s...', directory)
    source = os.path.basename(directory)

    # mutational_signatures_v2.filter.combined.tsv
    # mutational_signatures_v3_dbs.combined.tsv
    # mutational_signatures_v3_sbs.filter.combined.tsv
    # mutational_signatures_v3_id_strelka.filter.combined.tsv

    # targetted gene summary
    #fn = os.path.join(directory, 'out', 'aggregate', 'targetted_gene_summary.tsv')
    #if not os.path.isfile(fn):
    #  logging.info('skipping %s', directory)
    #else: 
    #  for row in csv.DictReader(open(fn, 'r'), delimiter='\t'):


    # signatures - v2
    if V2_SBS:
      fn = os.path.join(directory, 'out', 'aggregate', 'mutational_signatures_v2.filter.combined.tsv')
      if not os.path.isfile(fn):
        logging.info('skipping %s: no v2 sbs', directory)
        continue
      add_signature(fn, samples, header, source, 'v2')

    # v3 sbs signatures
    if V3_SBS:
      fn = os.path.join(directory, 'out', 'aggregate', 'mutational_signatures_v3_sbs_capture.filter.combined.tsv')
      if not os.path.isfile(fn):
        logging.info('skipping %s: no v3 sbs', directory)
        continue
      add_signature(fn, samples, header, source, 'SBS.')

    # v3 id signatures
    if V3_ID:
      fn = os.path.join(directory, 'out', 'aggregate', 'mutational_signatures_v3_id_strelka.filter.combined.tsv')
      #fn = os.path.join(directory, 'out', 'aggregate', 'mutational_signatures_v3_id_capture.filter.combined.tsv')
      if not os.path.isfile(fn):
        logging.info('skipping %s: no v3 id', directory)
        continue
      add_signature(fn, samples, header, source, 'ID.')

    if V3_DBS:
      fn = os.path.join(directory, 'out', 'aggregate', 'mutational_signatures_v3_dbs.filter.combined.tsv')
      #fn = os.path.join(directory, 'out', 'aggregate', 'mutational_signatures_v3_id_capture.filter.combined.tsv')
      if not os.path.isfile(fn):
        logging.info('skipping %s: no v3 id', directory)
        continue
      add_signature(fn, samples, header, source, 'DB.')

    # v3 sbs signatures
    if V31_SBS:
      fn = os.path.join(directory, 'out', 'aggregate', 'mutational_signatures_v3.1_sbs.combined.tsv')
      if not os.path.isfile(fn):
        logging.info('skipping %s: no v3.1 sbs', directory)
        continue
      add_signature(fn, samples, header, source, 'SBS.')

    # v3 id signatures
    if V31_ID:
      #fn = os.path.join(directory, 'out', 'aggregate', 'mutational_signatures_v3_id_strelka.filter.combined.tsv')
      #fn = os.path.join(directory, 'out', 'aggregate', 'mutational_signatures_v3.1_id.combined.tsv')
      fn = os.path.join(directory, 'out', 'aggregate', 'mutational_signatures_v3.1_id.combined.tsv')
      if not os.path.isfile(fn):
        logging.info('skipping %s no v3.1 id', directory)
        continue
      add_signature(fn, samples, header, source, 'ID.')


    # tmb
    fn = os.path.join(directory, 'out', 'aggregate', 'mutation_rate.tsv')
    for row in csv.DictReader(open(fn, 'r'), delimiter='\t'):
      sample = row['Filename'].split('/')[-1].split('.')[0]
      del row['Filename'] # include everything but
      samples['{}/{}'.format(sample, source)]['TMB'] = row['PerMB']
      logging.debug('add_signature: adding %s to %s', source, sample)
      header.add('TMB')

    # tmb with signature artefacts removed
    if TMB_CLEANED:
      fn = os.path.join(directory, 'out', 'aggregate', 'mutation_rate.artefact_filter.tsv')
      for row in csv.DictReader(open(fn, 'r'), delimiter='\t'):
        sample = row['Filename'].split('/')[-1].split('.')[0]
        del row['Filename'] # include everything but
        samples['{}/{}'.format(sample, source)]['TMB.cleaned'] = row['PerMB']
        header.add('TMB.cleaned')

    # msisensor
    fn = os.path.join(directory, 'out', 'aggregate', 'msisensor.tsv')
    try:
      for row in csv.DictReader(open(fn, 'r'), delimiter='\t'):
        sample = row['Sample']
        del row['Sample'] # include everything but
        samples['{}/{}'.format(sample, source)]['MSISensor'] = row['%']
        header.add('MSISensor')
    except:
      logging.error('failed to add msisensor')

    # mantis
    fn = os.path.join(directory, 'out', 'aggregate', 'mantis.tsv')
    try:
      for row in csv.DictReader(open(fn, 'r'), delimiter='\t'):
        sample = row['Sample'].split('/')[-1].split('.')[0]
        del row['Sample'] # include everything but
        samples['{}/{}'.format(sample, source)]['Mantis'] = row['Pct']
        header.add('Mantis')
    except:
      logging.error('failed to add mantis')

    # msiseq
    fn = os.path.join(directory, 'out', 'aggregate', 'msiseq.tsv')
    try:
      for row in csv.DictReader(open(fn, 'r'), delimiter='\t'):
        sample = row['Sample']
        del row['Sample'] # include everything but
        samples['{}/{}'.format(sample, source)]['msiseq'] = row['S.ind']
        header.add('msiseq')
    except:
      logging.error('failed to add msiseq')


    # ontarget coverage
    try:
      fn = os.path.join(directory, 'out', 'aggregate', 'ontarget.tsv')
      # first line can contain spaces in older pipeline
      lines = [re.sub('  *', '\t', x) for x in open(fn, 'r').readlines()]
      
      for row in csv.DictReader(lines, delimiter='\t'):
        sample = row['Filename'].split('/')[-1].split('.')[0]
        del row['Filename'] # include everything but
        samples['{}/{}'.format(sample, source)]['MeanOnTargetCoverage'] = row['Mean']
        header.add('MeanOnTargetCoverage')
    except:
      logging.error('failed to add ontarget')
      raise

    # loh
    try:
      fn = os.path.join(directory, 'out', 'aggregate', 'loh.genes.tsv')
      for row in csv.DictReader(open(fn, 'r'), delimiter='\t'): # sample gene accept
        if row['gene'] in GENES_OF_INTEREST:
          if 'LOH' not in samples['{}/{}'.format(row['sample'], source)]:
            samples['{}/{}'.format(row['sample'], source)]['LOH'] = set()
          samples['{}/{}'.format(row['sample'], source)]['LOH'].add(row['gene'])
          header.add('LOH')
    except:
      logging.error('failed to add loh')

  if phenotype is not None:
    header.add('Phenotype')
    header.add('Category')
    for row in csv.DictReader(open(phenotype, 'r'), delimiter='\t'):
      for key in samples.keys():
        if key.startswith('{}/'.format(row['Sample Name'])):
          samples[key]['Phenotype'] = row['Phenotype']
          samples[key]['Category'] = row['Category']
          logging.debug('adding phenotype for %s', key)
        else:
          pass #logging.info('skipping phenotype for %s', row['Sample Name'])
  else:
    logging.info('Not adding phenotype data')

  # now write everything to stdout
  sys.stdout.write('Sample\tSource\t{}\n'.format('\t'.join(sorted(list(header)))))
  for key in sorted(samples.keys()):
    if require is not None and require not in samples[key]:
      logging.debug('skipping %s...', key)
      continue
    logging.debug('writing %s...', key)
    sample, source = key.split('/')
    sys.stdout.write('{}\t{}\t{}\n'.format(sample, source, '\t'.join([display_value(samples, key, name)for name in sorted(list(header))])))

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='combine batch results')
  parser.add_argument('--directories', required=True, nargs='+', help='location of each batch')
  parser.add_argument('--phenotype', required=False, help='location of phenotype')
  parser.add_argument('--require', required=False, help='field to require')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.directories, args.phenotype, args.require)
