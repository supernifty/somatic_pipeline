#!/usr/bin/env python
'''
  given tumour and normal vcf pairs, explore msi status
'''

import argparse
import collections
import logging
import sys

import numpy as np

import cyvcf2

ORDER=('HIGH_transcript_ablation', 'HIGH_splice_acceptor_variant', 'HIGH_splice_donor_variant', 'HIGH_stop_gained', 'HIGH_frameshift_variant', 'HIGH_stop_lost', 'HIGH_start_lost', 'HIGH_transcript_amplification', 'MODERATE_inframe_insertion', 'MODERATE_inframe_deletion', 'MODERATE_missense_variant', 'MODERATE_protein_altering_variant')

#IDS=['0'..'9'] + ['a'..'z'] + ['A'..'Z']

def main(vcfs, lohs, genes, loci, plot):
  logging.info('starting...')

  results = collections.defaultdict(dict)

  variants = {}
  count = 0

  for vcf in vcfs:
    logging.info('processing %s...', vcf)
    sample = vcf.split('/')[-1].split('.')[0]
    result = results[sample] # dict
    
    vcf_in = cyvcf2.VCF(vcf)
    for variant in vcf_in:
      # assume vep format is Consequence|IMPACT|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|HGVSc|HGVSp|cDNA_position|CDS_position|HGVSc|HGVSp|cDNA_position|CDS_position|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MaxEntScan_alt|MaxEntScan_diff|MaxEntScan_ref|PICK
      veps = variant.INFO['CSQ'].split(',')
      # find PICK=1
      for vep in veps:
        fields = vep.split('|')
        if fields[-1] != '1':
          continue
        gene = fields[5]
        if gene not in genes:
          continue

        impact = fields[1]
        if impact not in ('MODERATE', 'HIGH'):
          continue

        consequence = fields[0].split('&')[0]
        #candidate = '{}_{}'.format(impact, consequence)
        candidate = consequence
        polyphen = fields[8]
        sift = fields[9]
        #if candidate not in ORDER:
        #  logging.warn('unexpected candidate: %s', candidate)
        if impact == 'HIGH' or 'damaging' in polyphen or 'deleterious' in sift: # include
          if gene not in result:
            result[gene] = {}
          if candidate not in result[gene]:
            result[gene][candidate] = list()
          change = '{}:{} {}>{}'.format(variant.CHROM, variant.POS, variant.REF, variant.ALT[0])
          result[gene][candidate].append(change)
          # keep track of found variants
          if change not in variants:
            variants[change] = '{}\t{}\t{}'.format(gene, fields[12], fields[13])

  for loh in lohs:
    logging.info('processing %s...', loh)
    sample = loh.split('/')[-1].split('.')[0]
    result = results[sample] # dict
    for line in open(loh, 'r'):
      fields = line.strip('\n').split('\t')
      chrom, start, finish = fields[0], int(fields[1]), int(fields[2]) # loh bed
      # find matching loci
      for idx, locus in enumerate(loci): # check each gene
        locus_chrom, rest = locus.split(':')
        locus_start, locus_finish = [int(x) for x in rest.split('-')]
        if locus_chrom == chrom and finish > locus_start and start < locus_finish:
          if genes[idx] not in result:
            result[genes[idx]] = {}
          if 'LOH' not in result[genes[idx]]:
            result[genes[idx]]['LOH'] = list()
          results[sample][genes[idx]]['LOH'].append('1')

  # print results
  sys.stdout.write('Sample\t{}\n'.format('\t'.join(sorted(genes))))
  for sample in sorted(results.keys()):
    output = []
    for gene in sorted(genes):
      if gene in results[sample]:
        cell = []
        for key in results[sample][gene]:
          #cell.append('{} {}'.format(key, len(results[sample][gene][key])))
          cell.append('{} {}'.format(key, ' '.join(results[sample][gene][key])))
        output.append(','.join(cell))
      else:
        output.append('None')
    sys.stdout.write('{}\t{}\n'.format(sample, '\t'.join(output)))

  # now print variants
  sys.stdout.write('\nVariants\tGene\tc-change\tp-change\n')
  for variant in sorted(variants.keys()):
    sys.stdout.write('{}\t{}\n'.format(variant, variants[variant]))

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--vcfs', required=True, nargs='+', help='tumour vcf')
  parser.add_argument('--lohs', required=True, nargs='+', help='tumour vcf')
  parser.add_argument('--genes', required=True, nargs='+', help='tumour vcf')
  parser.add_argument('--loci', required=True, nargs='+', help='tumour vcf')
  parser.add_argument('--plot', required=False, help='tumour vcf')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.vcfs, args.lohs, args.genes, args.loci, args.plot)

