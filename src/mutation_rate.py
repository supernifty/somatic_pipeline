#!/usr/bin/env python

'''
  Intersects a vcf file with a bed and annotates it with bed values

  Usage:
    $0 bed < vcf > annotated_filtered_vcf
'''

import argparse
import logging
import sys

import cyvcf2
import intervaltree

def main(vcfs, bed, min_dp, min_af, min_qual, indels):
  logging.info('parsing {}...'.format(bed))
  tree = {}
  size = overlaps = skipped = included = 0
  for line_count, line in enumerate(open(bed, 'r')):
    fields = line.strip('\n').split('\t')
    if len(fields) < 4:
      skipped += 1
      continue
    chrom, start, finish, annotation = fields[:4]
    if chrom.startswith('chr'):
      chrom = chrom[3:]
    if chrom not in tree:
      tree[chrom] = intervaltree.IntervalTree()
    s = int(start)
    f = int(finish)
    overlap = tree[chrom].search(s, f)
    if len(overlap) == 0:
      size += f - s
      tree[chrom][s:f] = True
      included += 1
    else:
      for item in overlap:
        new_begin = min(s, item.begin)
        new_end = max(f, item.end)
        tree[chrom].remove(item)
        tree[chrom][new_begin:new_end] = True
        overlaps += 1
        size += (new_end - new_begin) - (f - s)
        s = new_begin
        f = new_end
    if line_count % 100000 == 0:
      logging.debug('parsing {}: {} lines parsed. skipped {}. {} overlaps. size {}'.format(bed, line_count, skipped, overlaps, size))
  logging.info('parsing {}: done. lines skipped: {}. size: {}. count: {}'.format(bed, skipped, size, included))

  sys.stdout.write("Filename\tCount\tPerMB\tPerInterval\n")
  for vcf in vcfs:
    vcf_in = cyvcf2.VCF(vcf)
    accept = reject_exon = reject_filter = reject_indel = 0
    for variant in vcf_in:
      if variant.QUAL is not None and variant.QUAL < min_qual:
        reject_filter += 1
        continue
      #if variant.FILTER is not None:
      #  reject_filter += 1
      #  continue
      #dp = variant.format('DP')[0][0]
      #if dp < min_dp:
      #  reject_filter += 1
      #  continue
      #ad = variant.format('AD')[0][0]
      #af = ad / dp
      #if af < min_af:
      #  reject_filter += 1
      #  continue
      if indels and len(variant.REF) != len(variant.ALT[0]):
        reject_indel += 1
        continue

      if variant.CHROM.startswith('chr'):
        chrom = variant.CHROM[3:]
      else:
        chrom = variant.CHROM
      if chrom in tree:
        overlap = tree[chrom].search(variant.POS)
        if len(overlap) == 0:
          reject_exon += 1
        else:
          accept += 1
    sys.stdout.write("{}\t{}\t{:.2f}\t{:.2f}\n".format(vcf, accept, accept / size * 1e6, accept / included))
    logging.info('included %i variants. rejected %i non-exonic %i filtered %i non-indel variants.', accept, reject_exon, reject_filter, reject_indel)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Calculate mutation rate')
  parser.add_argument('--vcfs', nargs='+', help='list of vcfs')
  parser.add_argument('--bed', required=True, help='filter')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--indels_only', action='store_true', help='just indels')
  parser.add_argument('--min_dp', default=0, type=int, help='min dp')
  parser.add_argument('--min_af', default=0, type=float, help='min dp')
  parser.add_argument('--min_qual', default=0, type=float, help='min qual')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.vcfs, args.bed, args.min_dp, args.min_af, args.min_qual, args.indels_only)
