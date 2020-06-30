#!/usr/bin/env python
'''
  implement msiseq

  to use with hg38:
  wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz
  zcat simpleRepeat.txt.gz | cut -f2,3,4 > reference/hg38.repeats.bed
'''

import argparse
import collections
import csv
import gzip
import logging
import sys

import intervaltree

import cyvcf2

def add_intersect(tree, chrom, s, f):
  tree[chrom][s:f] = True

def bed_to_tree(bed):
  logging.info('parsing {}...'.format(bed))
  tree = {}
  size = overlaps = skipped = included = 0
  for line_count, line in enumerate(open(bed, 'r')):
    fields = line.strip('\n').split('\t')
    if len(fields) < 3:
      skipped += 1
      continue
    chrom, start, finish = fields[:3]
    if chrom.startswith('chr'):
      chrom = chrom[3:]
    if chrom not in tree:
      tree[chrom] = intervaltree.IntervalTree()
      logging.info('processing %s...', chrom)
    s = int(start)
    f = int(finish)
    # does it overlap with itself?
    overlap = tree[chrom].search(s, f)
    if len(overlap) == 0:
      size += f - s
      add_intersect(tree, chrom, s, f) 
      included += 1
    else:
      if len(overlap) > 1:
        logging.warn('unexpected multiple overlaps: %s', overlap)
      for item in overlap:
        #logging.debug('overlap %s:%i-%i with %s:%i-%i', chrom, item.begin, item.end, chrom, s, f)
        new_begin = min(s, item.begin)
        new_end = max(f, item.end)
        tree[chrom].remove(item)
        overlaps += 1

      add_intersect(tree, chrom, new_begin, new_end)
      size += (new_end - new_begin) - (f - s)
      s = new_begin
      f = new_end

    if line_count % 100000 == 0:
      logging.debug('parsing {}: {} lines parsed. skipped {}. {} overlaps. size {}'.format(bed, line_count, skipped, overlaps, size))
  logging.info('parsing {}: done. lines skipped: {}. size: {}. count: {}'.format(bed, skipped, size, included))
  return tree, size

def is_indel(v):
  return len(v.REF.replace('-', '')) != len(v.ALT[0].replace('-', ''))

def vcf_list(vcfs, is_maf):
  if is_maf:
    for maf in vcfs:
      logging.info('processing %s...', maf)
      for sample in maf_samples(maf):
        logging.info('processing %s with sample %s...', maf, sample)
        yield (sample, maf_to_vcf(maf, sample))
  else:
    for vcf in vcfs:
      logging.info('processing %s...', vcf)
      yield (vcf, cyvcf2.VCF(vcf))

def msiseq(vcfs, repeats, capture, threshold, capture_size, is_maf):
  capture_tree = None
  if capture is not None:
    capture_tree, capture_size = bed_to_tree(capture)
  repeats_tree, repeats_size = bed_to_tree(repeats)

  sys.stdout.write('Sample\tS.ind.count\tS.ind\tT.ind\tClass\n')
  for vcf, vcf_in in vcf_list(vcfs, is_maf):
    count = reject_capture = reject_repeat = 0
    for v in vcf_in:
      if is_indel(v): 
        if v.CHROM.startswith('chr'):
          chrom = v.CHROM[3:]
        else:
          chrom = v.CHROM
 
        if capture_tree is None or chrom in capture_tree and len(capture_tree[chrom].search(v.POS)) > 0:
          if chrom in repeats_tree and len(repeats_tree[chrom].search(v.POS)) > 0:
            # it's a simple repeat
            count += 1
          else:
            # it's an indel in the capture 
            reject_repeat += 1
        else:
          # it's an indel not in the capture
          reject_capture += 1

    logging.info('processing %s: count %i outside capture: %i in capture but not in a repeat: %i', vcf, count, reject_capture, reject_repeat)
    per_mb = count / capture_size * 1000000
    if per_mb > threshold:
      classification = 'MSI-H'
    else:
      classification = 'MSS'
    t_ind = (count + reject_repeat) / capture_size * 1000000
    sys.stdout.write('{}\t{}\t{:.3f}\t{:.3}\t{}\n'.format(vcf.split('/')[-1].split('.')[0], count, per_mb, t_ind, classification))
  logging.info('done')

def get_value(header, col, row):
  return row[header.index(col)]

def open_file(fn, is_gzipped):
  if is_gzipped:
    return gzip.open(fn, 'rt')
  else:
    return open(fn, 'rt')

def maf_samples(maf, sample_col='Tumor_Sample_Barcode'):
  samples = set()
  header = None
  for line, row in enumerate(csv.reader(open_file(maf, True), delimiter='\t')):
    if row[0].startswith('#'):
      continue
    if header is None:
      header = row
      continue
    row_sample = get_value(header, sample_col, row)
    samples.add(row_sample)

  return samples    

def maf_to_vcf(maf, sample, sample_col='Tumor_Sample_Barcode', chrom_col='Chromosome', pos_col='Start_Position', ref_col='Reference_Allele', alt_col='Tumor_Seq_Allele2'):

  Variant = collections.namedtuple('Variant', 'CHROM POS REF ALT')

  # enumeration a maf into a variant
  header = None
  for line, row in enumerate(csv.reader(open_file(maf, True), delimiter='\t')):
    if line % 1000 == 0:
      logging.debug('processed %i lines of %s...', line, maf)

    if row[0].startswith('#'):
      continue
    if header is None:
      header = row
      continue

    #Hugo_Symbol     Entrez_Gene_Id  Center  NCBI_Build      Chromosome      Start_Position  End_Position    Strand  Variant_Classification  Variant_Type    Reference_Allele        Tumor_Seq_Allele1       Tumor_Seq_Allele2       dbSNP_RS        dbSNP_Val_Status        Tumor_Sample_Barcode    Matched_Norm_Sample_Barcode     Match_Norm_Seq_Allele1  Match_Norm_Seq_Allele2  Tumor_Validation_Allele1        Tumor_Validation_Allele2        Match_Norm_Validation_Allele1   Match_Norm_Validation_Allele2   Verification_Status     Validation_Status       Mutation_Status Sequencing_Phase        Sequence_Source Validation_Method       Score   BAM_File        Sequencer       Tumor_Sample_UUID       Matched_Norm_Sample_UUID        HGVSc   HGVSp   HGVSp_Short     Transcript_ID   Exon_Number     t_depth t_ref_count     t_alt_count     n_depth n_ref_count     n_alt_count     all_effects     Allele  Gene    Feature Feature_type    One_Consequence Consequence     cDNA_position   CDS_position    Protein_position        Amino_acids Codons  Existing_variation      ALLELE_NUM      DISTANCE        TRANSCRIPT_STRAND       SYMBOL  SYMBOL_SOURCE   HGNC_ID BIOTYPE CANONICAL       CCDS    ENSP    SWISSPROT       TREMBL  UNIPARC RefSeq  SIFT    PolyPhen        EXON    INTRON  DOMAINS GMAF    AFR_MAF AMR_MAF ASN_MAF EAS_MAF EUR_MAF SAS_MAF AA_MAF  EA_MAF  CLIN_SIG        SOMATIC PUBMED  MOTIF_NAME      MOTIF_POS       HIGH_INF_POS    MOTIF_SCORE_CHANGE      IMPACT  PICK    VARIANT_CLASS   TSL     HGVS_OFFSET     PHENO   MINIMISED       ExAC_AF ExAC_AF_Adj     ExAC_AF_AFR     ExAC_AF_AMR     ExAC_AF_EAS     ExAC_AF_FIN     ExAC_AF_NFE     ExAC_AF_OTH     ExAC_AF_SAS     GENE_PHENO      FILTER  CONTEXT src_vcf_id      tumor_bam_uuid  normal_bam_uuid case_id GDC_FILTER      COSMIC  MC3_Overlap     GDC_Validation_Status

    row_sample = get_value(header, sample_col, row)
    if sample is not None and row_sample != sample:
      continue

    chrom = get_value(header, chrom_col, row).replace('chr', '')
    pos = int(get_value(header, pos_col, row))
    ref = get_value(header, ref_col, row)
    if ref == '-':
      pos += 1 # fix for TCGA mafs
    ref = ref.replace('-', '')
    alt = get_value(header, alt_col, row).replace('-', '')

    yield Variant(chrom, pos, ref, (alt,))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Classify MSI using msiseq algorithm')
  parser.add_argument('--vcfs', required=True, nargs='+', help='input vcf files')
  parser.add_argument('--is_maf', action='store_true', help='vcf is actually a maf')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--repeats', required=True, help='bed file of repeats')
  parser.add_argument('--capture', required=False, help='bed file of capture')
  parser.add_argument('--capture_size', required=False, type=int, help='capture size')
  parser.add_argument('--threshold', required=False, default=0.395, type=float, help='cutoff for msi-h')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
  msiseq(args.vcfs, args.repeats, args.capture, args.threshold, args.capture_size, args.is_maf)

