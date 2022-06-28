#!/usr/bin/env python
'''
  calculate af for strelka
'''

import collections
import logging
import sys

import numpy

import cyvcf2

def main(sample, vcf_fn):
  '''
    refCounts = Value of FORMAT column $REF + "U" (e.g. if REF="A" then use the value in FOMRAT/AU)
    altCounts = Value of FORMAT column $ALT + "U" (e.g. if ALT="T" then use the value in FOMRAT/TU)
    tier1RefCounts = First comma-delimited value from $refCounts
    tier1AltCounts = First comma-delimited value from $altCounts
    Somatic allele freqeuncy is $tier1AltCounts / ($tier1AltCounts + $tier1RefCounts)
  '''

  logging.info('reading %s...', vcf_fn)

  vcf_in = cyvcf2.VCF(vcf_fn)  
  vcf_in.add_info_to_header({'ID': 'AF', 'Description': 'Calculated allele frequency', 'Type':'Float', 'Number': '1'})

  sys.stdout.write(vcf_in.raw_header)

  if sample in ('0', '1'):
    sample_id = int(sample)
  else:
    sample_id = vcf_in.samples.index(sample)

  variant_count = multi = 0
  seen = set()
  for variant_count, variant in enumerate(vcf_in):
    if 'AD' in vcf_in:
      #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  STE0072C12I5EWT5
      #1       10159   .       A       ACCG    46      LowGQX;NoPassedVariantGTs       CIGAR=1M3I;RU=CCG;REFREP=0;IDREP=1;MQ=33        GT:GQ:GQX:DPI:AD:ADF:ADR:FT:PL  0/1:88:0:856:372,80:249,74:123,6:LowGQX:85,0,999
      #logging.info(variant.format('AD'))
      try:
        # handle multiallelic (could also consider ignoring)
        loci = '{}/{}/{}'.format(variant.CHROM, variant.POS, variant.REF)
        if loci in seen:
          #logging.info('skipped multiallelic at %s', loci)
          multi += 1
          continue
        else:
          refCount = variant.format('AD')[0][0]
          altCount = variant.format('AD')[0][1]
          seen.add(loci) # to deal with multiallelic
      except:
        logging.warn('bad AD on line %i at %s %s>%s: %s', variant_count, variant.POS, variant.REF, variant.ALT, variant.format('AD'))
        raise
    else:
      # GL000220.1      135366  .       T       C       .       LowEVS;LowDepth SOMATIC;QSS=1;TQSS=1;NT=ref;QSS_NT=1;TQSS_NT=1;SGT=TT->TT;DP=2;MQ=60.00;MQ0=0;ReadPosRankSum=0.00;SNVSB=0.00;SomaticEVS=0.71    DP:FDP:SDP:SUBDP:AU:CU:GU:TU    1:0:0:0:0,0:0,0:0,0:1,1 1:0:0:0:0,0:1,1:0,0:0,0
      tier1RefCounts = variant.format('{}U'.format(variant.REF))
      tier1AltCounts = variant.format('{}U'.format(variant.ALT[0])) # assume not multiallelic
      if len(variant.ALT) > 1:
        logging.warn('%s: variant %i is multi-allelic', vcf_fn, variant_count + 1)

      # just tier 1
      #refCount = tier1RefCounts[sample_id][0]
      #altCount = tier1AltCounts[sample_id][0]

      # both tiers
      refCount = sum(tier1RefCounts[sample_id])
      altCount = sum(tier1AltCounts[sample_id])

    if refCount + altCount == 0:
      af = 0.0
    else:
      af = 1.0 * altCount / (altCount + refCount)

    variant.INFO["AF"] = af
    sys.stdout.write(str(variant))
    if (variant_count + 1 ) % 100000 == 0:
      logging.info('reading %s: %i variants processed, skipped %i...', vcf_fn, variant_count + 1, multi)

  logging.info('reading %s: processed %i variants, skipped %i multi-allelic', vcf_fn, variant_count + 1, multi)

if __name__ == '__main__':
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  main(sys.argv[1], sys.argv[2])
