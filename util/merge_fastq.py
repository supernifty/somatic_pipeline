#!/usr/bin/env python
'''
  given tumour and normal fastqs, generate a script to merge them
'''

# ./Exome_fastq/1053-D-T-1_D13N8ACXX_AGTCAA_L002_R1.fastq.gz
#./Exome_fastq/1053-D-T-1_DOWTAACXX_AGTCAA_L003_R2.fastq.gz
#./Exome_fastq/1053-D-T-1_D13N8ACXX_AGTCAA_L002_R2.fastq.gz
#./Exome_fastq/1053-D-T-1_DOWTAACXX_AGTCAA_L004_R2.fastq.gz
#./Exome_fastq/1053-D-T-1_DOWTAACXX_AGTCAA_L003_R1.fastq.gz
#./Exome_fastq/1053-D-T-1_DOWTAACXX_AGTCAA_L004_R1.fastq.gz
#./fastq_normal/1053-D-N-1_DOWTAACXX_AGTTCC_L003_R1.fastq.gz
#./fastq_normal/1053-D-N-1_D13N8ACXX_AGTTCC_L002_R2.fastq.gz
#./fastq_normal/1053-D-N-1_DOWTAACXX_AGTTCC_L004_R2.fastq.gz
#./fastq_normal/1053-D-N-1_DOWTAACXX_AGTTCC_L003_R2.fastq.gz
#./fastq_normal/1053-D-N-1_D13N8ACXX_AGTTCC_L002_R1.fastq.gz
#./fastq_normal/1053-D-N-1_DOWTAACXX_AGTTCC_L004_R1.fastq.gz
# 9910279001_BC_HGNWHDSXY_CACTTCGA_L001_I2.fastq.gz

import argparse
import collections
import logging
import sys

def merge(files, outdir, sample_components, skip_components):
  logging.info('starting...')
  samples = collections.defaultdict(list)
  logging.info('considering input files...')
  for f in files:
    # extract sample and read
    fn = f.split('/')[-1].split('.')[0]
    components = fn.split('_')
    sample = '_'.join(components[skip_components:sample_components])
    readnum = components[-1]
    samples[(sample, readnum)].append(f)
    logging.debug('added %s to %s...', f, (sample, readnum))
  logging.info('generating merged files...')
  sys.stdout.write('#!/usr/bin/env bash\n# generated for {} for {}\n\n'.format(', '.join(files)[:1000], outdir))
  sys.stdout.write('echo "starting at $(date)"\n\n')
  sys.stdout.write('set -o errexit\nset -x\n')
  for s in samples:
    sample, readnum = s
    if len(samples[s]) == 1:
      # symlink
      sys.stdout.write('ln -s "{}" "{}"\n'.format(samples[s][0], outdir))
    else:
      # merge
      target = sorted([x[::-1] for x in samples[s]])[0] # take sample with lowest lane number
      target = target[::-1].split('/')[-1]
      sys.stdout.write('cat {} > {}/{}\n'.format(' '.join(samples[s]), outdir, target))

  sys.stdout.write('echo "finished at $(date)"\n')
 
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Merge multiple fastq lanes')
  parser.add_argument('--files', required=True, nargs='+', help='fastq input files')
  parser.add_argument('--outdir', required=True, help='fastq target')
  parser.add_argument('--skip_components', required=False, type=int, default=0, help='skip initial components')
  parser.add_argument('--sample_components', required=False, type=int, default=1, help='number of components in sample')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  merge(args.files, args.outdir, args.sample_components, args.skip_components)
