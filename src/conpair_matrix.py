#!/usr/bin/env python
'''
  run conpair on all combinations
'''

import argparse
import logging
import os
import sys

def execute(cmd):
  logging.debug('executing %s...', cmd)
  result = os.system(cmd)
  if result != 0:
    logging.warn('executing %s FAILED', cmd)
  logging.debug('executing %s: done', cmd)
  return result == 0

def main(tumours, normals, working, output, reference, skip_pileup, skip_comparison):

  if not skip_pileup:
    logging.info('generating pileups...')
    samples = set(tumours + normals)
    for s in samples:
      logging.info('generating pileup for %s...', s)
      name = s.split('/')[-1].split('.')[0]
      if not execute("python tools/Conpair/scripts/run_gatk_pileup_for_sample.py --reference {reference} --conpair_dir tools/Conpair --gatk tools/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -B {bam} -O {working}/{name}.pileup".format(reference=reference, bam=s, working=working, name=name)):
        sys.exit(1)

  if not skip_comparison:
    logging.info('running comparisons...')
    for t in tumours:
      for n in normals:
        if t == n:
          continue
        logging.info('comparing %s to %s...', t, n)
        tname = t.split('/')[-1].split('.')[0]
        nname = n.split('/')[-1].split('.')[0]
        if not execute("PYTHONPATH=tools/Conpair/modules CONPAIR_DIR=tools/Conpair python tools/Conpair/scripts/verify_concordance.py -T {working}/{tname}.pileup -N {working}/{nname}.pileup --outfile {working}/{tname}.{nname}.concordance --normal_homozygous_markers_only".format(working=working, nname=nname, tname=tname)):
          sys.exit(1)

  logging.info('merging results...')
  with open(output, 'w') as fh:
    fh.write('Tumour\tNormal\tConcordance\n')
    for t in tumours:
      for n in normals:
        if t == n:
          continue
        tname = t.split('/')[-1].split('.')[0]
        nname = n.split('/')[-1].split('.')[0]
        for line in open('{working}/{tname}.{nname}.concordance'.format(working=working, tname=tname, nname=nname)).readlines():
          if line.startswith('Concordance:'):
            fh.write('{}\t{}\t{}\n'.format(tname, nname, line.strip('\n').split(' ')[1].replace('%', '')))
            break

  #logging.info('plotting results...')
  #plotme.heatmap.plot_heat(open('{}.png'.format(output), 'w'), 
  logging.info('done')

#    "python src/conpair_matrix.py --tumours {input.tumours} --germlines {input.germlines} --reference {input.reference} --output tmp/conpair > {output} 2>{log.stderr}"
if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--tumours', required=True, nargs='+', help='tumour bams')
  parser.add_argument('--germlines', required=True, nargs='+', help='germline bams')
  parser.add_argument('--working', required=True, help='temp directory')
  parser.add_argument('--output', required=True, help='output file')
  parser.add_argument('--reference', required=True, help='temp directory')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--skip_pileup', action='store_true', help='pileup has already been run')
  parser.add_argument('--skip_comparison', action='store_true', help='comparison has already been run')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.tumours, args.germlines, args.working, args.output, args.reference, args.skip_pileup, args.skip_comparison)
