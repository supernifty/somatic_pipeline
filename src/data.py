#!/usr/bin/env python
# generates symlinks to data
# usage:
# python data.py < data.tsv

import glob
import os
import os.path
import sys

def run(cmd):
  sys.stderr.write('{}...\n'.format(cmd))
  result = os.system(cmd)
  if result != 0:
    sys.stderr.write('ERROR executing {}'.format(cmd))
    sys.exit(result)
  sys.stderr.write('{}: done\n'.format(cmd))

def main(source, dest):
  for line in sys.stdin:
    fields = line.strip('\n').split('\t') # A08978_25945    0151010101_BC
    for f in glob.glob('{source}/{prefix}_*'.format(source=source, prefix=fields[0])):
      target=os.path.basename(f).replace(fields[0], fields[1])
      run("ln -s {source_file} {dest}/{target}".format(source_file=f, dest=dest, target=target))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--source_dir', required=True, help='source of fastq')
  parser.add_argument('--target_dir', required=True, help='in directory')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
  main(source, dest)
