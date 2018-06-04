#!/usr/bin/env python
# generates symlinks to data
# usage:
# python data.py < data.tsv

import glob
import os
import os.path
import sys

SOURCE="/scratch/UOM0040/Data/AGRF_CAGRF17380_H5CNYDSXX"
DEST="/scratch/UOM0040/peter/AGRF_CAGRF17380_H5CNYDSXX/in"

def run(cmd):
  sys.stderr.write('{}...\n'.format(cmd))
  result = os.system(cmd)
  if result != 0:
    sys.stderr.write('ERROR executing {}'.format(cmd))
    sys.exit(result)
  sys.stderr.write('{}: done\n'.format(cmd))

def main():
  for line in sys.stdin:
    fields = line.strip('\n').split('\t') # A08978_25945    0151010101_BC
    for f in glob.glob('{source}/{prefix}_*'.format(source=SOURCE, prefix=fields[0])):
      target=os.path.basename(f).replace(fields[0], fields[1])
      run("ln -s {source_file} {dest}/{target}".format(source_file=f, dest=DEST, target=target))

if __name__ == '__main__':
  main()
