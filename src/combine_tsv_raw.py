#!/usr/bin/env python

import collections
import sys

def process():
  keys = set()
  vals = collections.defaultdict(dict)
  first_file = True
  for filename in sys.argv[1:]:
    first_line = True
    for line in open(filename, 'r'):
      # print header only if first file
      if first_file:
        first_file = False
        first_line = False
        sys.stdout.write('Filename\t{}'.format(line))
        continue
      elif first_line: # not first file
        first_line = False
        continue
      sys.stdout.write('{}\t{}'.format(filename, line))

if __name__ == '__main__':
  process()

