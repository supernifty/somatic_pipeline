#!/usr/bin/env python

import sys

def process(in_fh, out_fh):

  # write header
  out_fh.write('<html><head>' +
    '<link rel="stylesheet" href="https://unpkg.com/purecss@1.0.0/build/pure-min.css" integrity="sha384-nn4HPE8lTHyVtfCBi5yW9d20FjT8BJwUXyWZT9InLYax14RDjBj46LmSztkmNP9w" crossorigin="anonymous">' +
    '</head><body>\n')

  for line in in_fh:
    styled = line.replace('<table', '<table class="pure-table"')
    out_fh.write(styled)

  # write footer
  out_fh.write('\n</body></html>\n')

if __name__ == '__main__':
  process(sys.stdin, sys.stdout)
