#!/usr/bin/env python

import csv
import logging
import sys

delimiter = '\t'

fh=csv.DictReader(sys.stdin, delimiter=delimiter)
out = csv.DictWriter(sys.stdout, delimiter=delimiter, fieldnames=fh.fieldnames + ['VAF'])
out.writeheader()

for idx, row in enumerate(fh):
  if 'AD' not in row or ',' not in row['AD']:
    row['VAF'] = 0
  else:
  #sys.stderr.write('processing {}\n'.format(row))
    ref, alt = [int(x) for x in row['AD'].split(',')]
    if ref == 0 and alt == 0:
      row['VAF'] = 0
    else:
      row['VAF'] = '{:.2f}'.format(alt / (ref + alt))
  out.writerow(row)

