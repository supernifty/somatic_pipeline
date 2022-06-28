#!/usr/bin/env python

import csv
import sys

first = True
for filename in sys.argv[1:]:
  sys.stderr.write('processing {}...\n'.format(filename))
  reader = csv.DictReader(open(filename, 'r'), delimiter='\t')
  if first:
    writer = csv.DictWriter(sys.stdout, delimiter='\t', fieldnames=['Sample'] + reader.fieldnames)
    writer.writeheader()
    first = False

  for row in reader:
      row['Sample'] = filename.split('/')[-1].split('.')[0]
      #result['Sample'] = filename
      writer.writerow(row)
      # just write first row
      break

