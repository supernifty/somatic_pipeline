#!/usr/bin/env python

import csv
import sys

THRESHOLD = 0.4

writer = csv.DictWriter(sys.stdout, delimiter='\t', fieldnames=['Sample', 'Class', 'Pct'])
writer.writeheader()
for filename in sys.argv[1:]:
  reader = csv.DictReader(open(filename, 'r'), delimiter='\t')
  for row in reader:
    if row['Locus'] == 'Average':
      result = {}
      result['Sample'] = filename.split('/')[-1].split('.')[0]
      result['Pct'] = row['Difference']
      if float(row['Difference']) > THRESHOLD:
        result['Class'] = 'MSI-H'
      else:
        result['Class'] = 'MSS'
      writer.writerow(result)

