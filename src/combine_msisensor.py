#!/usr/bin/env python

import csv
import sys

# based on https://github.com/ding-lab/msisensor/issues/11 and their paper
CLASSES = (
  (0, 1, 'MSI-L'),
  (1, 3.5, 'MSS'),
  (3.5, 100, 'MSI-H')
)

writer = csv.DictWriter(sys.stdout, delimiter='\t', fieldnames=['Sample', 'Class', '%', 'Total_Number_of_Sites', 'Number_of_Somatic_Sites' ])
writer.writeheader()
for filename in sys.argv[1:]:
  reader = csv.DictReader(open(filename, 'r'), delimiter='\t')
  for row in reader:
    classification = ''
    val = float(row['%'])
    for c in CLASSES:
      if c[0] <= val < c[1]:
        classification = c[2]

    row['Class'] = classification
    row['Sample'] = filename.split('/')[-1].split('.')[0]
    writer.writerow(row)
