#!/usr/bin/env python

import sys

stats = { 'n': 0, 't': 0, 'max': -1e9, 'min': 1e9 }
for line in sys.stdin:
  v = float(line.strip('\n'))
  stats['n'] += 1
  stats['t'] += v
  stats['max'] = max(stats['max'], v)
  stats['min'] = min(stats['min'], v)
  if stats['n'] % 1000000 == 0:
    stats['mean'] = stats['t'] / stats['n']
    sys.stderr.write('processed {n}: min {min} mean {mean:.3f} max {max}\n'.format(**stats))

stats['mean'] = stats['t'] / stats['n']
sys.stdout.write('n\tMean\tMax\tMin\tTotal\n')
sys.stdout.write('{n}\t{mean:.3f}\t{max}\t{min}\t{t}\n'.format(**stats))
