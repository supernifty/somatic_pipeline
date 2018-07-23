#!/usr/bin/env python
'''
  plot bed coverage
'''

import argparse
import logging
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import matplotlib.cm as cm

def main(target, bedfiles, max_cov):
  logging.info('starting...')

  fig=plt.figure(figsize=(16,12))
  fig.show()
  ax=fig.add_subplot(111)
  ax.set_xscale("log", nonposx='clip')
  ax.set_xlabel("Coverage")
  ax.set_ylabel("Proportion of bases with at least this coverage")
  ax.set_title("Target Region Coverage")
  ax.grid()

  colors = cm.rainbow(np.linspace(0, 1, len(bedfiles)))
  markers = (".", ",", "o", "v", "^", "<", ">", "1", "2", "3", "4", "8", "s", "p", "P", "*", "h", "H", "+", "x", "X", "D", "d", "|", "_")
  
  for i, bedfile in enumerate(bedfiles):
    xs = []
    ys = []
    cumulative = 1.0
    for line in open(bedfile, 'r'):
      fields = line.strip('\n').split('\t')
      x = float(fields[1])
      if x > max_cov:
        break
      xs.append(x)
      cumulative -= float(fields[4])
      ys.append(cumulative)

    ax.plot(xs, ys, color=colors[i], marker=markers[i % len(markers)], label=bedfile.split('.')[0], markevery=0.1)

  plt.legend()
  plt.savefig(target)
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--target', required=True, help='more logging')
  parser.add_argument('--max', required=False, type=int, default=10000, help='more logging')
  parser.add_argument('--files', required=True, nargs='+', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.target, args.files, args.max)
