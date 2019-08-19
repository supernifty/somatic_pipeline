#!/usr/bin/env python
'''
  given tumour and normal vcf pairs, explore msi status
'''

import argparse
import logging
import sys

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FIGSIZE=(18,18)
#CLASSES=['None', 'LOH', 'splice_acceptor_variant', 'missense_variant', 'stop_gained', 'frameshift_variant'] # 0..5
CLASSES=['frameshift_variant', 'stop_gained', 'missense_variant', 'splice_donor_variant', 'splice_acceptor_variant', 'LOH', 'None'] # 0..6
COLORS=['#ff3030', '#ff60ff', '#f0d070', '#7070ff', '#9050ff', '#70f070', '#e0e0e0']

## matplotlib helpers from https://matplotlib.org/gallery/images_contours_and_fields/image_annotated_heatmap.html
def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Arguments:
        data       : A 2D numpy array of shape (N,M)
        row_labels : A list or array of length N with the labels
                     for the rows
        col_labels : A list or array of length M with the labels
                     for the columns
    Optional arguments:
        ax         : A matplotlib.axes.Axes instance to which the heatmap
                     is plotted. If not provided, use current axes or
                     create a new one.
        cbar_kw    : A dictionary with arguments to
                     :meth:`matplotlib.Figure.colorbar`.
        cbarlabel  : The label for the colorbar
    All other arguments are directly passed on to the imshow call.
    """

    if not ax:
        ax = plt.gca()

    logging.info('%i rows %i columns', len(row_labels), len(col_labels))
    logging.info('%s', data)

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar

def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=["black", "white"],
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Arguments:
        im         : The AxesImage to be labeled.
    Optional arguments:
        data       : Data used to annotate. If None, the image's data is used.
        valfmt     : The format of the annotations inside the heatmap.
                     This should either use the string format method, e.g.
                     "$ {x:.2f}", or be a :class:`matplotlib.ticker.Formatter`.
        textcolors : A list or array of two color specifications. The first is
                     used for values below a threshold, the second for those
                     above.
        threshold  : Value in data units according to which the colors from
                     textcolors are applied. If None (the default) uses the
                     middle of the colormap as separation.

    Further arguments are passed on to the created text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[im.norm(data[i, j]) > threshold])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts


def main(target):
  logging.info('reading from stdin...')
  header = sys.stdin.readline().strip('\r\n').split('\t')
  gene_names = header[1:]
  sample_names = []
  results = []
  for line in sys.stdin:
    fields = line.strip('\r\n').split('\t')
    if len(fields) <= 1:
      break
    sample_names.append(fields[0])
    results.append([CLASSES.index(x.split(' ')[0]) for x in fields[1:]]) # take the first word

  fig = plt.figure(figsize=(FIGSIZE[0], 1 + len(sample_names) / 2))
  ax = fig.add_subplot(111)
  #data = np.empty([len(sample_names), len(gene_names)])
  data = np.array(results)

  y = sample_names # each sample
  x = gene_names # each gene name

  norm = matplotlib.colors.BoundaryNorm(np.linspace(-0.5, len(CLASSES) - 0.5, len(CLASSES) + 1), len(CLASSES))
  fmt = matplotlib.ticker.FuncFormatter(lambda x, pos: CLASSES[norm(x)])
  #fmt = matplotlib.ticker.FuncFormatter(lambda x, pos: categories[::-1][norm(x)])

  # from most serious to least serious
  cmap = matplotlib.colors.ListedColormap(COLORS)

  im, _ = heatmap(data, y, x, ax=ax, cmap=cmap, norm=norm, cbar_kw=dict(ticks=np.arange(0, len(CLASSES)), format=fmt))

  plt.savefig(args.target)
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('--target', required=True, help='filename to generate')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.target)

