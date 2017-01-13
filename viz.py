#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')
import pandas
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1 import ImageGrid
import yaml
import pandas
import os

try:
    files = snakemake.input
    wc = snakemake.wildcards
    config = yaml.load(open('config.yaml'))

except NameError:
    wc = dict(
        state=4,
        subset='subset2',
        cell='umel',)

    files = dict(
        emissions='models/{subset}/{state}/emissions_{state}.txt'.format(**wc),
        transitions='models/{subset}/{state}/transitions_{state}.txt'.format(**wc),
        enrichment='models/{subset}/{state}/{cell}_enrichment.txt'.format(**wc),
        uniform_enrichment='models/{subset}/{state}/{cell}_enrichment.txt'.format(**wc),
        comparison='models/{subset}/{state}/{cell}_model_comparison_{state}.txt'.format(**wc),
    )


# load the subset sampletable to determine a uniform order. This could
# be forced by using the -holdcolumnorder param to LearnModel, but this way is
# a little more flexible.
subset = pandas.read_table(os.path.join('config', wc['subset'] + '.tsv'))
order = subset['mark'].values

dfs = {}
for k, v in files.items():
    df = pandas.read_table(v, index_col=0)
    if k == 'comparison':
        df[df < 0] = np.nan
    if 'enrichment' in k:
        df = df.fillna(0)
    dfs[k] = df

fig = plt.figure(figsize=(12, 6))
grid = ImageGrid(fig, 111, nrows_ncols=(1, 5),
                 axes_pad=0.1,
                 add_all=True,
                 label_mode="R")


cmaps = {
    'emissions': cm.Blues,
    'enrichment': cm.Greys,
    'uniform_enrichment': cm.Greys,
    'transitions': cm.Blues,
    'comparison': None,
}

for i, k in enumerate([
    'emissions',
    'enrichment',
    'uniform_enrichment',
    'transitions',
    'comparison'
]):
    ax = grid[i]
    v = dfs[k].copy()
    v.columns = [i.replace('.bed.gz', '').replace('.bed', '').replace('.txt', '') for i in v.columns]
    if k == 'uniform_enrichment':
        v = v.drop('Base')
        v = v.drop('Genome %', axis=1)
    if k == 'enrichment':
        v = v.drop('Base')
        v = v / v.max()
    if k == 'emissions':
        v = v.ix[:, order]
    sns.heatmap(v, ax=ax, cbar=False, cmap=cmaps[k], linewidths=0.5)
    ax.set_title(k)
    ax.patch.set_edgecolor('0.5')
    ax.patch.set_linewidth(2)
    if i > 0:
        for txt in ax.get_yticklabels():
            txt.set_visible(False)
        ax.yaxis.get_label().set_visible(False)
    else:
        ax.set_ylabel('State')
    for txt in ax.get_xticklabels():
        txt.set_rotation(90)

fig.tight_layout()
fig.suptitle('{subset}, {cell}, {state}-state'.format(**wc), weight='bold', size=20)
try:
    fig.savefig(snakemake.output.png)
    fig.savefig(snakemake.output.pdf)
except NameError:
    plt.show()
