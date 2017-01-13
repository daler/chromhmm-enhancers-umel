#!/usr/bin/env python

import os
import sys
import argparse
from textwrap import dedent
import yaml
import matplotlib as mpl
import pybedtools
from pybedtools import featurefuncs as ff

config_example = dedent(
    """
    # After chromhmm.snakefile has been run, view the output and make
    # a decision on on what model to use. Enter that info here:
    final_model:
        subset: subset1
        states: 6
        celltype: umel

    # Decide on the labels and colors you want to use. "E1" etc are the
    # existing labels and correspond to the emission heatmaps for the model.

    label_mapping:
        # existing label    new label       color (hex)
        # ---------------- ---------------- -----------
        - ['E1',           'Active_TSS',    '#5AA02C']
        - ['E2',           'CTCF_open',     '#4A80C8']
        - ['E3',           'Enhancer',      '#FA9B00']
        - ['E4',           'weak_enhancer', '#FAF400']
        - ['E5',           'no_signal',     '#808080']
        - ['E6',           'repressed',     '#323232']
    """)

usage = dedent(
    """
    Colorizes and re-labels output from ChromHMM.

    Colors and labels are configured via a YAML-format file. A documented
    example can be generated with the --example-config option.
    """)

ap = argparse.ArgumentParser(usage=usage)
ap.add_argument('--config', help='Config file to use')
ap.add_argument('--example-config', action='store_true',
                help='Print an example config YAML and exit')
args = ap.parse_args()

if args.config is None and args.example_config is None:
    ap.print_help()
    sys.exit(0)

if args.example_config:
    print(config_example)
    sys.exit(0)

def hex2rgb(x):
    """
    BED files need CSV RGB, like "128,34,234", so convert hex colors to this
    format.
    """
    return ','.join(
        map(
            str,
            (255 * mpl.colors.colorConverter.to_rgba_array(x))[0,:3].astype(int)
        )
    )

cfg = yaml.load(open(args.config))

# identify the BED file to use
bed = (
    'models/{c[subset]}/{c[states]}/{c[celltype]}_{c[states]}_segments.bed'
    .format(c=cfg['final_model'])
)
if not os.path.exists(bed):
    raise ValueError("Cannot find file '{0}'.".format(bed))


labels = {}
colors = {}
for existing, new, color in cfg['label_mapping']:
    labels[existing] = new
    colors[existing] = hex2rgb(color)

def transform(f):
    f = ff.extend_fields(f, 9)
    label = labels[f[3]]
    color = colors[f[3]]
    f[8] = color
    f[3] = label
    return f

print(pybedtools.BedTool(bed).each(transform))
