"""
This workflow expects `prepare_data.snakefile` to have been completed
successfully such that all data specified in DATA_TABLE have been downloaded,
merged, and lifted over as appropriate.
"""

import pandas as pd
from snakemake.utils import makedirs
import os
import yaml
shell.prefix("set -o pipefail; set -e; ")

config.update(yaml.load(open('config.yaml')))

# See config.yaml for details on these options.
ASSEMBLY = config['TARGET_ASSEMBLY']
DATA_TABLE = config['DATA_TABLE']
MEM_PER_THREAD = config['MEM_PER_THREAD']
BINSIZE = config['BINSIZE']
CHROMHMM_THREADS = config['CHROMHMM_THREADS']

localrules: dummy_symlink, copy_models_to_compare, aggregate_plots

# ----------------------------------------------------------------------------
# Data table prep
#
# Read the data table
df = pd.read_table(DATA_TABLE, comment='#')

# Extract just the data used for the models
sampletable = df[df['used_for'] == 'model']

# Identify the unique (cell, mark) values and generate underscore-delimited
# keys.
cells_marks = [tuple(i) for i in sampletable[['cell', 'mark']].drop_duplicates().values]
cells = list(sampletable['cell'].unique())
marks = list(sampletable['mark'].unique())
keys = ['_'.join(i) for i in cells_marks]

# We use the chromosomes to construct expected filenames. We don't really need
# to generate ALL chrom names, just enough to detect if a run has been
# completed. If you're using a different assembly, adjust as appropriate.
chroms = [line.split()[0] for line in open(ASSEMBLY)][:3]

# We expect these files to be merged
merged = expand('data/model_data/{key}.bam', key=keys)

# BAMs are converted to bed
beds = expand('data/beds_from_bams/{key}.bed', key=keys)

# Construct the list of files we will use for enrichment testing.
enriched_df = df[df['used_for'] == 'enrich']
enriched_beds = list(set(expand(
    'data/prepare_data_lifted_over/{_type}/{accession}.{_type}', zip,
    _type=enriched_df['type'], accession=enriched_df['accession']
)))

# ChromHMM enrichment program needs all files in a single directory. This
# copies everything over, or gzips as needed.
linked = []
for i in enriched_beds:
    if i.endswith('.bed'):
        i = i + '.gz'
    linked.append(os.path.join('data/symlinked_for_enrich', os.path.basename(i)))


def subset_expand(subset, pattern):
    celltypes = pd.read_table(os.path.join('config', subset + '.tsv')).cell
    states = config['SUBSETS'][subset]
    return sorted(set(expand(pattern, subset=subset, cell=celltypes, chrom=chroms, state=states)))

_patterns = {
    'binarized': 'binarized/{subset}/{cell}_{chrom}_binary.txt',
    'compare_models': [
        'models/{subset}/{state}/{cell}_model_comparison_{state}.txt',
        'compare_models/{subset}/emissions_{state}.txt'],
    'learn_model': 'models/{subset}/{state}/{cell}_{state}_segments.bed',
    'enrichment': 'models/{subset}/{state}/{cell}_enrichment.txt',
    'agg': [
        'models/{subset}/{state}/{cell}_aggregated_plots.pdf',
        'models/{subset}/{state}/{cell}_aggregated_plots.png'
    ],
}

targets_dict = {}
for k, v in _patterns.items():
    t = []
    if isinstance(v, str):
        v = [v]
    for subset in config['SUBSETS'].keys():
        for pattern in v:
            t.extend(subset_expand(subset, pattern))
    targets_dict[k] = t


def make_relative_symlink(target, linkname):
    linkdir = os.path.dirname(linkname)
    relative_target = os.path.relpath(target, start=linkdir)
    linkbase = os.path.basename(linkname)
    shell('cd {linkdir}; ln -sf {relative_target} {linkbase}')


def fill_filetable(subset_filename):
    """
    Given the filename of a defined subset, fill in the BED files for each
    cell/mark.
    """
    df = pd.read_table(subset_filename).fillna("")

    def add_control(x):
        if x.endswith('_'):
            return ""
        return x + '.bed'

    df['mark_filename'] = (df['cell'] + '_' + df['mark']).apply(lambda x: x + '.bed')
    df['control_filename'] = (df['cell'] + '_' + df['control']).apply(add_control)
    return df[['cell', 'mark', 'mark_filename', 'control_filename']]


rule all:
    input: targets_dict.values()


# ----------------------------------------------------------------------------
# Convert BAM to bed
rule bam_to_bed:
    input: 'data/model_data/{cell}_{mark}.bam'
    output: 'data/beds_from_bams/{cell}_{mark}.bed'
    shell:
        '''
        bedtools bamtobed -i {input} > {output}.tmp
        mv {output}.tmp {output}
        '''


# ----------------------------------------------------------------------------
# Binarize the bed files for each subset. This is run once for each defined
# subset.

def beds_for_subset(wc):
    beds = []
    df = pd.read_table(os.path.join('config', subset + '.tsv'))
    for _, row in df.iterrows():
        beds.append('data/beds_from_bams/{0.cell}_{0.mark}.bed'.format(row))
    return beds


rule binarize_bed:
    input:
        beds=beds_for_subset,
        cellmarkfiletable='config/{subset}.tsv'
    output:
        binarized=expand('binarized/{{subset}}/{{cell}}_{chrom}_binary.txt', chrom=chroms)
    log: 'binarized/{subset}/BinarizeBed.log'
    run:
        fill_filetable(
            input.cellmarkfiletable).to_csv(
                input.cellmarkfiletable + '.filled',
                header=False, index=False, sep='\t')
        shell(
            'ChromHMM.sh -Xmx{MEM_PER_THREAD}M '
            'BinarizeBed '
            '-b {BINSIZE} '
            '{ASSEMBLY} '
            'data/beds_from_bams '
            '{input.cellmarkfiletable}.filled '
            'binarized/{wildcards.subset} > {log} 2>&1'
        )


# ----------------------------------------------------------------------------
# Learn model. This is run once for each defined subset for each defined state
# number.

def _learn_model_inputs(wc):
    return subset_expand(wc.subset, 'binarized/{subset}/{cell}_{chrom}_binary.txt')

rule learn_model:
    input: _learn_model_inputs
    output:
        png='models/{subset}/{state}/emissions_{state}.png',
        txt='models/{subset}/{state}/emissions_{state}.txt',
        transitions_txt='models/{subset}/{state}/transitions_{state}.txt',

        # Since this output file also contains "{cell}", its wildcards do not
        # match the other output files. See the `dummy_symlink` rule below.
        #segments='models/{subset}/{state}/{cell}_{state}_segments.bed'

    log: 'models/{subset}/{state}/LearnModel.log'
    threads: config['CHROMHMM_THREADS']
    shell:
        "ChromHMM.sh -Xmx{MEM_PER_THREAD}M "
        "LearnModel "
        "-b {BINSIZE} "
        "-p {threads} "
        "-noenrich "
        "-nobrowser "
        "-printposterior "
        "binarized/{wildcards.subset} "
        "$(dirname {output.txt}) "
        "{wildcards.state} "
        "{ASSEMBLY} "
        "> {log} 2>&1"


rule dummy_symlink:
    input: 'models/{subset}/{state}/emissions_{state}.txt'
    output: 'models/{subset}/{state}/{cell}_{state}_segments.bed'
    shell:
        'touch {output}'

# ----------------------------------------------------------------------------
# Copy enrichment files to directory gzipping if needed.
# Note: symlinking would be a good option if your filesystem supports the touch
# -h option.
rule copy_to_enrichment_dir:
    input: enriched_beds
    output: linked
    run:
        for target, linkname in zip(input, output):
            if target.endswith('.gz'):
                shell("cp {target} {linkname}")
            else:
                shell("gzip -c {target} > {linkname}")


# ----------------------------------------------------------------------------
# Create heatmaps of enrichment of each state. Runs once for each subset for
# each state number.
#
rule overlap_enrichment:
    input:
        segments='models/{subset}/{state}/{cell}_{state}_segments.bed',
        beds=linked,
    output:
        txt='models/{subset}/{state}/{cell}_enrichment.txt',
        png='models/{subset}/{state}/{cell}_enrichment.png',
        svg='models/{subset}/{state}/{cell}_enrichment.svg',
    log: 'models/{subset}/{state}/{cell}_OverlapEnrichment.log'
    params: prefix='models/{subset}/{state}/{cell}_enrichment'
    shell:
        "ChromHMM.sh -Xmx{MEM_PER_THREAD}M "
        "OverlapEnrichment {input.segments} "
        "data/symlinked_for_enrich "
        "{params.prefix} "
        "> {log} 2>&1"


# ----------------------------------------------------------------------------
# Copy over the models for comparison. Runs once after all models have been
# created.
rule copy_models_to_compare:
    input: 'models/{subset}/{state}/emissions_{state}.txt'
    output: 'compare_models/{subset}/emissions_{state}.txt'
    shell:
        'cp {input} {output}'


# ----------------------------------------------------------------------------
# Create heatmaps of model comparisons
rule compare_models:
    input:
        reference='models/{subset}/{state}/emissions_{state}.txt',
        to_compare=rules.copy_models_to_compare.output
    output: 'models/{subset}/{state}/{cell}_model_comparison_{state}.txt'
    params: prefix='models/{subset}/{state}/{cell}_model_comparison_{state}'
    shell:
        "ChromHMM.sh -Xmx{MEM_PER_THREAD}M "
        "CompareModels {input.reference} $(dirname {input.to_compare}) {params.prefix} "

# ----------------------------------------------------------------------------
# Aggregate plots
rule aggregate_plots:
    input:
        emissions='models/{subset}/{state}/emissions_{state}.txt',
        transitions='models/{subset}/{state}/transitions_{state}.txt',
        enrichment='models/{subset}/{state}/{cell}_enrichment.txt',
        uniform_enrichment='models/{subset}/{state}/{cell}_enrichment.txt',
        comparison='models/{subset}/{state}/{cell}_model_comparison_{state}.txt',
        script='viz.py',
    output:
        png='models/{subset}/{state}/{cell}_aggregated_plots.png',
        pdf='models/{subset}/{state}/{cell}_aggregated_plots.pdf'
    script:
        'viz.py'



# vim: ft=python
