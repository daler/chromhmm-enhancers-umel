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
states = config['STATES']
subsets = config['SUBSETS']


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
enriched_beds = expand(
    'data/prepare_data_lifted_over/{_type}/{accession}.{_type}', zip,
    _type=enriched_df['type'], accession=enriched_df['accession']
)

# ChromHMM enrichment program needs all files in a single directory. This
# copies everything over, or gzips as needed.
linked = []
for i in enriched_beds:
    if i.endswith('.bed'):
        i = i + '.gz'
    linked.append(os.path.join('data/symlinked_for_enrich', os.path.basename(i)))
linked = list(set(linked))

# ChromHMM model comparison needs all models' emissions file in the same
# directory, so this copies everything over.
linked_emissions = expand('compare_models/{subset}/emissions_{state}.txt', state=states, subset=subsets)

# Target files for model comparison
compared = expand('models/{subset}/{state}/model_comparison_{state}.txt', state=states, subset=subsets)

# binarized output files. Example output:
# e14.5-liver  chr1
# h3k27ac  h3k27me3 h3k36me3 h3k4me1 h3k4me2 h3k4me3 h3k9ac h3k9me3
# 0        0        0        0       0       0       0      0
# 0        0        0        0       0       0       0      0
# 0        0        0        0       0       0       0      0
binarized = expand('binarized/{subset}/{cell}_{chrom}_binary.txt', cell=cells, chrom=chroms, subset=subsets)

# Expected files for model output (here we just detect the emissions PNG)
models = expand('models/{subset}/{state}/emissions_{state}.png', subset=subsets, state=states)

# Expected files for enrichment output
enrich = expand('models/{subset}/{state}/{cell}_enrichment.txt', subset=subsets, state=states, cell=cells)
enrich += expand('models/{subset}/{state}/{cell}_uniformscale_enrichment.txt', subset=subsets, state=states, cell=cells)


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
    input:  models + linked + enrich + linked_emissions + compared


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
rule binarize_bed:
    input:
        beds=beds,
        cellmarkfiletable='config/{subset}.tsv'
    output:
        binarized=expand('binarized/{{subset}}/{cell}_{chrom}_binary.txt', cell=cells, chrom=chroms)
    log: 'binarized/{subset}/BinarizeBed.log'
    run:
        fill_filetable(
            input.cellmarkfiletable).to_csv(
                input.cellmarkfiletable + '.filled',
                header=False, index=False, sep='\t')
        shell(
            '''
            ChromHMM.sh -Xmx{MEM_PER_THREAD}M \\
            BinarizeBed \\
            -b {BINSIZE} \\
            {ASSEMBLY} \\
            data/beds_from_bams \\
            {input.cellmarkfiletable}.filled \\
            binarized/{wildcards.subset} > {log} 2>&1''')

# ----------------------------------------------------------------------------
# Learn model. This is run once for each defined subset for each defined state
# number.
rule learn_model:
    input: rules.binarize_bed.output
    output:
        'models/{subset}/{state}/emissions_{state}.png',
        'models/{subset}/{state}/emissions_{state}.txt',
        expand('models/{{subset}}/{{state}}/{cell}_{{state}}_segments.bed', cell=cells)
    log: 'models/{subset}/{state}/LearnModel.log'
    threads: config['CHROMHMM_THREADS']
    run:
        outdir = os.path.dirname(output[0])
        shell(
            """
            ChromHMM.sh -Xmx{MEM_PER_THREAD}M \\
            LearnModel \\
            -b {BINSIZE} \\
            -p {threads} \\
            -noenrich \\
            -nobrowser \\
            binarized/{wildcards.subset} \\
            {outdir} \\
            {wildcards.state} \\
            {ASSEMBLY} \\
            > {log} 2>&1
            """)

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
# This creates both column-scaled and uniform color scale versions.
rule overlap_enrichment:
    input:
        segments='models/{subset}/{state}/{cell}_{state}_segments.bed',
        beds=linked,
    output:
        txt='models/{subset}/{state}/{cell}_enrichment.txt',
        png='models/{subset}/{state}/{cell}_enrichment.png',
        svg='models/{subset}/{state}/{cell}_enrichment.svg',
        txtu='models/{subset}/{state}/{cell}_uniformscale_enrichment.txt',
        pngu='models/{subset}/{state}/{cell}_uniformscale_enrichment.png',
        svgu='models/{subset}/{state}/{cell}_uniformscale_enrichment.svg',
    log: 'models/{subset}/{state}/{cell}_OverlapEnrichment.log'
    run:
        inputcoorddir = 'data/symlinked_for_enrich'

        outfileprefix = output.txt.replace('.txt', '')
        shell(
            """
            ChromHMM.sh -Xmx{MEM_PER_THREAD}M \\
            OverlapEnrichment {input.segments} {inputcoorddir} {outfileprefix} \\
            > {log} 2>&1
            """)

        outfileprefix = output.txtu.replace('.txt', '')
        shell(
            """
            ChromHMM.sh -Xmx{MEM_PER_THREAD}M \\
            OverlapEnrichment -uniformscale {input.segments} {inputcoorddir} {outfileprefix} \\
            > {log} 2>&1
            """)


# ----------------------------------------------------------------------------
# Copy over the models for comparison. Runs once after all models have been
# created.
rule copy_models_to_compare:
    input: expand('models/{{subset}}/{state}/emissions_{state}.txt', state=states)
    output: expand('compare_models/{{subset}}/emissions_{state}.txt', state=states)
    run:
        for target, linkname in zip(input, output):
            shell('cp {target} {linkname}')


# ----------------------------------------------------------------------------
# Create heatmaps of model comparisons
rule compare_models:
    input:
        reference='models/{subset}/{state}/emissions_{state}.txt',
        to_compare=rules.copy_models_to_compare.output
    output: 'models/{subset}/{state}/model_comparison_{state}.txt'
    run:
        compdir = os.path.dirname(input.to_compare[0])
        outprefix = output[0].replace('.txt', '')

        shell(
            """
            ChromHMM.sh -Xmx{MEM_PER_THREAD}M \\
            CompareModels {input.reference} {compdir} {outprefix}
            """
        )
# vim: ft=python
