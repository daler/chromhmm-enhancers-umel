"""
Prepares data for Spectacle/ChromHMM model building and inspection
"""

import os
import pandas
import pandas as pd
import pybedtools
from pybedtools import featurefuncs as ff
from gffutils.helpers import asinterval
import gffutils
import yaml
from snakemake.utils import makedirs

shell.prefix("set -o pipefail; set -e; source activate enhancers; ")

config.update(yaml.load(open('config.yaml')))

df = pandas.read_table(config['DATA_TABLE'], comment='#')

to_merge = dict(
    (k, list(v['accession'])) for k, v in df.groupby(['cell', 'mark'])
)


def raw_file_from_row(row):
    """
    Given a row from the data table, return the path to where its raw data
    should be stored.
    """
    return (
        'data/prepare_data/{0[assembly]}/{0[type]}/{0[accession]}.{0[type]}'
        .format(row)
    )


def make_relative_symlink(target, linkname):
    """
    Helper function to create a relative symlink.

    Changes to the dirname of the linkname and figures out the relative path to
    the target before creating the symlink.
    """
    linkdir = os.path.dirname(linkname)
    relative_target = os.path.relpath(target, start=linkdir)
    linkbase = os.path.basename(linkname)
    shell('cd {linkdir}; ln -sf {relative_target} {linkbase}')


# Set up some lookup dictionaries
raw_files = []
lifted_over_files = []
raw_file_lookup = {}
accession_lookup = {}
model_data = []
expression_data = []
for _, row in df.iterrows():
    raw_file = raw_file_from_row(row)
    raw_file_lookup[(row['assembly'], row['type'], row['accession'])] = row
    raw_files.append(raw_file)
    accession_lookup[row['accession']] = row
    lifted_over_files.append(
        'data/prepare_data_lifted_over/{0[type]}/{0[accession]}.{0[type]}'
        .format(row)
    )
    if row['used_for'] == 'model':
        if row['type'] != 'bam':
            raise ValueError('non-bam data specified for model: %s' % row)
        model_data.append(
            'data/model_data/{0[cell]}_{0[mark]}.bam'.format(row))


targets = [
            'data/{0}.fa'.format(config['TARGET_ASSEMBLY']),
            'data/prepare_data/enhancers/enhancer.lbl.gov.bed',
            'data/prepare_data/gene_annotations/genes.bed',
            'data/prepare_data/gene_annotations/intergenic.bed',
            'data/prepare_data/gene_annotations/first_introns.bed',
            'data/prepare_data/gene_annotations/tsses.bed',
            'data/erythroid_genes.txt',
        ] + lifted_over_files + model_data + expression_data


# Create all output directories needed
dirs = list(set([os.path.dirname(i) for i in targets]))
logdirs = [os.path.join('logs', i) for i in dirs]
makedirs(dirs + logdirs)

rule all:
    input: targets


# ----------------------------------------------------------------------------
# Most of the BAM files to be used for the model have replicates. So we merge
# together replicates for model building.
rule merge:
    input:
        lambda wc: expand('data/prepare_data_lifted_over/bam/{accession}.bam', accession=to_merge[(wc.cell, wc.mark)])
    output: 'data/model_data/{cell}_{mark}.bam'
    run:
        if len(input) == 1:
            make_relative_symlink(input[0], output[0])
        else:
            inputs = ' '.join(['INPUT=%s' % i for i in input])
            #shell('conda run -n enhancers -- samtools merge {output} {input}')
            shell('picard MergeSamFiles {inputs} OUTPUT={output} ')


# ----------------------------------------------------------------------------
# Download or make symlinks as needed
#
# If it's a local file, then it must exist, and therefore we use it as an input
# depdendency. For example, for the introns and genes files this will trigger
# those files to be created in their respective rules.
#
def local_or_remote(wildcards):
    """
    If url starts with "http", return empty list, otherwise assume it's
    a relative path and return it.
    """
    row = raw_file_lookup[
        (
            wildcards.assembly,
            wildcards.type,
            wildcards.accession
        )
    ]
    if row['url'].startswith('http'):
        return []
    return row['url']

rule local_files:
    input: local_or_remote
    output: 'data/prepare_data/{assembly}/{type}/{accession}.{type}'
    run:
        row = raw_file_lookup[
            (
                wildcards.assembly,
                wildcards.type,
                wildcards.accession)]
        url = row['url']
        if url.startswith('http'):
            shell('wget -nv -O {output} {url}')

        else:
            full_output_path = os.path.abspath(output[0])
            full_raw_path = os.path.abspath(url)
            shell('ln -sf {full_raw_path} {full_output_path}')


# ----------------------------------------------------------------------------
# Perform liftovers if a file isn't in the target assembly; otherwise simply
# create symlink.
#
# Final files go in data/prepare_data_lifted_over/{type}/{accession}.{type}
def _liftovers_inputs(wildcards):
    """
    Input-handling function for `liftovers` rule.

    If the file is already in the target assembly, simply return the raw file.
    Otherwise return the raw file and the required chainfile.

    The `liftovers` rule will then figure out whether to simply symlink or to
    perform the liftover.
    """
    row = accession_lookup[wildcards.accession]
    raw_file = raw_file_from_row(row)
    if row['assembly'] != config['TARGET_ASSEMBLY']:
        return (
            raw_file,
            'data/prepare_data/chainfiles/'
            '{0[assembly]}To{1}.over.chain.gz'
            .format(row, config['TARGET_ASSEMBLY'].title())
        )
    return raw_file


rule liftovers:
    input: _liftovers_inputs
    output: 'data/prepare_data_lifted_over/{type}/{accession}.{type}'
    run:

        # CrossMap has different ways of sending to stdout depending on
        # subprogram (here, bam or bed). For bam, it tries to add an extra .bam
        # extension to the output filename so we use an output file of "-". For bed, trying to
        # redirect actually dumps the full before-and-after on each line, so we
        # need to call it differently.
        #
        # Also, it does not gzip the output, so we have to handle that as well.
        if len(input) == 2:
            if wildcards.type in ('bed', 'bed.gz'):
                if wildcards.type == 'bed.gz':
                    outfile = output[0] + '.tmp'
                else:
                    outfile = output[0]
                shell(
                    'source activate crossmap; CrossMap.py '
                    'bed {input[1]} {input[0]} {outfile}'
                )
                if wildcards.type == 'bed.gz':
                    shell('gzip -c {outfile} > {output}')

            elif _type == 'bam':
                shell(
                    'source activate crossmap; CrossMap.py '
                    '{_type} {input[1]} {input[0]} - > {output}'
                )

        # if no chainfile, then just make the symlink
        else:
            make_relative_symlink(input[0], output[0])



# ----------------------------------------------------------------------------
# Download any necessary chainfiles for doing liftovers.
rule download_chainfiles:
    output: 'data/prepare_data/chainfiles/{from_assembly}To{to_assembly}.over.chain.gz'
    shell:
        'wget -nv -O {output} http://hgdownload.cse.ucsc.edu/goldenPath/'
        '{wildcards.from_assembly}/liftOver/'
        '{wildcards.from_assembly}To{wildcards.to_assembly}.over.chain.gz'


# ----------------------------------------------------------------------------
# Download FASTA file from VISTA enhancer browser
rule download_vista_enhancers:
    output: 'data/prepare_data/enhancers/enhancer.lbl.gov.fa'
    shell:
        'wget -nv -O {output} '
        '"http://enhancer.lbl.gov/cgi-bin/imagedb3.pl?page_size=100;show=1;'
        'search.result=yes;form=search;search.form=no;action=search;search.sequence=1"'


# ----------------------------------------------------------------------------
# Download GENCODE gtf file
rule download_gencode:
    output: 'data/prepare_data/gene_annotations/' + os.path.basename(config['GTF'])
    shell:
        'wget -nv -O {output} {config[GTF]}'


# ----------------------------------------------------------------------------
# Strip the gene version off of gene_ids
rule strip_gene_version:
    input: rules.download_gencode.output[0]
    output: 'data/prepare_data/gene_annotations/{0}'.format(os.path.basename(config['GTF']).replace('.gtf.gz', '.fixed.gtf'))
    run:
        with open(str(output), 'w') as fout:
            for i, f in enumerate(gffutils.iterators.DataIterator(str(input))):
                f['gene_id'][0] = f['gene_id'][0].split('.')[0]
                fout.write(str(f) + '\n')


# ----------------------------------------------------------------------------
# Make a db from the GTF file.
rule prepare_gtfdb:
    input: rules.strip_gene_version.output[0]
    output: rules.strip_gene_version.output[0] + '.db'
    run:
        gffutils.create_db(
            str(input[0]),
            output[0],
            disable_infer_genes=True,
            disable_infer_transcripts=True,
            force=True,
            id_spec={
                'gene': 'gene_id',
                'transcript': 'transcript_id'},
            verbose=True,
        )


# ----------------------------------------------------------------------------
# Create gene, intron, intergenic regions from the GENCODE annotation
rule prepare_gene_annotations:
    input: rules.prepare_gtfdb.output[0]
    output:
        genes     ='data/prepare_data/gene_annotations/genes.bed',
        intergenic='data/prepare_data/gene_annotations/intergenic.bed'
    run:
        db = gffutils.FeatureDB(str(input))
        def gen():
            for gene in db.features_of_type('gene'):
                yield ff.gff2bed(asinterval(gene))

        tmp = pybedtools.BedTool(gen()).saveas()

        pybedtools.BedTool(tmp).sort(output=output.genes, faidx=config['TARGET_ASSEMBLY'])\
            .complement(g=config['TARGET_ASSEMBLY'], output=str(output.intergenic))


# ----------------------------------------------------------------------------
# Create TSS file.
rule prepare_tsses:
    input: rules.prepare_gtfdb.output[0]
    output: 'data/prepare_data/gene_annotations/tsses.bed'
    run:
        db = gffutils.FeatureDB(str(input))
        def gen():
            for gene in db.features_of_type('gene'):
                for child in db.children(gene, level=1):
                    yield ff.TSS(ff.gff2bed(asinterval(child)), upstream=100, downstream=100)
        pybedtools.BedTool(gen()).sort().merge(output=output[0])


# ----------------------------------------------------------------------------
# Identify the 5'-most intron 
rule prepare_first_introns:
    input: rules.prepare_gtfdb.output[0]
    output: 'data/prepare_data/gene_annotations/first_introns.bed'
    run:
        db = gffutils.FeatureDB(str(input))
        def first_intron_generator():
            """
            Uses gffutils and db to identify the 5'-most intron of a gene.
            """
            for gene in db.features_of_type('gene'):

                # get the first intron of each transcript
                transcripts = db.children(gene, level=1)
                introns = []
                for transcript in transcripts:

                    # get the first 2 exons of the transcript
                    exons = list(db.children(transcript, featuretype='exon', level=1, order_by='start'))
                    if transcript.strand == '-':
                        first_exons = exons[-2:]
                    else:
                        first_exons = exons[:2]

                    inter = list(
                        db.interfeatures(
                            first_exons,
                            new_featuretype='intron',
                            merge_attributes=False,
                            update_attributes={'gene_id': gene['gene_id']}
                        )
                    )
                    # no introns for this transcript....
                    if len(inter) == 0:
                        continue

                    # ...otherwise there had better be just one.
                    assert len(inter) == 1
                    introns.append(inter[0])

                # No introns at all for this gene.
                if len(introns) == 0:
                    continue

                five_prime_most = True
                if five_prime_most:
                    # Get the 5'-most intron
                    introns = sorted(introns, key=lambda x: x.start)
                    if gene.strand == '-':
                        first_intron = introns[-1]
                    else:
                        first_intron = introns[0]
                    yield asinterval(first_intron)

                else:
                    for intron in introns:
                        yield ff.gff2bed(asinterval(first_intron), 'gene_id')
        pybedtools.BedTool(first_intron_generator()).each(ff.gff2bed).saveas(output[0])


# ----------------------------------------------------------------------------
# Extracts the coordinates from the downloaded FASTA file from the VISTA
# enhancer browser; lifts over and concatenates mouse (mm9) and human (hg19)
# coords to one final mm10 file.
rule prepare_vista_enhancers:
    input:
        enh=rules.download_vista_enhancers.output[0],
        hg19ToMm10_chainfile='data/prepare_data/chainfiles/hg19ToMm10.over.chain.gz',
        mm9ToMm10_chainfile='data/prepare_data/chainfiles/mm9ToMm10.over.chain.gz'
    output: 'data/prepare_data/enhancers/enhancer.lbl.gov.bed'
    run:
        def gen():
            """
            FASTA file contains coordinates and metadata headers, e.g.

                >Human|chr16:80372593-80373755 | element 4 | positive  | neural
                tube[6/10] | hindbrain (rhombencephalon)[10/10] | midbrain
                (mesencephalon)[10/10]

            This extracts the information and converts to BED format, putting
            the extra non-coordinate info into the name field.
            """
            for line in open(input.enh):
                line = line.replace('<pre>', '')
                if not line.startswith('>'):
                    continue
                toks = line.replace('>', '', 1).strip().split('|')
                toks = [i.strip() for i in toks]
                genome, coords, element, result = toks[:4]
                if result != 'positive':
                    continue
                if len(toks) > 4:
                    name = "|".join(toks[4:] + [genome])
                else:
                    name = '.'
                chrom, startstop = coords.split(':')
                start, stop = startstop.split('-')
                name = name.replace('\t', ' ').replace(' ', '_')
                yield pybedtools.create_interval_from_list([chrom, start, stop, '"' + name + '"'])

        x = pybedtools.BedTool(gen()).saveas()

        # Save incremental output for troubleshooting (e.g., to confirm an
        # enhancer missing in the final output was unmapped)
        human_fn = os.path.join(
            os.path.dirname(output[0]),
            'enhancer.lbl.gov.hg19.bed'
        )
        mouse_fn = os.path.join(
            os.path.dirname(output[0]),
            'enhancer.lbl.gov.mm9.bed'
        )
        hg19_to_mm10 = os.path.join(
            os.path.dirname(output[0]),
            'enhancer.lbl.gov.hg19-to-mm10.bed'
        )
        mm9_to_mm10 = os.path.join(
            os.path.dirname(output[0]),
            'enhancer.lbl.gov.mm9-to-mm10.bed'
        )

        # Sort enhancers into mouse or human
        human = open(human_fn, 'w')
        mouse = open(mouse_fn, 'w')
        for i in x:
            if 'Human' in i[-1]:
                human.write(str(i))
            elif 'Mouse' in i[-1]:
                mouse.write(str(i))
        human.close()
        mouse.close()

        # Liftover human and mouse to mm10
        shell('liftOver -bedPlus=4 {human.name} '
              '{input.hg19ToMm10_chainfile} '
              '{hg19_to_mm10} {hg19_to_mm10}.unmapped')

        shell('liftOver -bedPlus=4 {mouse.name} '
              '{input.mm9ToMm10_chainfile} '
              '{mm9_to_mm10} '
              '{mm9_to_mm10}.unmapped')

        # combine mm10 coords into final file.
        shell('cat {mm9_to_mm10} {hg19_to_mm10} > {output}')



# ----------------------------------------------------------------------------
# Download genome fasta file
rule download_genome_fasta:
    output: 'data/{0}.fa'.format(config['TARGET_ASSEMBLY'])
    shell:
        'wget -O {output}.gz {config[GENOME_FASTA]}; '
        'gunzip {output}'


# ----------------------------------------------------------------------------
# Download the transcriptome fasta file
rule transcriptome_sequence:
    output: 'data/{0}.transcriptome.fa.gz'.format(config['TARGET_ASSEMBLY'])
    shell:
        'wget -O {output} {config[TRANSCRIPTOME]} && touch {output}'


# ----------------------------------------------------------------------------
# Create a kallisto index
rule kallisto_index:
    input: rules.transcriptome_sequence.output[0]
    output: 'data/{0}.transcriptome.kallisto.idx'.format(config['TARGET_ASSEMBLY'])
    shell:
        'kallisto index -i {output} {input}'


# ----------------------------------------------------------------------------
# Quantify expression with kallisto
rule kallisto_quant:
    input:
        fastqs=lambda wc: expand('data/prepare_data_lifted_over/fastq/{accession}.fastq', accession=to_merge[(wc.cell, wc.mark)]),
        idx=rules.kallisto_index.output[0]
    output: 'data/expression/{cell}_{mark}/abundance.tsv'
    shell:
        'kallisto quant -i {input.idx} '
        '-o $(dirname {output}) '
        '--single '
        '--fragment-length=200 '
        '--sd 20 '
        '--plaintext '
        '{input.fastqs}'


# ----------------------------------------------------------------------------
# Aggregate kallisto output by summing transcript abundance by gene
rule kallisto_summarize:
    input: rules.kallisto_quant.output
    output: 'data/expression/{cell}_{mark}/abundance.gene.tsv'
    run:
        df = pandas.read_table(input[0])
        df['gene'] = df.target_id.apply(lambda x: x.split('|')[1].split('.')[0])
        df.groupby('gene').agg(sum)[['tpm']].to_csv(output[0], sep='\t')


# ----------------------------------------------------------------------------
# Run Jupyter notebook to generate the lists of erythroid genes
rule erythroid_genes:
    input: 'erythroid-genes/erythroid-gene-lists.ipynb'
    output: 'erythroid-genes/gene-lists/union-ensembl.gene.txt'
    shell:
        'jupyter nbconvert --execute {input}'


# ----------------------------------------------------------------------------
# Symlink the union erythroid genes over to the data dir
rule symlink_erythroid_genes:
    input: rules.erythroid_genes.output
    output: 'data/erythroid_genes.txt'
    run:
        make_relative_symlink(input[0], output[0])


# ----------------------------------------------------------------------------
# Download interaction data from Schoenfelder et al
rule download_interactions:
    output: 'data/mm9/E-MTAB-2414.additional.1.zip'
    shell:
        'wget -O {output} {config[ZIP_URL]}'


# ----------------------------------------------------------------------------
# Unzip interaction data from Schoenfelder et al
rule unzip_interactions:
    input: 'data/mm9/E-MTAB-2414.additional.1.zip'
    output:
        'data/mm9/ESC_promoter_other_significant_interactions.txt',
        'data/mm9/ESC_promoter_promoter_significant_interactions.txt',
        'data/mm9/FLC_promoter_other_significant_interactions.txt',
        'data/mm9/FLC_promoter_promoter_significant_interactions.txt'
    run:
        shell('cd $(dirname {input}); unzip $(basename {input})')
        shell('touch {output}')


# ----------------------------------------------------------------------------
# Convert data for promoter-other interactions to BEDPE format and paired-end
# BAM format.
rule prepare_promoter_other:
    input: 'data/mm9/{celltype}_promoter_other_significant_interactions.txt'
    output:
        bedpe=temp('data/mm9/{celltype}_promoter_{distal}_significant_interactions.bedpe'),
        bam='data/mm9/{celltype}_promoter_{distal}_significant_interactions.bam',
        bai='data/mm9/{celltype}_promoter_{distal}_significant_interactions.bam.bai'
    run:
        if 'promoter_other' in str(input[0]):
            bedpe_coord_cols = ['chr bait', 'start bait', 'end bait', 'chr', 'start', 'end']
        elif 'promoter_promoter' in str(input[0]):
            bedpe_coord_cols = ['chr', 'start', 'end', 'chr.1', 'start.1', 'end.1']

        df = pandas.read_table(str(input[0]))
        df['strand'] = '.'
        df['name'] = df.index
        cols = bedpe_coord_cols + ['name', 'log(observed/expected)', 'strand', 'strand']
        df[cols].to_csv(str(output.bedpe), sep='\t', index=False, header=False)
        bt = pybedtools.BedTool(str(output.bedpe))
        bt = bt.bedpe_to_bam(g=config['TARGET_ASSEMBLY'])
        shell('samtools sort {bt.fn} -o {output.bam}')
        shell('samtools index {output.bam}')



# vim: ft=python
