## Overview

This repository contains all information and code needed to call enhancers from
histone modification ChIP-seq data in uninduced MEL cells. It boostraps the
installation of requirements, downloads all data needed, and runs multiple
ChromHMM models in parallel. It can be run on a high-performance cluster or
a laptop, all using the same code.

The analysis uses about 25GB of storage and takes a few hrs to run a machine
with 6 cores and 24GB RAM, using a single set of chromatin marks and from 2 to
12 states. The time will vary depending on internet connection speed and
available hardware.

There are three major stages:

1. Download and prepare data (`prepare_data.snakefile`)

    - downloads BAM, BED, and FASTQ files from ENCODE
    - downloads GTF and FASTA from GENCODE
    - downloads enhancer sequences from VISTA Enhancer Browser; lifts over
      any human enhancers to mouse coordinates
    - downloads from any other configured locations
    - merges replicates as needed
    - lifts over to configured target genome as needed
    - prepares annotation data
    - downloads lists of erythroid genes from Li et al and Hughes et al,
      resolves conflicts and typos in the lists, and creates the union of the
      lists

Here is the DAG of jobs for `prepare_data.snakefile`:

![img](_images/prepare_data_rulegraph.png)

The output of this stage is all configured data downloaded and lifted over to
the target assembly, ready for use by ChromHMM.

2. Run ChromHMM (`chromhmm.snakefile`)

    - build models for multiple configured states
    - build models for multiple configured subsets of data
    - run enrichment over all models
    - compute pairwise comparsions of all models


![img](_images/chromhmm_rulegraph.png)

The output of this stage is multiple ChromHMM models using different sets of
input and different numbers of states.

3. Final model preparation (`finalize.py`)

    - add labels to segments
    - color segments by specified colors
    - create final BED files

This is a single script that converts ChromHMM output into a BED file ready for
downstream work.


## Setup

We use [bioconda](https://bioconda.github.io/) to handle software installation
dependencies.

1. Follow the bioconda installation instructions at https://bioconda.github.io to enable bioconda.

2. Create two new conda environments. These create isolated environments that
are completely independent of anything you might already have installed on your
machine. The names are important:

```
conda create -n chromhmm-enhancers --file requirements.txt python=3
conda create -n crossmap crossmap python=2
```

Activate the new `chromhmm-enhancers` environment when following the instructions below:

```
source activate chromhmm-enhancers
```

When you're done, reset back to normal with:

```
source deactivate
```

## Configuration

The main configuration points are the following. If you want to generate
enhancer calls for UMEL cells in the mm9 assembly, you do not need to make any
changes, except perhaps for `CHROMHMM_THREADS` and `MEM_PER_THREAD`, which
depend on your available CPU and RAM.

- `config.yaml`. This file contains lots of documentation to describe exactly
what is configured.

- `data.tsv`. The format of this file is described in `config.yaml`. It defines
URLs or local paths from which to acquire data as well as how they should be
labeled.

- `config/` directory. See the `config.yaml` for a description.

## Run

Activate the environment:

```
source activate chromhmm-enhancers
```

Do a dry-run of the `prepare_data.snakefile`:

```
snakemake --dryrun --snakefile prepare_data.snakefile
```

After lots of output, if you have no config errors, you'll see something like
this, reporting how many jobs will be run:

```
Job counts:
        count   jobs
        1       all
        3       download_chainfiles
        1       download_gencode
        1       download_genome_fasta
        1       download_vista_enhancers
        30      liftovers
        30      local_files
        7       merge
        1       prepare_first_introns
        1       prepare_gene_annotations
        1       prepare_gtfdb
        1       prepare_tsses
        1       prepare_vista_enhancers
        1       strip_gene_version
        1       symlink_erythroid_genes
        81
```

Now you can run the workflow, and provide as many cores as you'd like for
parallel processing where possible. Here for example we're running on 6 cores:

```
snakemake --printshellcmds --snakefile prepare_data.snakefile --cores 6
```

After this completes, do the same with the `chromhmm.snakefile` -- do a dry run, and if no errors then run the workflow:

```
snakemake --dryrun --snakefile chromhmm.snakefile
snakemake --printshellcmds --snakefile chromhmm.snakefile --cores 6
```

## Understanding the output

Below is an annotated directory tree after running both files. Files of
particular note are indicated with a `!`.

The easiest way to get a feel for the different state numbers is to go to
a subset's directory and view all the output images across states. For example,
open these files in an image viewer to compare emission probabilities across
all state numbers. These show which states have high probabilities of having
a particular mark:

```
eog models/subset1/*/emissions_*.png
```


Model comparisons have been automatically performed. These comparison show
which states are shared with other states. Lighter colors indicate a state in
this model (row label) not seen in a model with a different state number
(column label):

```
models/subset1/*/model_comparison_*.png
```

For interpreting states, view the enrichment plots, which show which states are
enriched for the various factors held out of the model:

```
models/subset1/*/*_enrichment.png
```

## Final output

Decide on a model to use. For example, here 6-state model for `subset1` had
biologically interpretable states that were non-redundant.

```
.
├── requirements.txt                              ! # Programs required for this analysis, installed with conda
├── data.tsv                                      ! # CONFIGURE DATA SOURCES HERE
├── config.yaml                                   ! # Top level config file to be edited by user
├── prepare_data.snakefile                        ! # Data prep workflow
├── chromhmm.snakefile                            ! # ChromHMM workflow
├── README.md                                     ! # This file
├── slurm_cluster_config.yaml                     ! # Cluster config for running on SLURM cluster
├── config                                          # Subsets for ChromHMM to run are stored here
│   ├── subset1.tsv                               !   # Original file written by user
│   └── subset1.tsv.filled                            # (Version filled in by the workflow, pointing to the lifted-over, converted-to-BED files.)
├── colorize_lookup                               ! # Colors and labels to use for segmentation classes (used by colorize_model.py, written by user)
├── make_chromhmm_track.sh                        ! # Create a labeled BED12 file of segmentations
├── models                                        ! # Output from ChromHMM
│   └── subset1                                       # (Each configured subset has its own subdirectory)
│       ├── 10                                          # (Each configured number of states has its own sub-subdirectory)
│       │   ├── emissions_10.png                  !       # Emissions, critical for interpreting model
│       │   ├── model_comparison_10.png           !       # Compares other states in other models
│       │   ├── transitions_10.png                !       # Transitions, useful for interpreting model
│       │   ├── umel_10_segments.bed                      # BED file of segmentations, labeled by state number
│       │   ├── umel_enrichment.png               !       # Enrichment with files configured in `data.tsv`. Each column's colorscale is scaled independently. Critical for interpreting model
│       │   ├── umel_uniformscale_enrich ment.png         # Enrichment with files configured in `data.tsv`. Single colorscale for all.
│       │   └── webpage_10.html                           # Aggregates output for the model
│       ├── 11                                         # 11-state model for subset1
│       │   └── ...                                      # files for 11-state model
│       └── ...

# The remaining files documented below are either source files or intermediate files created while generating the above output.


├── binarized                           # binarized files used by ChromHMM
│   └── subset1                           # (each configured subset has its own subdirectory)
│       ├── umel_chr10_binary.txt           # (there is one file for each chromosome in the subdirectory)
│       ...
│       ...
├── colorize_model.py               # Colors the segmentation. Used by make_chromhmm_track.sh.
├── compare_models                  # ChromHMM emissions from all state numbers, used for comparison of states.
│   └── subset1                       # (each configured subset has its own subdirectory)
│       ├── emissions_10.txt
│       ├── emissions_11.txt
│       ...
│       ...
├── data
│   ├── beds_from_bams              # BED files created from lifted-over BAMs, ready for ChromHMM
│   │   ├── umel_control.bed
│   │   ├── umel_dnase.bed
│       ...
│       ...
│   ├── model_data                         # symlinks to lifted-over BAMs, relabeled and ready for use by ChromHMM
│   │   ├── umel_control.bam
│   │   ...
│   │   ...
│   ├── prepare_data                                 # Staging area where data are downloaded and prepared. Can be from any assembly.
│   │   ├── chainfiles                               # chainfiles downloaded for lifting over to target assemblies
│   │   │   └── ...
│   │   ├── enhancers                                # Preparation of files from the VISTA enhancer browser
│   │   │   └── ...
│   │   ├── gene_annotations                         # Preparation of annotations
│   │   │   ├── first_introns.bed
│   │   │   ├── gencode.vM1.annotation.fixed.gtf     # "chr" added to chroms
│   │   │   ├── gencode.vM1.annotation.fixed.gtf.db  # gffutils database
│   │   │   ├── gencode.vM1.annotation.gtf.gz        # downloaded file (source configured in `config.yaml`)
│   │   │   ├── genes.bed                            # just the genes
│   │   │   ├── intergenic.bed                       # just the intergenic space
│   │   │   └── tsses.bed                            # just the TSSes
│   │   └── mm9
│   │       ├── bam                                  # Downloaded BAM files. URLs configured in `data.tsv`
│   │       │   ├── wgEncodeLicrHistoneMelH3k04me1MImmortalC57bl6StdAlnRep1.bam
│   │       │   └── ...
│   │       ├── bed                                  # Downloaded or symlinked BED files. Configured in `data.tsv`
│   │       │   ├── first_introns.bed
│   │       │   └── ...
│   │       ├── bed.gz                               # Downloaded or symlinked gzipped BED files. URLs configured in `data.tsv`
│   │       │   ├── chd1_umel_encode.bed.gz
│   │       │   └── ...
│   ├── prepare_data_lifted_over                     # Mirrors `prepare_data` directory, but data here are have been lifted-over if needed; otherwise symlinked.
│   │   ├── bam
│   │   ├── bed
│   │   ├── bed.gz
│   └── symlinked_for_enrich                  # Directory of files to test for enrichment with ChromHMM models. Configured in `data.tsv`.
│       ├── chd1_umel_encode.bed.gz
│       └── ...
├── erythroid-genes                           # Generation of erythroid gene lists from published data
│   ├── data
│   │   ├── hughes-2014
│   │   │   └── GSE47758_Captured_Regions+500bp.bed.gz
│   │   └── li-2013
│   │       └── TableS3.xlsx
│   ├── erythroid-gene-lists.ipynb            # Jupyter notebook to generate the lists
│   └── gene-lists                            # Final gene lists
│       ├── union-ensembl.gene.txt            # (e.g., union of Li and Hughes genes, in Ensembl accession format)
│       └── ...
├── helpers.py                                                # Helpful functions
├── include                                                   # Data included in the repository -- mostly because it can't be easily downloaded in script
    ├── MEL_Ldb1_induced.bed                                  # LDB1 peaks in IMEL, downloaded from PSU Genome Browser mirror
    └── MEL_Ldb1_uninduced.bed                                # LDB1 peaks in UMEL, downloaded from PSU Geneome Browser mirror

```
