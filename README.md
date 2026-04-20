# remapping_workflow

This is a nextflow workflow for remapping reads to sets of virus sequences.  In our lab we use it to validate virus sequences identified by metagenomic sequencing.  It also outputs mapping statistics and plots that can be used for downstream analyses.  It could be used as a stand-alone pipeline or as a subworkflow.  We designed this to be used for virus sequences but it could also be used for other purposes. 

- [Running the pipeline](#Running-the-pipeline)
    - [Running from github](#Running-from-github)
    - [Running test datasets](#Running-test-datasets)
    - [Running a specific version:](#Running-a-specific-version)
    - [Making sure that you are running the latest version](#Making-sure-that-you-are-running-the-latest-version)
    - [Running by cloning the pipeline's repository](#Running-by-cloning-the-pipeline-repository)
- [Inputs](#Inputs)
    - [Input fastq](#Input-fastq)
    - [Mapping samplesheet](#Mapping-samplesheet)
    - [Reference sequence fasta files](#Reference-sequence-files).
- [Outputs](#Outputs)
    - [Output directory](#Output-directory)
- [Configurable parameters](#Configurable-parameters)
- [Software Dependencies](#Software-dependencies)
    - [Nextflow](#Nextflow)
    - [Singularity ](#Singularity)
    - [Conda](#Conda)
    - [R libraries](#R-libraries)


## Running the pipeline

See the [dependencies section](#dependencies) below for information about the main dependencies required for running this pipeline(including nextflow and singularity).

### Running from github

A simple way to run the pipeline is directly from github, like this:

```
nextflow run stenglein-lab/remapping_workflow -resume   --mapping_samplesheet mapping_samplesheet.txt --fastq_dir /path/to/fastq/directory -profile singularity
```

You must specify two required inputs to the pipeline: the path to a mapping samplesheet and the path to a directory containing input fastq.  See [this section](#inputs) for more information on required inputs.

### Running test datasets

The pipeline includes several small datasets (<= 1000 read pairs) that are derived from real known positive (or known negative) datasets.  These are included in the [test directory](./test/) of the repository.  These datasets serve as positive and negative controls and allow you to test that the pipeline is working as expected.  To use these test datasets, run with the test profile, for instance:

```
nextflow run stenglein-lab/remapping_workflow -profile singularity,test
```

Or to run with conda:
```
nextflow run stenglein-lab/remapping_workflow -profile conda,test
```

The results of the test run will be placed in a `test/results` sub-directory.

### Running a specific version

To run a specific version of the pipeline, use the -r option, for example:

```
nextflow run stenglein-lab/remapping_workflow -profile singularity,test -r v1.0
```

### Making sure that you are running the latest version

Nextflow will cache the pipeline in `$HOME/.nextflow/assets/` and continue to use the cached version, even if the pipeline has newer versions on github.  To remove the locally cached version, which will force nextflow to download and cache the latest version, run:

```
nextflow drop stenglein-lab/remapping_workflow
nextflow run stenglein-lab/remapping_workflow -profile singularity,test
```

Alternatively, you can just delete the cached pipeline directory:
```
rm -rf ~/.nextflow/assets/stenglein-lab/remapping_workflow/
```

Running nextflow pipelines from github is [described in more detail here](https://www.nextflow.io/docs/latest/sharing.html).

### Running by cloning the pipeline repository

It is also possible to download the pipeline code to a directory of your choosing.  This can be useful if, for instance, you want to modify or debug the code.  You can do this by cloning the repository (or a fork of the repository):

```
git clone https://github.com/stenglein-lab/remapping_workflow.git
cd remapping_workflow
nextflow run main.nf -resume --mapping_samplesheet mapping_samplesheet.txt --fastq_dir /path/to/fastq/directory -profile singularity
```

## Inputs

The pipeline requires three main inputs:

1. [Sequence datasets in fastq format](#Input-fastq).

2. [A mapping samplesheet](#Mapping-samplesheet).

3. [Reference sequence fasta files](#Reference-sequence-files).


### Input fastq

Input sequence data should single-end or paired-end Illumina data.  Paired-end reads should be in separate read1 and read2 files. Files must be .gz compressed (names should end in `.fastq.gz`).  

This pipeline uses bowtie2 to map reads to reference sequences so is not designed to work with long reads.  

The expected names of the fastq files are defined by the parameter `fastq_pattern`, whose default value is defined in nextflow.config as `*_R[12]_*.fastq*`.  This pattern can be overridden on the nextflow command line using the `--fastq_pattern` parameter.

The location of the fastq files is specified by the required `fastq_dir` parameter.  The pipeline will detect all files in the `fastq_dir` whose filename matches `fastq_pattern`.

### Mapping samplesheet 

A mapping samplesheet must be input to the pipeline.  This is a 2-column plain-text tab-delimited file with columns:

1. The first row should contain column names: sampleID referenceFasta separated by a tab.
2. The first column should contain sample IDs, i.e. sample IDs input to the Illumina samplesheet.  These will match the base name of the fastq files.
3. The second column should contain the path to a fasta format file containing reference sequences to be mapped to for each dataset.  Not all datasets in fastq_dir need to have a row in this mapping samplesheet.  

An example of a working mapping samplesheet [can be found here](./mapping_samplesheet.txt)

### Reference sequence files

These consist of fasta format files containing sequences that will be mapped to.  The location (path) of these files should match the corresponding entry in the mapping [samplesheet](#Mapping-samplesheet).  Multiple datasets can map to the same reference sequence file.

## Read preprocessing

This pipeline does not perform any adapter or quality trimming of sequences.  We use our lab's [read_preprocessing pipeline](https://github.com/stenglein-lab/read_preprocessing/tree/main) to do this prior to running this workflow.

## Outputs

Some of the main outputs of the pipeline are:

1. Bams of reads mapped to reference sequences: both one bam per dataset and one bam per dataset/refseq. 
2. A tidy-format file of samtools stats output summarizing mapping output
3. A tidy-format file of samtools depth output summarizing per-base depth for all samples and reference sequences
4. A tidy-format file of samtools coverage output summarizing per-refseq coverage for all samples and reference sequences
5. A tidy-format file of the number of reads mapping to each sequence to sense and anti-sense strands (requires strand-specific RNA library prep/data)
6. A tidy-format file of insert sizes for all samples and reference sequences, parsed from samtools stats output
7. Coverage plots for each dataset/refseq.

These can be used as is, input to downstream processing scripts (e.g. R scripts), or input to further processing nextflow steps by cloning and modifying this pipeline.

### Output directory

Output files are placed by default in a `results` directory (or `test/results` when running with `-profile test`).  The output directory location can be overridden using the `--outdir` parameter.  For instance:

```
nextflow run stenglein-lab/remapping_worfklow ... --outdir some_other_results_directory_name
```

## Configurable parameters

Some of the main configurable parameters include:

#### Directories and paths
- fastq_dir: the path to a directory containing fastq 
- fastq_pattern: a pattern to match [fastq input](#Input-fastq)
- outdir: the path of an output directory that will be created and populated with results.
- mapping_sample_sheet: the path to a [mapping samplesheet](#Mapping-samplesheet)

#### Parameters to run optional pipeline steps
- quantify_strand_bias: a boolean to turn on or off strand bias quantification [default: true].  Turn off if you don't have strand-specific data.
- R1_antisense_orientation: a boolean indicating that read 1s are in the antisense orientation relative to the original RNA molecule for stand-specific RNA libraries [default: true]
- skip_per_base_coverage: a boolean that can be turned on to turn off per-base coverage output [default: false].  
- tabulate_insert_sizes: a boolean to turn on tabulation of insert sizes [default: true].  Calculating insert sizes requires paired-end data.  Turn off if you are analyzing single-end data. 
- generate_consensus_sequences: a boolean to turn on generation of new consensus sequences from mapped reads [default: off].

See [nextflow.config](nextflow.config) for more information.

## Software dependencies

This pipeline has two main dependencies: nextflow and singularity.  These programs must be installed on your computer to run this pipeline.

### Nextflow

To run the pipeline you will need to be working on a computer that has nextflow installed. Installation instructions can be [found here](https://www.nextflow.io/docs/latest/getstarted.html#installation).  To test if you have nextflow installed, run:

```
nextflow -version
```

This pipeline requires nextflow version > 24.10

### Singularity

The pipeline uses singularity containers to run programs like cutadapt, BLAST, and R.  To use these containers you must be running the pipeline on a computer that has [singularity](https://sylabs.io/singularity) [installed](https://sylabs.io/guides/latest/admin-guide/installation.html).  To test if singularity is installed, run:

```
singularity --version
```

There is no specified minimum version of singularity, but older versions of singularity (<~3.9) may not work.  the pipeline has been tested with singularity v3.9.5.

Singularity containers will be automatically downloaded and stored in a directory named `singularity_cacheDir` in your home directory.  They will only be downloaded once.

### Conda

It is possible to run this pipeline using [conda](https://docs.conda.io/en/latest/) to handle dependencies. But it is strongly recommended to use singularity instead of conda.

### R libraries

Some of the pipeline code is implemented in [R scripts](./bin/).  Some of these scripts require R packages like [patchwork](https://patchwork.data-imaginist.com/) for creating plot output.  When running with singularity, these packages are installed locally, on top of a [Rocker tidyverse singularity image](https://rocker-project.org/images/).  This occurs in this [nextflow process](modules/stenglein-lab/setup_R_dependencies/main.nf)

