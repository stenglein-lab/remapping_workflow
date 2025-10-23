# remapping_workflow
This is a nextflow workflow for remapping reads to sets of virus sequences.  In our lab we use it to validate virus sequences identified by metagenomic sequencing.  It could be used as a stand-alone pipeline or as a subworkflow.

### Software dependencies

These analyses are implemented in [nextflow](https://www.nextflow.io/docs/latest/).  Dependencies are handled through singularity so installation of other software besides nextflow and singularity shouldn't be necessary.

### virus_evol_paper branch

This branch contains the version of this pipeline that was used to remap adapter- and quality-trimmed reads to recovered galbut virus reference sequences for the paper "Constrained Reassortment and Genotype-Specific Traits Shape the Evolutionary Landscape of Galbut Virus" by Keene-Snickers et al published in 2025 in Virus Evolution.

