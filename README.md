# remapping_workflow
This is a nextflow workflow for remapping reads to sets of virus sequences.  In our lab we use it to validate virus sequences identified by metagenomic sequencing.  It could be used as a stand-alone pipeline or as a subworkflow.

### Software dependencies

These analyses are implemented in [nextflow](https://www.nextflow.io/docs/latest/).  Dependencies are handled through singularity so installation of other software besides nextflow and singularity shouldn't be necessary.


