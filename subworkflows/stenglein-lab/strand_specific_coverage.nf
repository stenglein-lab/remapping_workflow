/*
 tabulate strand-specific depth of coverage
 */
workflow STRAND_SPECIFIC_COVERAGE {

  take:
    bam                   // channel [val (meta), path (bam)]
    min_mapq

  main:
    BAM_TO_STRAND_SPECIFIC_COVERAGE(bam, min_mapq)

  emit: 
    strand_specific_coverage = BAM_TO_STRAND_SPECIFIC_COVERAGE.out.strand_specific_coverage  

}

/*
 The invoked python script uses pysam to tabulate and output strand-specific coverage at all positions of all
 refseqs in bam.
 */
process BAM_TO_STRAND_SPECIFIC_COVERAGE {
  label 'lowmem_non_threaded'

  // singularity info for this process
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.0--py314ha0fe93d_0' :
        'biocontainers/pysam:0.24.0--py314ha0fe93d_0' }"

  input:
  tuple val(meta), path(bam)
  val(min_mapping_quality)

  output:
  tuple val(meta), path ("*.strand_specific_coverage.txt") , emit: strand_specific_coverage

  when:
  task.ext.when == null || task.ext.when

  script:
  """
  bam_to_strand_specific_coverage_depth.py \
     --prefix ${meta.id} \
     --min_mapq ${min_mapping_quality} \
     $bam \
   > ${meta.id}.strand_specific_coverage.txt
  """

}

