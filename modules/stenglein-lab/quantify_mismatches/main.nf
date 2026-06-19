process QUANTIFY_MISMATCHES {
  label 'lowmem_non_threaded'

  // singularity info for this process
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.0--py314ha0fe93d_0' :
        'biocontainers/pysam:0.24.0--py314ha0fe93d_0' }"

  input:
  tuple val(meta), path(bam), val(refseq), path (refseq_fasta), path (fai), path(bai)
  val(min_depth)
  val(min_base_quality)
  val(min_mapping_quality)

  output:
  tuple val(meta), path ("*.mismatches.txt")  , emit: txt
  // path "versions.yml"  , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:

  def args             = task.ext.args ?: ''

  """
  # pipe to awk to prepend sample ID
  tabulate_mismatches_from_bam.py \
     $args \
     --ref $refseq_fasta \
     --bam $bam \
     --min-depth $min_depth \
     --min-base-quality $min_base_quality \
     --min-mapping-quality $min_mapping_quality \
  | awk '{print "${meta.id}" "\t" \$0}' > ${meta.id}.${refseq.id}.mismatches.txt
  """

}

