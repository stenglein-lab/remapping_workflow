process QUANTIFY_MISMATCHES {
  label 'process_medium'

  // singularity info for this process
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.0--py314ha0fe93d_0' :
        'biocontainers/pysam:0.24.0--py314ha0fe93d_0' }"

  input:
  tuple val(meta), path(bam), val(refseq), path (refseq_fasta), path (fai), path(bai)
  val(min_depth)
  val(max_depth)
  val(min_base_quality)
  val(min_mapping_quality)

  output:
  tuple val(meta), path ("*mismatches.txt")             , emit: mismatches
  tuple val(meta), path ("*.mismatches.txt")             , emit: per_refseq_per_position_mismatches
  tuple val(meta), path ("*.per_refseq_mismatches.txt")  , emit: per_refseq_mismatches
  tuple val(meta), path ("*.by_read_position_mismatches.txt") , emit: by_position_mismatches
  tuple val(meta), path ("*.total_mismatches.txt")       , emit: total_mismatches
  // path "versions.yml"  , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:

  def args             = task.ext.args ?: ''

  """
  tabulate_mismatches_from_bam.py \
     $args \
     --ref $refseq_fasta \
     --bam $bam \
     --prefix ${meta.id} \
     --mismatch_types_per_ref_out ${meta.id}.${refseq.id}.per_refseq_mismatches.txt \
     --mismatch_types_by_pos ${meta.id}.${refseq.id}.by_read_position_mismatches.txt \
     --mismatch_types_total_out ${meta.id}.${refseq.id}.total_mismatches.txt \
     --min_depth $min_depth \
     --max_depth $max_depth \
     --min_base_quality $min_base_quality \
     --min_mapping_quality $min_mapping_quality \
   > ${meta.id}.${refseq.id}.mismatches.txt
  """

}

