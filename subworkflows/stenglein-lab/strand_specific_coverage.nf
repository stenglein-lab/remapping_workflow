/*
 Use bamtocov utility to tabulate strand-specific depth of coverage
 */
workflow STRAND_SPECIFIC_COVERAGE {

  take:
  bam                   // channel [val (meta), path (bam)]

  main:
  ch_versions = Channel.empty()

  BAM_TO_COV(bam)
  ch_versions.mix(BAM_TO_COV.out.versions)

  BED_TO_PER_BASE_COV(BAM_TO_COV.out.bed_coverage)
  ch_versions.mix(BED_TO_PER_BASE_COV.out.versions)

  emit: 
  coverage_report          = BAM_TO_COV.out.coverage_report             // channel [val(meta), path(coverage_report)
  strand_specific_coverage = BED_TO_PER_BASE_COV.out.per_base_coverage  // channel [val(meta), path(per_base_coverage)
  bed_coverage             = BAM_TO_COV.out.bed_coverage                // channel [val(meta), path(bamtocov_bed_output)
  versions                 = ch_versions                                // channel path(versions.yml)

}

/*
  This process uses a perl scrip tto convert bamtocov BED coverage output into tidy per-base coverage output
 */
process BED_TO_PER_BASE_COV {
  tag "$meta.id"
  label 'lowmem_threaded'

  conda (params.enable_conda ? 'conda-forge::perl=5.32.1' : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "https://depot.galaxyproject.org/singularity/perl:5.26.2"
  } else {
    container "quay.io/biocontainers/perl:5.26.2"
  }

  input:
  tuple val(meta), path (bed)

  output:
  tuple val (meta), path ("*.per_base.txt")         , emit: per_base_coverage
  path "versions.yml"                               , emit: versions            


  when:
  task.ext.when == null || task.ext.when

  script:
  def args             = task.ext.args ?: ''
  """
  # perl script in workflow bin dir
  stranded_bed_to_per_base_coverage $bed > ${bed}.per_base.txt

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
  END_VERSIONS
  """

}

process BAM_TO_COV {
  tag "$meta.id"
  label 'lowmem_threaded'

  conda (params.enable_conda ? 'bioconda::bamtocov=2.7.0' : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "https://depot.galaxyproject.org/singularity/bamtocov:2.7.0--hbd632db_1"
  } else {
    container "quay.io/biocontainers/bamtocov:2.7.0--hbd632db_1"
  }

  input:
  tuple val(meta), path (bam)

  output:
  tuple val (meta), path ("*.coverage_report.txt")  , emit: coverage_report
  tuple val (meta), path ("*.coverage.txt")         , emit: bed_coverage
  path "versions.yml"                               , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args             = task.ext.args ?: ''

  """
  bamtocov --stranded --threads $task.cpus -o ${meta.id}.coverage_report.txt $bam  > ${meta.id}.coverage.txt

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      bamtocov: 2.7.0 
  END_VERSIONS
  """

}
