include { VIRAL_CONSENSUS       } from '../../modules/stenglein-lab/viral_consensus/main'

/*
  Generate one new consensus sequence from a refseq and mapped reads
  using the viral_consensus tool
 */ 

workflow GENERATE_CONSENSUS_SEQUENCE {

 take:
  bam_fasta       // [meta, bam, fasta]
  min_depth
  min_qual
  min_freq

 main:

  // define some empty channels for keeping track of stuff
  ch_versions         = Channel.empty()

  // generate consensus 
  VIRAL_CONSENSUS(bam_fasta, min_qual, min_depth, min_freq)
  ch_versions = ch_versions.mix ( VIRAL_CONSENSUS.out.versions )

  // concatenate consensus sequences for each dataset
  // group consensus fasta by meta
  ch_merged_fasta = VIRAL_CONSENSUS.out.fasta.groupTuple() 
  CONCATENATE_CONSENSUS(ch_merged_fasta)

 emit:
  merged_fasta      = CONCATENATE_CONSENSUS.out.fasta
  refseq_and_new    = VIRAL_CONSENSUS.out.refseq_and_new
  new_fasta         = VIRAL_CONSENSUS.out.fasta
  position_counts   = VIRAL_CONSENSUS.out.position_counts
  refseq            = VIRAL_CONSENSUS.out.refseq
  versions          = ch_versions

}

process CONCATENATE_CONSENSUS {
  tag "$meta.id"

  // just a basic container that includes cat command
  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/viral_consensus:1.0.1--hcf1f8c1_0':
    'biocontainers/viral_consensus:1.0.1--hcf1f8c1_0' }"

  input:
  tuple val(meta), path(fastas)

  output:
  tuple val(meta), path("*.consensus.fasta", includeInputs: false)                   , emit: fasta

  when:
  task.ext.when == null || task.ext.when

  script:
  """
  cat ${fastas} > "${meta.id}.consensus.fasta"
  """
}


