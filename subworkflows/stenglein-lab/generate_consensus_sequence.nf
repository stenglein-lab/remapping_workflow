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

 emit:
  versions          = ch_versions
  refseq_and_new    = VIRAL_CONSENSUS.out.refseq_and_new
  new_fasta         = VIRAL_CONSENSUS.out.fasta
  position_counts   = VIRAL_CONSENSUS.out.position_counts
  refseq            = VIRAL_CONSENSUS.out.refseq

}
