include { SAMTOOLS_STATS                         } from '../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_COVERAGE                      } from '../../modules/nf-core/samtools/coverage/main'
include { SAMTOOLS_DEPTH                         } from '../../modules/nf-core/samtools/depth/main'
include { PREPEND_TSV_WITH_ID as PREPEND_STATS_WITH_ID } from '../../modules/stenglein-lab/prepend_tsv_with_id'
include { PREPEND_TSV_WITH_ID as PREPEND_COV_WITH_ID   } from '../../modules/stenglein-lab/prepend_tsv_with_id'
include { PREPEND_TSV_WITH_ID as PREPEND_DEPTH_WITH_ID } from '../../modules/stenglein-lab/prepend_tsv_with_id'

/*
  Calculate basic mapping summary statistics from bam
 */ workflow MAPPING_STATS {

 take:
  bam_fasta       // [meta, bam, fasta]
  per_base_depth  // boolean: calculate per-base coverage depth values using samtools depth?

 main:

  // define some empty channels for keeping track of stuff
  ch_versions     = Channel.empty()                                               

  // split up input channel into separate bam and fasta channels
  bam          = bam_fasta.map{meta, bam, fasta -> [meta, bam]   }
  genome_fasta = bam_fasta.map{meta, bam, fasta -> [meta, fasta] }

  // ------------------
  // samtools stats
  // ------------------

  SAMTOOLS_STATS(bam, genome_fasta)
  ch_versions = ch_versions.mix ( SAMTOOLS_STATS.out.versions )      

  PREPEND_STATS_WITH_ID(SAMTOOLS_STATS.out.stats)

  // ------------------
  // samtools coverage
  // ------------------

  SAMTOOLS_COVERAGE(bam)
  ch_versions = ch_versions.mix ( SAMTOOLS_COVERAGE.out.versions )      

  PREPEND_COV_WITH_ID(SAMTOOLS_COVERAGE.out.coverage)

  // ------------------
  // samtools depth
  // ------------------
  // this outputs per-based coverage depth values
  // only run if requested to do so

  ch_depth           = Channel.empty()
  ch_prepended_depth = Channel.empty()

  if (per_base_depth) {

    SAMTOOLS_DEPTH(bam)
    
    // do these extra mix calls because possibility of empty channel  
    ch_depth           = ch_depth.mix           ( SAMTOOLS_DEPTH.out.tsv )

    PREPEND_DEPTH_WITH_ID(SAMTOOLS_DEPTH.out.tsv)

    ch_prepended_depth = ch_prepended_depth.mix ( PREPEND_DEPTH_WITH_ID.out.tsv)

    ch_versions = ch_versions.mix ( SAMTOOLS_DEPTH.out.versions )      
  }

 emit: 
  versions      = ch_versions
  stats         = SAMTOOLS_STATS.out.stats
  coverage      = SAMTOOLS_COVERAGE.out.coverage
  depth         = ch_depth
  prepended_stats    = PREPEND_STATS_WITH_ID.out.tsv
  prepended_coverage = PREPEND_COV_WITH_ID.out.tsv
  prepended_depth    = ch_prepended_depth

}
