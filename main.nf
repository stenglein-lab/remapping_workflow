#!/usr/bin/env nextflow

include { REMAPPING_WORKFLOW   } from './subworkflows/stenglein-lab/remapping_workflow'

workflow {

  main: 
    REMAPPING_WORKFLOW ()

  publish:
     bam            = REMAPPING_WORKFLOW.out.bam
     bowtie2_log    = REMAPPING_WORKFLOW.out.bowtie2_log
     consensus      = REMAPPING_WORKFLOW.out.consensus
     coverage       = REMAPPING_WORKFLOW.out.coverage
     coverage_plots = REMAPPING_WORKFLOW.out.coverage_plots
     depth          = REMAPPING_WORKFLOW.out.depth
     insert_sizes   = REMAPPING_WORKFLOW.out.insert_sizes
     samples        = REMAPPING_WORKFLOW.out.samples
     stats          = REMAPPING_WORKFLOW.out.stats
     strand_bias    = REMAPPING_WORKFLOW.out.strand_bias
    
}

// specifye where main output files will go
output {
    bam {
        path 'mapping/bam'
        mode 'link'
    }
    bowtie2_log {
        path 'mapping/mapping_logs'
        mode 'link'
    }
    consensus {
        path 'consensus_sequences'
        mode 'link'
    }
    coverage {
        path 'mapping/coverage'
        mode 'link'
    }
    coverage_plots {
        path 'coverage_plots'
        mode 'link'
    }
    depth {
        path 'mapping/depth'
        mode 'link'
    }
    insert_sizes {
        path 'mapping/insert_sizes'
        mode 'link'
    }
    samples {
        path 'sample_info'
        index { path 'sample_info.json' }
    }
    stats {
        path 'mapping/stats'
        mode 'link'
    }
    strand_bias {
        path 'mapping/strand_bias'
        mode 'link'
    }
}
