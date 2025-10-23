#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { REMAPPING_WORKFLOW   } from './subworkflows/stenglein-lab/remapping_workflow'
include { QUANTIFY_STRAND_BIAS } from './subworkflows/stenglein-lab/quantify_strand_bias'

// main named workflow
workflow MAIN_WORKFLOW {

    // main remapping workflow
    REMAPPING_WORKFLOW ()

    // run optional workflow to quantify strand bias
    if (params.quantify_strand_bias) {
       QUANTIFY_STRAND_BIAS(REMAPPING_WORKFLOW.out.sam)
    }
}

//
// entry workflow
// https://www.nextflow.io/docs/latest/reference/syntax.html#workflow
//
workflow {
    MAIN_WORKFLOW ()
}

