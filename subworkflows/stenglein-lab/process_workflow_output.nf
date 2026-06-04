include { SETUP_R_DEPENDENCIES   } from '../../modules/stenglein-lab/setup_R_dependencies'
include { PLOT_REFSEQ_COVERAGE   } from '../../modules/stenglein-lab/plot_refseq_coverage'
include { PLOT_MISMATCHES        } from '../../modules/stenglein-lab/plot_mismatches'
                                                                                
workflow PROCESS_WORKFLOW_OUTPUT {
                                                                                
  take:
   // these are all paths to output or metadata files
   coverage     
   stats        
   depth        
   insert_sizes 
   mismatches 
   strand_bias  

  main:                                                                         

   // setup (install) some R packages that will be needed on top of tidyverse
   SETUP_R_DEPENDENCIES(params.R_packages)

   PLOT_REFSEQ_COVERAGE(depth, SETUP_R_DEPENDENCIES.out.R_lib_dir)

   if (params.quantify_mismatches) {
      PLOT_MISMATCHES(mismatches, SETUP_R_DEPENDENCIES.out.R_lib_dir)
   }

  emit:

   coverage_plots = PLOT_REFSEQ_COVERAGE.out.pdf
   mismatch_plots = PLOT_MISMATCHES.out.pdf
 
}         
