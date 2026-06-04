process PLOT_MISMATCHES {
  label 'process_single'
  tag   "${mismatches}"

  // singularity info for this process
  if (workflow.containerEngine == 'singularity'){
      container "docker://rocker/tidyverse:4.3.2"
  }     

  input:
  path (mismatches)
  path (R_lib_dir)

  output:
  path "*.pdf"                       , emit: pdf, optional: true
  path "*.txt"                       , emit: txt, optional: true
  // path "collected_coverage_plot.pdf" , emit: coverage_plot

  when:
  task.ext.when == null || task.ext.when

  script:

  def args             = task.ext.args ?: ''

  """
   plot_mismatches.R $mismatches $R_lib_dir
  """

}
