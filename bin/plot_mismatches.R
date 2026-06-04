#!/usr/bin/env Rscript

# This code block sets up input arguments to either come from the command line
# (if running from the pipeline, so not interactively) or to use expected default values
# (if running in interactive mode, for example in RStudio for troubleshooting
# or development).
#
if (!interactive()) {
  # if running from Rscript
  args = commandArgs(trailingOnly=TRUE)
  # lib_dir=args[1]
  mismatch_input        = args[1]
  R_lib_dir             = args[2]
  output_dir            = "./"
} else {
  # if running via RStudio
  mismatch_input        = "../results/mapping/mismatches/collected_mismatches.txt"
  R_lib_dir             = NA
  output_dir            = "../results/plot/"
}

# this library will be available in the tidyverse singularity image we are using 
# (or analogous conda env)
library(tidyverse)

# these libraries are not part of the standard tidyverse, so may have to load it
# from a specified path
# either from pipeline's R lib dir or from R environment
if (!is.na(R_lib_dir)) {
  library(rstatix, lib.loc=R_lib_dir)
  library(ggpubr, lib.loc=R_lib_dir)
  library(patchwork, lib.loc=R_lib_dir)

} else {
  # in this case assuming these will be installed
  library(rstatix)
  library(ggpubr)
  library(patchwork)
}

debug <- 0
if (debug) {
  mismatch_input        = "../results/mapping/mismatches/cm"
}

mismatches <- read.delim(mismatch_input, header=F, sep="\t")

# damage profiler columns in misincorporation.txt output file: 
# Chr     End     Std     Pos     A       C       G       T       Total   
# G>A     C>T     A>G     T>C     A>C     A>T     C>G     C>A    
# T>G     T>A     G>C     G>T     
# A>-     T>-     C>-     G>-     
# ->A     ->T     ->C     ->G     S
colnames(mismatches) <- c(
 "sample_id",
 "refseq",
 "end",
 "strand",
 "position",
 "a_count", "c_count", "g_count", "t_count",
 "total_count",
 "g_a", "c_t", "a_g", "t_c", "a_c", "a_t", "c_g", "c_a", "t_g", "t_a", "g_c", "g_t",  
 "del_a", "del_t", "del_c", "del_g",
 "ins_a", "ins_t", "ins_c", "ins_g", 
 "s")
 
# calculate frequencies
mismatches <- 
  mismatches %>% mutate(
    g_a_f = g_a / total_count,
    c_t_f = c_t / total_count,
    a_g_f = a_g / total_count,
    t_c_f = t_c / total_count,
    a_c_f = a_c / total_count,
    a_t_f = a_t / total_count,
    c_g_f = c_g / total_count,
    c_a_f = c_a / total_count,
    t_g_f = t_g / total_count,
    t_a_f = t_a / total_count,
    g_c_f = g_c / total_count,
    g_t_f = g_t / total_count,
    del_a_f = del_a / total_count,
    del_t_f = del_t / total_count,
    del_c_f = del_c / total_count,
    del_g_f = del_g / total_count,
    ins_a_f = ins_a / total_count,
    ins_t_f = ins_t / total_count,
    ins_c_f = ins_c / total_count,
    ins_g_f = ins_g / total_count
  )

mismatches_long <- 
  mismatches %>% 
  select( sample_id,
          refseq,
          end,
          g_a_f,
          c_t_f,
          a_g_f,
          t_c_f,
          a_c_f,
          a_t_f,
          c_g_f,
          c_a_f,
          t_g_f,
          t_a_f,
          g_c_f,
          g_t_f) %>%
  pivot_longer(cols=ends_with("_f"), names_to = "substitution", values_to = "substitution_frequency")

ggplot(mismatches_long) +
  geom_boxplot(aes(x=substitution, y = substitution_frequency)) +
  facet_grid(sample_id~end) + 
  scale_y_log10() + 
  theme_bw()

ggplot(mismatches) +
  geom_point(aes(x=position, y = g_a_f)) +
  facet_grid(sample_id~end) + 
  xlim(c(0,25)) +
  theme_bw()

ggplot(mismatches) +
  geom_boxplot(aes(x=position, y = g_a_f)) +
  facet_grid(sample_id~end) + 
  xlim(c(0,25)) +
  theme_bw()

plot_all_refseqs <- function(dataset_refseqs){
  
  # TODO: make this configurable via a command-line arg
  plots_per_page <- 8
  
  page_number <- 1
  
  pdf_list <- c()
  
  # iterate through the datasets, doing up to datasets_per_page per page
  for (i in seq(1, nrow(dataset_refseqs), plots_per_page)) {
    
    # plots_per_page at a time
    subset_datasets <- dataset_refseqs %>% filter(row_number() >= i & row_number() <= (i+(plots_per_page-1)))
    
    pdf_name <- paste0(output_dir, "/coverage_plot_page_", page_number, ".pdf")
    
    pdf_list <- c(pdf_list, pdf_name)
    
    # generate & print plot
    plot_datasets(subset_datasets, pdf_name, plots_per_page)
    
    page_number = page_number + 1
  }
  
  # return a list of PDF filenames
  pdf_list
}

# a function to create coverage plots for a certain number of datasets
plot_datasets <- function(dataset_refseq_names, pdf_name = "test.pdf", plots_per_page){
  
  # make the plots
  plots <- apply(dataset_refseq_names, 1, plot_one_refseq)
  
  page_p <- NULL

  # this uses patchwork to add the plot to the page
  for (plot in plots) {
    page_p <- page_p / plot 
  }

  # add a single x axis label at the bottom
  # since all x axes are the same
  page_p <- page_p + xlab("genome position (nt)")
  
  # size of page 
  page_h <- 11
  page_w <- 8.5
  page_u <- "in"

  # size of normal page margins in inches (page_u)
  top_margin    <- 0.25
  left_margin   <- 0.25
  right_margin  <- 0.25
  bottom_margin <- 1

  # add blank spacing so that the size of plots stays consistent on the last page
  # when the last page has fewer plots
  # increase margin if missing plots compared to expected # (e.g. last page)
  missing_plots <- plots_per_page - nrow(dataset_refseq_names) 
  if (missing_plots > 0) {
    extra_margin  <- (missing_plots / plots_per_page) * (page_h - top_margin - bottom_margin)
    bottom_margin <- bottom_margin + extra_margin
  }
  
  # add margins, leaving space at bottom for fig legends
  page_p <- page_p + plot_annotation(theme = theme(plot.margin = margin(t = top_margin, r = right_margin, l = left_margin, b = bottom_margin, unit = page_u)))
  
  # output plot to console
  # print(page_p)

  # save a 1-page PDF
  ggsave(pdf_name, page_p, height=page_h, width=page_w, units=page_u)
}

# a function to create on coverage plot
plot_one_refseq <- function (dataset_refseq_names, number_plot_cols = 1) {
  
  dataset_to_plot <- dataset_refseq_names[1]
  refseq_to_plot  <- dataset_refseq_names[2] 
  
  debug <- 0
  if (debug) {
    dataset_to_plot <- "1004277"
    refseq_to_plot  <- "OR820562_Galbut_virus_RNA_1"
  }

  # subset the main dataframes to get the data just for this dataset/reference_sequence
  subset_df <- df_windowed %>% 
    filter(dataset  == dataset_to_plot & reference_sequence == refseq_to_plot)
  
  # convert any depth of 0 to 1, since plotting on a log10 y scale...
  # subset_df <- subset_df %>% mutate(depth = if_else(depth == 0, 1, depth))
  subset_df <- subset_df %>% mutate(depth = if_else(depth == 0, NA, depth))
  
  # select only necessary columns
  subset_df <- subset_df %>% select(reference_sequence, depth, position)
  
  p <- ggplot(subset_df) + 
    geom_line(aes(x=position, y=depth), linewidth=0.5) +
    geom_ribbon(aes(x=position, ymin = 1, ymax = depth), fill="lightgrey", alpha=0.5) +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    scale_y_log10(limits=c(1,NA)) +
    xlab("") +
    ylab ("coverage depth") +
    ggtitle(NULL, subtitle = paste0(dataset_to_plot, "-", refseq_to_plot)) +
    theme(strip.text.y = element_text(angle = 0)) 
    
  # return the plot
  p
}

# plot all reference sequences
# plot_all_refseqs(dataset_refseqs)
#
