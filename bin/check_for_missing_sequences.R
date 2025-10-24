library(tidyverse)

# double check that we've accounted for all galbut virus positive samples

if (!interactive()) {
  # if running from Rscript
  args = commandArgs(trailingOnly=TRUE)
  mapping_sample_sheet=args[1]
  galbut_tallies_files=args[2]
  outdir="./"
} else {
  # if running via RStudio
  mapping_sample_sheet="../mapping_samplesheet.txt"
  galbut_tallies_files="../../metagenomic/results/taxa_matrices/virus_matrix.txt"
  outdir="."
}

# read in mapping sample sheet for galbut virus sequences Lexi has recovered
samples <- read.delim(mapping_sample_sheet, sep="\t", header=T)
colnames(samples) <- c("sample_id", "refseq_fasta")

# read in matrix of virus-mapping reads
virus_matrix <- read.delim(galbut_tallies_files, sep="\t", skip=1, header=T)
virus_matrix <- virus_matrix %>% rename(filename = "X") 
virus_matrix <- 
  virus_matrix %>% 
  mutate(sample_id = str_replace(filename, ".trim_contigs_and_singletons.fasta.bn_nt.tally", "")) %>%
  relocate(sample_id) %>%
  select(-filename)

# of galbut virus-mapping reads per dataset
galbut_reads <- virus_matrix %>% select(sample_id, Galbut.virus)

# datasets with >100 galbut-mapping reads
galbut_gt_100 <- filter(galbut_reads, Galbut.virus > 100)

# how many datasets with >100 galbut-mapping reads and with galbut sequences?
paste0("There are ", 
       nrow(galbut_gt_100),
       " datasets with >100 galbut-mapping reads, and ",
       nrow(samples), 
       " datasets with galbut virus sequences.")

# But are they the same?  
# check if there are any datasets with > 100 galbut-mapping reads that Lexi hasn't accounted for
for (sample_id in galbut_gt_100$sample_id) {
  if (!sample_id %in% samples$sample_id){
    error_msg <- paste0("Error: sample ", sample_id, "has > 100 galbut-mapping reads ",
                   "but no galbut sequences defined in mapping sample sheet")
    stop(error_msg)
  }
  else {
    message(paste0("sample ", sample_id, "has > 100 galbut-mapping reads and is accounted for."))
  }
}

