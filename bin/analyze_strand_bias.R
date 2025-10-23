library(tidyverse)
library(rstatix)
library(patchwork)
library(ggpubr)

strand_bias <- read.delim("../results/save/collected_strand_bias.tsv", sep="\t", header=F)

sample_accession_map <- read.delim("../virus_refseq/sample_accession_map.txt", sep="\t", header=F)
colnames(sample_accession_map) <- c("sample_id", "segment", "accession", "full_name")


colnames(strand_bias) <- c("sample_id", 
                           "sam_file", 
                           "accession", 
                           "plus_strand_mapping_reads", 
                           "minus_strand_mapping_reads",
                           "plus_strand_fraction")

# merge in accession metadata
strand_bias <- left_join(strand_bias, sample_accession_map)

# calculate +strand copies per -strand copy
strand_bias <- strand_bias %>% mutate(plus_strand_copies_per_minus_strand = plus_strand_mapping_reads / minus_strand_mapping_reads)

# number of host-mapping reads in each dataset
# read in host mapping read counts
host_mapping_counts <- read.delim("../../DiverseCollections/analyses/data/host_mapping_read_counts.txt", sep="\t", header=F)
colnames(host_mapping_counts) <- c("sample_id", "host_mapping_reads")

# read in RpL32 mRNA mapping reads counts
rpl32_mapping_counts <- read.delim("../../DiverseCollections/analyses/data/RpL32_mapping_reads_per_dataset.txt", sep="\t", header=F)
# output format from samtools coverage
colnames(rpl32_mapping_counts) <- c("id", "accession", "startpos", "endpos", 
                                    "numreads", "covbases", "coverage", 
                                    "meandepth",  "meanbaseq", "meanmapq")

# simplify
rpl32_mapping_counts <- rpl32_mapping_counts %>% select(id, numreads) %>% rename(rpl32_mapping_counts = numreads)

# merge in host-mapping read counts
strand_bias <- left_join(strand_bias, host_mapping_counts, by = )

# normalize read-mapping counts to host-mapping reads counts
strand_bias <- strand_bias %>% mutate(plus_strand_rpm  = 1e6 * plus_strand_mapping_reads  / host_mapping_reads,
                                      minus_strand_rpm = 1e6 * minus_strand_mapping_reads / host_mapping_reads,
                                      log_plus_strand_rpm  = log10(plus_strand_rpm), 
                                      log_minus_strand_rpm = log10(minus_strand_rpm)
                                      )

# reorder levels for plotting
strand_bias$segment <- fct_relevel(strand_bias$segment, "Chaq", after=3)

# run stats

# first, test for normality of data
# this Shapiro tests indicates fractional counts are not normally distributed
# TODO: should use log-normalized values?
df_shapiro     <- strand_bias %>% group_by(segment) %>% shapiro_test(log_plus_strand_rpm)
filter(df_shapiro, p<0.05)
df_shapiro     <- strand_bias %>% group_by(segment) %>% shapiro_test(minus_strand_rpm)
filter(df_shapiro, p<0.05)


# run wilcoxon tests 
wilcox_plus_strand_fraction <- strand_bias %>% wilcox_test(plus_strand_fraction ~ segment)
wilcox_plus_strand_rpm      <- strand_bias %>% wilcox_test(plus_strand_rpm ~ segment)
wilcox_minus_strand_rpm     <- strand_bias %>% wilcox_test(minus_strand_rpm ~ segment)
wilcox_plus_per_minus       <- strand_bias %>% wilcox_test(plus_strand_copies_per_minus_strand ~ segment)

wilcox_plus_strand_fraction_p <- wilcox_plus_strand_fraction %>% add_xy_position(x="segment")
wilcox_plus_strand_rpm_p      <- wilcox_plus_strand_rpm %>% add_xy_position(x="segment")
wilcox_minus_strand_rpm_p     <- wilcox_minus_strand_rpm %>% add_xy_position(x="segment")
wilcox_plus_per_minus_p       <- wilcox_plus_per_minus %>% add_xy_position(x="segment")


shared_theme <- function()
{
  theme_bw(base_size = 11) 
}

# plot fraction +strand
fraction_plus_strand_p <- ggplot(strand_bias) +
  geom_jitter(aes(x=segment, y=plus_strand_fraction), 
              shape=21, fill="slateblue", color="black", stroke=0.25,
              height=0, width=0.25, alpha=0.5) +
  geom_boxplot(aes(x=segment, y=plus_strand_fraction), 
               fill=NA, width=0.5, outlier.shape = NA,
               color="slateblue", size = 0.5) +
  scale_y_continuous(limits=c(0,1.15), breaks=c(0, 0.5, 1)) +
  ylab("Fraction plus-strand\nmapping reads") +
  xlab("") +
  shared_theme() 

fraction_plus_strand_p <- 
  fraction_plus_strand_p + 
  stat_pvalue_manual(filter(wilcox_plus_strand_fraction_p, p.adj.signif != "ns"), 
                     y.position = c(1.05, 1.10, 1.15, 1.05),
                     tip.length = 0.01, bracket.shorten=0.1, size=3, label="p.adj.signif")

fraction_plus_strand_p

# plus copies per minus strand
plus_copies_per_minus_p <- ggplot(strand_bias) +
  geom_jitter(aes(x=segment, y=plus_strand_copies_per_minus_strand), 
              shape=21, fill="slateblue", color="black", stroke=0.25,
              height=0, width=0.25, alpha=0.5) +
  geom_boxplot(aes(x=segment, y=plus_strand_copies_per_minus_strand), 
               fill=NA, width=0.5, outlier.shape = NA,
               color="slateblue", size = 0.5) +
  scale_y_log10() + 
  xlab("") +
  ylab("Plus strand mapping reads\nper minus strand mapping read") +
  shared_theme()

plus_copies_per_minus_p <-
  plus_copies_per_minus_p +
  stat_pvalue_manual(filter(wilcox_plus_per_minus_p, p.adj.signif != "ns"), 
                     y.position = log10(c(700, 900, 1100, 700, 700)),
                     tip.length = 0.01, bracket.shorten=0.1, size=3, label="p.adj.signif")

plus_copies_per_minus_p

# plot minus strand mapping reads
minus_strand_p <- ggplot(strand_bias) +
  geom_jitter(aes(x=segment, y=minus_strand_mapping_reads), 
              shape=21, fill="slateblue", color="black", stroke=0.25,
              height=0, width=0.25, alpha=0.5) +
  geom_boxplot(aes(x=segment, y=minus_strand_mapping_reads), 
               fill=NA, width=0.5, outlier.shape = NA,
               color="slateblue", size = 0.5) +
  scale_y_log10() +
  xlab("") +
  ylab("Minus strand mapping reads\nper million host-mapping reads") +
  shared_theme() 

# all non-significant
minus_strand_p

# plot plus strand mapping reads
plus_strand_p <- ggplot(strand_bias) +
  geom_jitter(aes(x=segment, y=plus_strand_mapping_reads), 
              shape=21, fill="slateblue", color="black", stroke=0.25,
              height=0, width=0.25, alpha=0.5) +
  geom_boxplot(aes(x=segment, y=plus_strand_mapping_reads), 
               fill=NA, width=0.5, outlier.shape = NA,
               color="slateblue", size = 0.5) +
  scale_y_log10() +
  xlab("") +
  ylab("Plus strand mapping reads\nper million host-mapping reads") +
  shared_theme() 

plus_strand_p

combined_p <- 
  fraction_plus_strand_p + 
  plus_strand_p + 
  minus_strand_p + 
  plus_copies_per_minus_p + 
  plot_layout(ncol = 1) + 
  plot_annotation(tag_levels = "A")

combined_p

ggsave("../results/Supplemental_figure_6_segment_strand_levels.pdf", units="in", width=7.5, height=10)

# calculate averages
df_avg <- strand_bias %>% 
  group_by(segment) %>%
  summarize(median_plus      = median(plus_strand_rpm, na.rm = T),
            median_minus     = median(minus_strand_rpm, na.rm = T),
            median_ppm       = median(plus_strand_copies_per_minus_strand, na.rm = T),
            median_frac_plus = median(plus_strand_fraction, na.rm = T))

# output text for paper

paste0("RNA 1 was the most common coding-incomplete sequence (Fig. X). ",
       "We had prepared libaries using a strand-specific protocol, ",
       "so could assess levels of plus (+) and minus (-) strand mapping reads separately (Supp. Fig. X). ",
       "Galbut virus genomic RNA is presumed to be double stranded, ",
       "so levels of -strand RNA can be used as a proxy for genomic RNA levels, ",
       "whereas +strand RNA reflect levels of genomic RNA and mRNA. ",
       "Levels of RNA 1 -strand RNA were not significantly lower than -strand levels of other segments (Supp. Fig. X). ",
       "However, RNA 1 exhibited the lowest ratio of +strand RNA copies per -strand RNA (Supp. Fig. X). ",
       "This is consistent with decreased transcription of RNA 1 or a shorter half life of RNA 1 transcripts, ",
       " which contributed to recovery of fewer coding complete RNA 1 sequences. "
)
       
and exhibited lower ",
       "average coverage levels than RNA 2, RNA 3, or chaq virus (Supp. Fig. 6). ",
       "The median level of RNA 1 mapping reads per million host-mapping reads (RPM) was ",
       sprintf("%0.0f", filter(df_avg, segment == "RNA1") %>% pull(median_plus)),
       ", which was s"
       
       )


       # 193, which was significantly lower than the median RPM values for RNA 2 and RNA 3 of 428 and 448 respectively (Wilcoxon p adjusted = 0.025 and 0.026; Supp. Fig. 6).
       # 193, which was significantly lower than the median RPM values for RNA 2 and RNA 3 of 428 and 448 respectively (Wilcoxon p adjusted = 0.025 and 0.026; Supp. Fig. 6).

