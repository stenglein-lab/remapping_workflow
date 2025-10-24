library(tidyverse)

# read in data
df <- read.delim("virus_refseq/per_sample_segment_counts.txt", header=F, sep="\t")
colnames(df) <- c("sample_id", "segment", "number_genotypes")

ggplot(df) +
  geom_tile(aes(y=sample_id, x=segment, fill=number_genotypes)) +
  scale_fill_continuous (low="white", high="black")  +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("") +
  ylab("")

ggsave("segments_per_sample.pdf", units="in", width=7.5, height=10)
