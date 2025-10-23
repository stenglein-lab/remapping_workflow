library (tidyverse)
library (openxlsx)
library (pdftools)
library (patchwork)

depth_df <- read.delim("../results/save/collected_per_base_depth.tsv", header=F, sep="\t")

colnames(depth_df) <- c("dataset", "reference_sequence", "position", "depth")

# read in dataset sizes file
total_reads <- read.delim("../../DiverseCollections/analyses/data/all_read_counts.txt")
total_reads <- total_reads %>% mutate (dataset = str_replace(sample_id, ".trim", ""))
total_reads <- total_reads %>% filter(count_type == "post_collapse") %>% select(dataset, count) %>% rename(total_reads = count)

# read in map of refseq accessions to segment, etc
segment_map <- read.delim("../../DiverseCollections/analyses/data/sample_accession_map.txt", header=F)
colnames(segment_map) <- c("dataset", "segment", "reference_sequence", "description")
segment_map <- segment_map %>% select(-description)
# segment_map$segment <- fct_relevel(segment_map$segment, "RNA1", "RNA2", "RNA3", "Chaq")

depth_df <- left_join(depth_df, segment_map, by = join_by(dataset, reference_sequence))
# depth_df$segment <- fct_relevel(depth_df$segment, "RNA1", "RNA2", "RNA3", "Chaq")

# create new variable - combined label for plotting
depth_df <- depth_df %>% mutate(acc_seg = paste0(segment, "_", reference_sequence))

# calculated median depth of (total) coverage for each ref seq in each dataset and store it in a new df
median_depths <- depth_df %>% 
  group_by(dataset, reference_sequence) %>% 
  summarize(median_depth = median(depth),
            mean_depth = mean(depth))

# what are lowest coverage samples?
median_depths %>% arrange (median_depth)

# output text related to coverage levels for paper
min_median    <- min(median_depths$median_depth)
max_median    <- max(median_depths$median_depth)
median_median <- median(median_depths$median_depth)

paste0(
  "The median coverage depth across all segments was ",
  sprintf("%0.0fx ", median_median),
  "(range: ",
  sprintf("%0.0fx", min_median),
  "-",
  sprintf("%0.0fx", max_median),
  ")."
  )

# ----------------------------------
# calculate average depth in windows
# ----------------------------------
# %/% is the integer division operator
window_size = 10
depth_df <- depth_df %>% mutate (window = position %/% window_size)

# calculate average coverage depth in each window
df_windowed <- depth_df %>% 
  group_by(dataset, reference_sequence, window, acc_seg, segment)  %>% 
  summarize(depth = mean(depth), .groups = "drop") %>% 
  mutate(position = (window*window_size) + 1) %>% 
  ungroup()

##now plot coverage data on multiple pdf pages

# these are the datasets
datasets <- depth_df %>% group_by(dataset) %>% summarize(.groups="drop") %>% pull()

plot_all_datasets <- function(datasets){
  
  datasets_per_page <- 8
  
  page_number <- 1
  
  pdf_list <- c()
  
  # iterate through the datasets, doing up to 12 per page
  for (i in seq(1, length(datasets), datasets_per_page)) {
    
    # plots_per_page at a time
    subset_datasets <- datasets[i:(i+(datasets_per_page-1))]
    
    # if we've gone out of bounds would have NAs, get rid of them
    subset_datasets <- subset_datasets[!is.na(subset_datasets)]
    
    # output to console which ones we're doing
    # print(paste0(subset_datasets))
    
    pdf_name <- paste0("../results/coverage_plot_page_", page_number, ".pdf")
    
    pdf_list <- c(pdf_list, pdf_name)
    
    # generate & print plot
    plot_datasets(subset_datasets, pdf_name)
    
    page_number = page_number + 1
  }
  
  # return a list of PDF filenames
  pdf_list
}

plot_datasets <- function(dataset_names, pdf_name = "test.pdf"){
  
  # dataset_names <- c("1020-A7", "ME-M-3")
  
  # dataset_names <- head(datasets, n=3)
  # dataset_names
  plots <- lapply(dataset_names, plot_one_dataset)
  
  page_p <- NULL
  for (plot in plots) {
    page_p <- page_p / plot 
  }
  page_p <- page_p + xlab("genome position (nt)")
  # page_p
  
  # output plot to console
  print(page_p)
  
  # save a 1-page PDF
  ggsave(pdf_name, page_p, height=10.5, width=7.5, units="in")
}

plot_one_dataset <- function (dataset_to_plot, number_plot_cols = 5) {
  
  # dataset_to_plot <- "ME-M-3"
  # number_plot_cols = 5
  # subset the main dataframes to get the data just for these reference_sequence
  subset_df <- df_windowed %>% filter(dataset  == dataset_to_plot)
  
  # convert any depth of 0 into 1, since plotting on a log10 y scale...
  subset_df <- subset_df %>% mutate(depth = if_else(depth == 0, 1, depth))
  
  # create numeric sorting variable so segments will plot in order
  subset_df <- 
    subset_df %>% mutate(fac_order = case_when (
      segment == "RNA1" ~ 1,
      segment == "RNA2" ~ 2,
      segment == "RNA3" ~ 3,
      .default = 4
  ))
  
  subset_df <- subset_df %>% select(acc_seg, depth, position, fac_order)
  
  acc_segs <- subset_df %>% group_by(acc_seg) %>% summarize() %>% pull()
  
  # make more cols if necessary
  if (length (acc_segs) > number_plot_cols ) { 
    number_plot_cols = length(acc_segs)
  }
  
  # create fake data rows if necessary
  if (length (acc_segs) < number_plot_cols ) { 
    for (space_i in (length(acc_segs)+1):number_plot_cols) {
      new_rows = tibble(acc_seg = space_i, depth = 1, position = c(0,1500), fac_order = space_i * 10)
      subset_df <- rbind(subset_df, new_rows)
    }
  }

  # reorder for plotting
  subset_df$acc_seg <- fct_reorder(subset_df$acc_seg, subset_df$fac_order)
  
  p <- ggplot(subset_df) + 
    geom_line(aes(x=position, y=depth), size=0.5) +
    geom_area(aes(x=position, y=depth), fill="lightgrey", alpha=0.5) +
    scale_color_manual(values = c("red", "black")) +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    scale_y_log10() +
    # xlab ("genome position (nt)") +
    xlab("") +
    ylab ("coverage depth") +
    ggtitle(NULL, subtitle = dataset_to_plot) +
    facet_wrap(~acc_seg, scales="free_x", ncol = number_plot_cols) + 
    theme(strip.text.y = element_text(angle = 0)) 
 
   
  # return the plot
  p
}

# first_30 <- head(datasets, n=30)
# plot_all_datasets(first_30)

pdf_list <- plot_all_datasets(datasets)

# output a combined PDF
# this pdf will have blank facets, which help everything stay aligned
# manually remove these blank facets in Affinity Designer 
pdf_combine(pdf_list, output="../results/Supplemental_Figure_1_coverage_plot.pdf")
