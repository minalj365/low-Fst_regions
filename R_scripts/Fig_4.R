rm(list=ls())
library(gggenomes)
library(tidyverse)
library(tidyverse)
library(ggthemes)
library(cowplot)
library(zoo)
library(purrr)
library(stringr)
library(readr)
library(gridExtra)

setwd("C:/Users/minal/OneDrive - Texas A&M University/U-Drive/Immune_response_genes/IRG_2.0/Writing/IgHV_figure")

#### gggenome ###
chr1_positions <- read.delim("chr1_positions.tsv")

# preparing gene track
chr1_GeneTrack <- chr1_positions %>% 
  group_by(seq_id) %>%
  mutate(value = ifelse(strand == "+", min(start), max(end)),
         start = ifelse(strand == "+", (start - value), abs(start - value)),
         end = ifelse(strand == "+", (end - value), abs(end - value))) %>%
  gggenomes:::swap_if(strand == "-", start, end) %>%
  select(seq_id, start, end, strand) 

### order for the SeqTrack table
pretended_order <- c("Reference", "CS4_h1", "CS4_h2", "CS5_h1", "CS5_h2")

# preparing length track
chr1_SeqTrack <- chr1_GeneTrack %>%
  group_by(seq_id) %>%
  summarise(length=max(end)) %>%
  arrange(factor(seq_id, levels=pretended_order))

# plot
chr1_gggenome <- gggenomes(seqs=chr1_SeqTrack, genes=chr1_GeneTrack)  
chr1_gggenome +
  geom_seq(size=1) +         # draw contig/chromosome lines
  geom_bin_label(size = 4) +   # label each sequence 
  geom_gene(fill="black", color = "black", size=4, shape = 0) +        # draw genes as arrow
  theme(axis.text.x = element_text(size=12), legend.position = "none") +
  scale_x_bp(suffix = "bp", sep=" ") + 
  labs(title = "IGHV cluster", tag = "(A)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.tag = element_text())


############# Coverage ###############
# read files
filenames <- list.files(pattern = ".cov")
coverage_files <- purrr::map_df(filenames,
                                ~read.delim(.x, stringsAsFactors = FALSE, row.names = NULL, header = FALSE) %>%
                                  dplyr::mutate(filename = .x))

# tidy the dataframe
coverage_data <- coverage_files %>%
  separate(col = "filename", into = c("samples", "cov"), sep = "\\.") %>%
  select(-cov) %>%
  `colnames<-`(c("chr", "position", "coverage", "samples")) %>%
  select(samples, position, coverage) %>%
  add_column(genome_coverage = rep(c(31, 28), times = c(180451, 186471)))

## calculate rolling mean by 10000 bp
cov_rolling_10000 <- coverage_data %>%
  group_by(samples) %>%
  mutate(rolling_cov = rollapplyr(coverage, 10000, mean, partial = TRUE))

# positions excluding 50kb flanking
xmin <- min(cov_rolling_10000$position)+50000
xmax <- max(cov_rolling_10000$position)-50000

# rolling plot ##
p10000 <- ggplot(cov_rolling_10000) +
  annotate("rect", xmin=xmin, xmax=xmax, ymin=min(cov_rolling_10000$rolling_cov), ymax=max(cov_rolling_10000$rolling_cov), fill="grey", alpha=0.7) +
  geom_point(mapping=aes(x = position, y = rolling_cov), size=1) +
  geom_hline(aes(yintercept = genome_coverage), color = "red", size = 1, linetype = "longdash") +
  facet_wrap(~samples, scales = "free") +
  scale_x_continuous(labels = function(x)round(x/1000000, 2)) +
  theme_classic(base_size = 15) +
  labs(x = "Position on the reference assembly (Mb)",
       y = "Coverage (x)",
       tag = "(B)") +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.tag = element_text())



plots <- grid.arrange(chr1_gggenome, p10000, nrow = 2)

ggsave(filename = "IgHV_fig.png", plot = plots)


