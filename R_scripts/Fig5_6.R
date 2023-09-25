rm(list=ls())
library(gggenomes)
library(tidyverse)
library(msa)
library(seqinr)
library(ape)
library(pegas)
library(phytools)
library(ggtree)
library(purrr)
library(stringr)
library(readr)
library(zoo)
library(ggrepel)
library(tidyverse)
library(gridExtra)
library(ggforce)
library(Cairo)

#################################### Figure 5 ##########################################
dir <- "C:/Users/minal03/OneDrive - Texas A&M University/U-Drive/Immune_response_genes/IRG_2.0/Writing/CLM2_IFIT_analysis_final/"

################ genomic oranization ################
positions <- read.delim(paste0(dir, "gggenome/chr22_positions.txt"))

# preparing gene track
GeneTrack <- positions %>% 
  group_by(seq_id) %>%
  mutate(value = ifelse(strand == "+", min(start), max(end)),
         start = ifelse(strand == "+", (start - value), abs(start - value)),
         end = ifelse(strand == "+", (end - value), abs(end - value))) %>%
  gggenomes:::swap_if(strand == "-", start, end) %>%
  select(seq_id, start, end, gene) 

### order for the SeqTrack table
pretended_order <- c("Reference", "CS4_h1", "CS4_h2", "CS5_h1", "CS5_h2")

# preparing length track
SeqTrack <- GeneTrack %>%
  group_by(seq_id) %>%
  summarise(length=max(end)) %>%
  arrange(factor(seq_id, levels=pretended_order))

## for extracting sequences using samtools ##
contigs_positions <- positions %>% 
  group_by(seq_id) %>%
  summarise(min = min(start), max=max(end))

# plot
gggenome <- gggenomes(seqs=SeqTrack, genes=GeneTrack)  

#colorPalette <- c("#D55E00", "#E69F00", "#0072B2", "#56B4E9", "#009E73", "#CC79A7", 
#"#999999", "#000000")
colorPalette <- hcl.colors(9, palette = "Dynamic")

Fig5A <- gggenome +
  geom_seq(size=1) +         # draw contig/chromosome lines
  geom_bin_label(size = 4) +   # label each sequence 
  geom_gene(aes(fill=gene), size=3) +        # draw genes as arrow
  geom_gene_label(aes(label=gene, color = gene), fontface = "italic", size = 4, angle = 32, nudge_y = 0.20, nudge_x = -0.2) +
  #geom_text(aes(label=gene)) +
  scale_fill_manual(values= colorPalette) +
  scale_color_manual(values= colorPalette) +
  theme(axis.text.x = element_text(size=12), legend.position = "none") +
  scale_x_bp(suffix = "b", sep=" ", limits=c((0-8000), max(SeqTrack$length)+200)) +
  ggtitle('(A)')
#labs(tag = "(A)")


################## tree #############
in_file <- readDNAStringSet(paste0(dir, "CLM2_CDS.fa"))
align <- msa(in_file, "ClustalO", order = "input")
#Biostrings::writeXStringSet(DNAStringSet(align), "CLM2_aligned.fa")
class(align) <- "DNAMultipleAlignment" 
align_seqinr <- msaConvert(align, type = "seqinr::alignment")
dist <- seqinr::dist.alignment(align_seqinr, matrix = "identity")
tree <- bionj(dist)
write.tree(tree, paste0(dir, "clm2_bionj.nwk"))
clm2_tree <- read.tree(paste0(dir, "clm2_bionj_figtree.nwk"))


# metadata for colors
CLM2_labels <- data.frame(seq = rownames(align),
                            category = c("A", "B1", "B2", "C",
                                         "B2_null", "B1", "C", "A",
                                         "B2", "C", "B2_null", "B1",
                                         "C", "B1", "B2", "C",
                                         "outgroup"))

Fig5B <- ggtree(clm2_tree) %<+% CLM2_labels + 
  geom_tiplab(aes(color = category), hjust = -0.05) + 
  scale_color_manual(values=c(A="#DB9D85",B1="#86B875", B2="#38BDBB", C="#ACA4E2",
                              B2_null="#6CB4D9",
                              outgroup="#000000")) +
  xlim(0, 0.7) +
  geom_treescale() +
  theme(legend.position='none') + geom_rootedge() +
  geom_rootpoint(size=1) +
  ggtitle('(B)')

par(mfrow=c(2,1))

Fig5 <- Fig5A / Fig5B

#Fig5 <- grid.arrange(Fig5A, Fig5B, ncol = 1, nrow = 2)

#plot_grid(Fig5A, Fig5B, nrow=2, ncol=1)

ggsave(paste0(dir, "gggenome/Fig5.png"), Fig5, width = 22, height = 30, dpi = 300, units = "cm")
ggsave(paste0(dir, "gggenome/Fig5.pdf"), Fig5, device = cairo_pdf, width = 22, height = 30, dpi = 300, units = "cm")








#################################### Figure 6 ##########################################

positions <- read.delim(paste0(dir, "gggenome/chr23_positions.txt"))

# preparing gene track
GeneTrack <- positions %>% 
  group_by(seq_id) %>%
  mutate(value = ifelse(strand == "+", min(start), max(end)),
         start = ifelse(strand == "+", (start - value), abs(start - value)),
         end = ifelse(strand == "+", (end - value), abs(end - value))) %>%
  gggenomes:::swap_if(strand == "-", start, end) %>%
  select(seq_id, start, end, gene) 

### order for the SeqTrack table
pretended_order <- c("Reference", "CS4_h1", "CS4_h2", "CS5_h1", "CS5_h2")

# preparing length track
SeqTrack <- GeneTrack %>%
  group_by(seq_id) %>%
  summarise(length=max(end)) %>%
  arrange(factor(seq_id, levels=pretended_order))

## for extracting sequences using samtools ##
contigs_positions <- positions %>% 
  group_by(seq_id) %>%
  summarise(min = min(start), max=max(end))

# plot
gggenome <- gggenomes(seqs=SeqTrack, genes=GeneTrack)  

#colorPalette <- c("#D55E00", "#E69F00", "#0072B2", "#56B4E9", "#009E73", "#CC79A7", 
#"#999999", "#000000")
colorPalette <- hcl.colors(4, palette = "Dynamic")

Fig6A <- gggenome +
  geom_seq(size=1) +         # draw contig/chromosome lines
  geom_bin_label(size = 4) +   # label each sequence 
  geom_gene(aes(fill=gene), size=3) +        # draw genes as arrow
  geom_gene_label(aes(label=gene, color = gene), fontface = "italic", size = 4, angle = 32, nudge_y = 0.20, nudge_x = -0.2) +
  #geom_text(aes(label=gene)) +
  scale_fill_manual(values= colorPalette) +
  scale_color_manual(values= colorPalette) +
  theme(axis.text.x = element_text(size=12), legend.position = "none") +
  scale_x_bp(suffix = "b", sep=" ", limits=c((0-7000), max(SeqTrack$length)+8000)) +
  ggtitle('(A)')
#labs(tag = "(A)")




###### tree #########
in_file <- readDNAStringSet(paste0(dir, "IFIT10_CDS.fa"))
align <- msa(in_file, "ClustalO", order = "input")
#Biostrings::writeXStringSet(DNAStringSet(align), "CLM2_aligned.fa")
class(align) <- "DNAMultipleAlignment"
align_seqinr <- msaConvert(align, type = "seqinr::alignment")
dist <- seqinr::dist.alignment(align_seqinr, matrix = "identity")
tree <- bionj(dist)
write.tree(tree, "ifit10_bionj.nwk")
ifit10_tree <- read.tree(paste0(dir, "ifit10_bionj_figtree.nwk"))

# metadata for colors
IFIT10_labels <- data.frame(seq = rownames(align),
                            category = c("A", "B", "A", "B", rep("A", 3),
                                         "B_null", "outgroup"))

Fig6B <- ggtree(ifit10_tree) %<+% IFIT10_labels + 
  geom_tiplab(aes(color = category), hjust = -0.05) + 
  scale_color_manual(values=c(A="#DB9D85",B="#86B875", B_null="#4CB9CC", outgroup="#000000")) +
  xlim(0, 0.7) +
  geom_treescale(offset = -0.5) +
  #ggplot2::labs(tag = "(B)") +
  theme(legend.position='none') + geom_rootedge() +
  geom_rootpoint(size=1) +
  ggtitle('(B)')


Fig6 <- Fig6A / Fig6B

ggsave(paste0(dir, "gggenome/Fig6.png"), Fig6, width = 20, height = 23, dpi = 300, units = "cm")


ggsave(paste0(dir, "gggenome/Fig6.pdf"), Fig6, device = cairo_pdf, width = 20, height = 23, dpi = 300, units = "cm")

