############ Multi sequence alignment, nucleotide diversity, phylogenetic tree #########
rm(list=ls())
setwd("C:/Users/minal03/OneDrive - Texas A&M University/U-Drive/Immune_response_genes/IRG_2.0/Writing/Seq_analysis")
library(msa)
library(seqinr)

files=c("CLM2a", "CLM2b1", "CLM2b2", "CLM2c", "IFIT10a", "IFIT10b")

## CDS pi calculations ##

CDS_pi_table = data.frame()   # creating an empty table to write pi values of CDS sequences

for (i in files) {
  in_file <- readDNAStringSet(paste0("CDS/", i, ".fa"))   # reading files from CDS folder
  align <- msa(in_file, "ClustalOmega", order = "input")
  CDS_pi <- round(nuc.div(as.DNAbin(align)), digits = 3)
  CDS_pi_table <- rbind(CDS_pi_table, CDS_pi)
}

rownames(CDS_pi_table) <- files
colnames(CDS_pi_table) <- c("pi")
CDS_pi_table
write.table(CDS_pi_table, "CDS/R_outputs/CDS_pi.tsv", sep="\t")

