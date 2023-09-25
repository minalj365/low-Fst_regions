#!/bin/sh
######## getting gene sequences from annotation files #####

gff="/data1/CluHar/PacBio/CCS_data/Immune_response_genes/IGHV/HiCanu_Hifiasm/Annotations/making_HiCanu_gtf_file/gff/hifiasm_gff/"
out="/data1/CluHar/PacBio/CCS_data/Immune_response_genes/IGHV/seq_analysis/alignments/genes/"
hifiasm_haps="/data1/CluHar/PacBio/CCS_data/Immune_response_genes/IGHV/HiCanu_Hifiasm/hifiasm_haplotypes/"
samples="CS4_hap1 CS4_hap2 CS5_hap1 CS5_hap2"



for i in ${samples}; do
        cat $gff${i}.gff | grep -v "#" | awk '{print $1,$4,$5,$7;}' >${i}.bed
done;

#### convert gene positions into samtools format ####
for i in ${samples}; do
    cat ${i}.bed | awk 'BEGIN {OFS=""; ORS=" "}{print $1,":",$2,"-",$3;}' > ${i}.positions
    echo "samtools faidx $hifiasm_haps${i}.fa `cat ${i}.positions` > ${i}.fa" >> samtools.sh
done;

## get gene sequences
bash samtools.sh
