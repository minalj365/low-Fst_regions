#!/bin/sh 
######### getting CDS and protein sequences ##########
assemblies="/data1/CluHar/PacBio/CCS_data/denovo/hifiasm/assemblies/"
liftoff="/data1/CluHar/PacBio/CCS_data/Annotation/ENSEMBL/"
out="/data1/CluHar/PacBio/CCS_data/Immune_response_genes/chr22_18Mb/seq_analysis/CDS_and_AA/"

samples="ref CS4_hap1 CS4_hap2 CS5_hap1 CS5_hap2"

# second gene has two transcripts - ENSCHAT00000006955 and ENSCHAT00000006972. Only taking first one.
transcript_IDs="ENSCHAT00000006887 ENSCHAT00000006955  ENSCHAT00000007031 ENSCHAT00000007093"
new_names="CLM2a CLM2b1 CLM2b2 CLM2c"

for i in ${samples}; do
	#### get CDS sequences corresponding to transcript_IDs
	samtools faidx ${liftoff}CDS/${i}.fa ${transcript_IDs} > ${out}${i}_CDS.fa

	#### get protein sequences corresponding to transcript_IDs
	samtools faidx ${liftoff}AA/${i}.fa ${transcript_IDs} > ${out}${i}_AA.fa
done;

## change header (already need to have a file with new gene names)
################## how to write file in a bash script?
cd ${out}
for i in ${samples}; do
	### for CDS
    cat ${i}_CDS.fa | grep ">" | sed 's/>//g' > old_name_${i}.txt
    paste -d "=" old_name_${i}.txt new_name.txt > list_${i}.txt
    awk -F= 'FNR==NR {f2[$1]=$2;next} $2 in f2 {$2=f2[$2]}1' list_${i}.txt FS='>' OFS='>' ${i}_CDS.fa > ${i}_new_name.fa
    rm ${i}_CDS.fa
    mv ${i}_new_name.fa ${i}_CDS.fa
    sed -i "s/>/>${i}_/g" ${i}_CDS.fa

    ### for AA
    cat ${i}_AA.fa | grep ">" | sed 's/>//g' > old_name_${i}.txt
    paste -d "=" old_name_${i}.txt new_name.txt > list_${i}.txt
    awk -F= 'FNR==NR {f2[$1]=$2;next} $2 in f2 {$2=f2[$2]}1' list_${i}.txt FS='>' OFS='>' ${i}_AA.fa > ${i}_new_name.fa
    rm ${i}_AA.fa
    mv ${i}_new_name.fa ${i}_AA.fa
    sed -i "s/>/>${i}_/g" ${i}_AA.fa
done;
cat ref_CDS.fa CS*_CDS.fa > CDS_all.fa
cat ref_AA.fa CS*_AA.fa > AA_all.fa
#rm list* *bed *positions old* sam*

### rename fasta headers from 'hap' to 'h' ###
sed -i 's/hap/h/g' CDS_all.fa
sed -i 's/hap/h/g' AA_all.fa

##### LliftOver annotation names are not in order for CS4_hap1 and CS5_hap1. Hence renaming those sequences manually #####
### The changes are 
# CS4_hap1_CLM2b to CS4_hap1_CLM2c
# CS4_hap1_CLM2c to CS4_hap1_CLM2b
# CS5_hap1_CLM2b to CS5_hap1_CLM2c
# CS5_hap1_CLM2b to CS5_hap1_CLM2c

samtools faidx CDS_all.fa

### extracting functional CDS
samtools faidx CDS_all.fa ref_CLM2a ref_CLM2b1 ref_CLM2b2 ref_CLM2c CS4_h1_CLM2b2 CS4_h1_CLM2c CS4_h2_CLM2a CS4_h2_CLM2b2 CS4_h2_CLM2c CS5_h1_CLM2b2 CS5_h1_CLM2c CS5_h2_CLM2b1 CS5_h2_CLM2b2 CS5_h2_CLM2c > CDS_functional.fa 
