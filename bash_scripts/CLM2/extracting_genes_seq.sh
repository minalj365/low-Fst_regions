#!/bin/sh 
######### getting gene sequences ########

## get coordinates from gff file for ENSEMBL IDs
liftoff="/data1/CluHar/PacBio/CCS_data/Annotation/ENSEMBL/"
assemblies="/data1/CluHar/PacBio/CCS_data/denovo/hifiasm/assemblies/"
out="/data1/CluHar/PacBio/CCS_data/Immune_response_genes/chr22_18Mb/seq_analysis/genes/"
samples="ref CS4_hap1 CS4_hap2 CS5_hap1 CS5_hap2"
gene_IDs="ENSCHAG00000003799|ENSCHAG00000003851|ENSCHAG00000003891|ENSCHAG00000003927"

cd ${out}
### get gene positions in bed format ####
for i in ${samples}; do
        cat $liftoff${i}.gff3 | grep -E $gene_IDs | awk '{if($3=="gene") print $1,$4,$5,$7;}' >${i}.bed
done;

# for the genes on the opposite strand
cat CS4_hap1.bed | sort -r -k 2 > CS4_hap1_new.bed
rm CS4_hap1.bed
mv CS4_hap1_new.bed CS4_hap1.bed

cat CS5_hap2.bed | sort -r -k 2 > CS5_hap2_new.bed
rm CS5_hap2.bed
mv CS5_hap2_new.bed CS5_hap2.bed

#### convert gene positions into samtools format ####
for i in ${samples}; do
    cat ${i}.bed | awk 'BEGIN {OFS=""; ORS=" "}{print $1,":",$2,"-",$3;}' > ${i}.positions
    echo "samtools faidx $assemblies${i}.fa `cat ${i}.positions` > ${i}.fa" >> samtools.sh
done;

## get gene sequences
bash samtools.sh

## change header (already need to have a file with new gene names)
for i in ${samples}; do
        cat ${i}.fa | grep ">" | sed 's/>//g' > old_name_${i}.txt
        paste -d "=" old_name_${i}.txt new_name.txt > list_${i}.txt
        awk -F= 'FNR==NR {f2[$1]=$2;next} $2 in f2 {$2=f2[$2]}1' list_${i}.txt FS='>' OFS='>' ${i}.fa > ${i}_new_name.fa
        rm ${i}.fa
        mv ${i}_new_name.fa ${i}.fa
        sed -i "s/>/>${i}_/g" ${i}.fa
done;

## RT
seqtk seq -r CS4_hap1.fa > CS4_hap1_RT.fa
rm CS4_hap1.fa
mv CS4_hap1_RT.fa CS4_hap1.fa

seqtk seq -r CS5_hap2.fa > CS5_hap2_RT.fa
rm CS5_hap2.fa
mv CS5_hap2_RT.fa CS5_hap2.fa

cat ref.fa CS*.fa > genes_all.fa
#rm list* *bed *positions old* sam*

### rename fasta headers from 'hap' to 'h' ###
sed -i 's/hap/h/g' genes_all.fa

### extracting functional CDS
samtools faidx genes_all.fa
samtools faidx genes_all.fa ref_CLM2a ref_CLM2b1 ref_CLM2b2 ref_CLM2c CS4_h1_CLM2b2 CS4_h1_CLM2c CS4_h2_CLM2a CS4_h2_CLM2b2 CS4_h2_CLM2c CS5_h1_CLM2b2 CS5_h1_CLM2c CS5_h2_CLM2b1 CS5_h2_CLM2b2 CS5_h2_CLM2c > genes_functional.fa 
