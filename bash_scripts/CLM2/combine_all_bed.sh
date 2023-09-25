#!/bin/sh
dir="/data1/CluHar/PacBio/CCS_data/Immune_response_genes/chr22_18Mb/seq_analysis/genes/"
samples="ref CS4_hap1 CS4_hap2 CS5_hap1 CS5_hap2"

### get the gene positions from all haplotypes into one file ###
for i in ${samples}; do
        awk "{\$1=\"$i\"; print;}" ${dir}$i.bed > $i.bed
done;

cat ref.bed CS*.bed > all.bed
sed -i 's/hap/h/g' all.bed

### adding 5th column as names ###
for i in {1..5}; do echo -e "CLM1\nCLM2\nCLM3\nCLM4" >> names.txt; done
paste all.bed names.txt > genes.bed
