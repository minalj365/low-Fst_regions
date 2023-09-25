#!/bin/sh 

### renaming fasta header ####
samples="CS4_hap1 CS4_hap2 CS5_hap1 CS5_hap2"
for i in ${samples}; do
 sed -i "s/>.*/>${i}/" ${i}.fa
done;
