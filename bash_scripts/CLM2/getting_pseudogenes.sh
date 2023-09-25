#!/bin/sh

assemblies="/data1/CluHar/PacBio/CCS_data/denovo/hifiasm/assemblies/"
out="/data1/CluHar/PacBio/CCS_data/Immune_response_genes/chr22_18Mb/seq_analysis/genes/getting_pseudogenes/"

### get all the sequences ###
samtools faidx /data1/CluHar/Genomes/Assembly_v2.0.2/Ch_v2.0.2.fasta chr22:18428730-18429514 > ${out}ref.fa
samtools faidx ${assemblies}CS4_hap1.fa h1tg000569l:277150-277938 > ${out}CS4_hap1.fa
samtools faidx ${assemblies}CS4_hap1.fa h1tg000569l:258167-258950 >> ${out}CS4_hap1.fa
samtools faidx ${assemblies}CS4_hap1.fa h1tg000569l:243153-243936 >> ${out}CS4_hap1.fa

samtools faidx ${assemblies}CS4_hap2.fa h2tg000391l:450539-451345 > ${out}CS4_hap2.fa

samtools faidx ${assemblies}CS5_hap1.fa h1tg000287l:310352-311129 > ${out}CS5_hap1.fa
samtools faidx ${assemblies}CS5_hap1.fa h1tg000287l:327662-328442 >> ${out}CS5_hap1.fa
samtools faidx ${assemblies}CS5_hap1.fa h1tg000287l:344274-345051 >> ${out}CS5_hap1.fa
samtools faidx ${assemblies}CS5_hap1.fa h1tg000287l:366183-366971 >> ${out}CS5_hap1.fa

samtools faidx ${assemblies}CS5_hap2.fa h2tg001147l:313166-313943 > ${out}CS5_hap2.fa
samtools faidx ${assemblies}CS5_hap2.fa h2tg001147l:287933-288716 >> ${out}CS5_hap2.fa
samtools faidx ${assemblies}CS5_hap2.fa h2tg001147l:279787-280570 >> ${out}CS5_hap2.fa
samtools faidx ${assemblies}CS5_hap2.fa h2tg001147l:267127-267915 >> ${out}CS5_hap2.fa

RT="CS4_hap1 CS5_hap2"
for rt in ${RT}; do
        seqtk seq -r ${rt}.fa > ${rt}_RT.fa
        rm ${rt}.fa
        mv ${rt}_RT.fa ${rt}.fa
done;

cat ref.fa CS*fa > pseudogenes.fa
