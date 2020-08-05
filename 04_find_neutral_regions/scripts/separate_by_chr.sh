#!/bin/bash

for i in {1..22} X Y
do

mkdir /agavescratch/tphung3/Turkana_Demography/04_find_neutral_regions/download/chr$i
done

for i in {1..22} X Y
do

grep -w chr$i /agavescratch/tphung3/Turkana_Demography/04_find_neutral_regions/download/GRCh38_gencodev32_genome_knownGene | awk '{print$1"\t"$2"\t"$3}' > /agavescratch/tphung3/Turkana_Demography/04_find_neutral_regions/download/chr$i/GRCh38_gencodev32_genome_knownGene_chr$i.bed


grep -w chr$i /agavescratch/tphung3/Turkana_Demography/04_find_neutral_regions/download/GRCh38_phastConsElements100way | awk '{print$1"\t"$2"\t"$3}' > /agavescratch/tphung3/Turkana_Demography/04_find_neutral_regions/download/chr$i/GRCh38_phastConsElements100way_chr$i.bed

grep -w chr$i /agavescratch/tphung3/Turkana_Demography/04_find_neutral_regions/download/GRCh38_repeats | awk '{print$1"\t"$2"\t"$3}' > /agavescratch/tphung3/Turkana_Demography/04_find_neutral_regions/download/chr$i/GRCh38_repeats_chr$i.bed

grep -w chr$i /agavescratch/tphung3/Turkana_Demography/04_find_neutral_regions/download/GRCh38_cpgislands | awk '{print$1"\t"$2"\t"$3}' > /agavescratch/tphung3/Turkana_Demography/04_find_neutral_regions/download/chr$i/GRCh38_cpgislands_chr$i.bed

done

for i in {1..22} X
do
grep -w chr$i /agavescratch/tphung3/Turkana_Demography/04_find_neutral_regions/download/bscores_chrfmt_high.tsv > /agavescratch/tphung3/Turkana_Demography/04_find_neutral_regions/download/chr$i/bscores_chrfmt_high_chr$i.bed
done
