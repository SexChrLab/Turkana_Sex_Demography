#!/bin/bash

# Date: 11/27/2019
# Obtain putatively neutral regions by removing genic regions, conserved regions, repeats, CpG islands and keep sites where b values is greater or equal to 0.9

for i in {1..22} X
do

bedtools subtract -a /agavescratch/tphung3/Turkana_Demography/04_find_neutral_regions/download/chr${i}/chr${i}.bed -b /agavescratch/tphung3/Turkana_Demography/04_find_neutral_regions/download/chr${i}/GRCh38_gencodev32_genome_knownGene_chr${i}.bed | bedtools subtract -a stdin -b /agavescratch/tphung3/Turkana_Demography/04_find_neutral_regions/download/chr${i}/GRCh38_phastConsElements100way_chr${i}.bed | bedtools subtract -a stdin -b  /agavescratch/tphung3/Turkana_Demography/04_find_neutral_regions/download/chr${i}/GRCh38_repeats_chr${i}.bed | bedtools subtract -a stdin -b /agavescratch/tphung3/Turkana_Demography/04_find_neutral_regions/download/chr${i}/GRCh38_cpgislands_chr${i}.bed | bedtools intersect -a stdin -b /agavescratch/tphung3/Turkana_Demography/04_find_neutral_regions/download/chr${i}/bscores_chrfmt_high_chr${i}_hg38liftover.bed > /agavescratch/tphung3/Turkana_Demography/04_find_neutral_regions/results/chr${i}/chr${i}_putatively_neutral_regions_112719.bed

done
