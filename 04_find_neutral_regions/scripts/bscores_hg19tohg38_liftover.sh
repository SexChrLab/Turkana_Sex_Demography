#!/bin/bash

for i in {1..22} X
do

~/softwares/liftOver ~/scratch/Turkana_Demography/04_find_neutral_regions/download/chr${i}/bscores_chrfmt_high_chr${i}.bed /home/tphung3/softwares/chain_files_for_liftover/hg19ToHg38.over.chain.gz ~/scratch/Turkana_Demography/04_find_neutral_regions/download/chr${i}/bscores_chrfmt_high_chr${i}_hg38liftover.bed ~/scratch/Turkana_Demography/04_find_neutral_regions/download/chr${i}/unMapped

done
