#!/bin/bash

# # Step 1: convert VCF files to plink file format
# plink2 --vcf /scratch/tphung3/Turkana_Demography/02_genotype_calling/genotyped_vcfs/autosomes/gatk.called.raw.autosomes.females.males.vcf.gz --make-bed --out /scratch/tphung3/Turkana_Demography/03_data_cleaning/relatedness/plink_file_format/autosomes_WGS --double-id

# Step 2: run KING
# /home/tphung3/softwares/KING_software/king -b /scratch/tphung3/Turkana_Demography/03_data_cleaning/relatedness/plink_file_format/autosomes_WGS.bed --unrelated --degree 3
# /home/tphung3/softwares/KING_software/king -b /scratch/tphung3/Turkana_Demography/03_data_cleaning/relatedness/plink_file_format/chr21_WGS.bed --related

/home/tphung3/softwares/KING_software/king -b /scratch/tphung3/Turkana_Demography/03_data_cleaning/relatedness/plink_file_format/autosomes_WGS.bed --kinship
