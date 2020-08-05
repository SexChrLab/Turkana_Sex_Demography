#!/bin/bash
#SBATCH --job-name=bgzip  # Job name
#SBATCH --mail-type=END,FAIL           # notifications for job done & fail
#SBATCH --mail-user=tphung3@asu.edu # send-to address
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 24:00:00

bgzip -c /scratch/tphung3/Turkana_Demography/02_genotype_calling/genotyped_vcfs/autosomes/gatk.called.raw.autosomes.females.males.vcf > /scratch/tphung3/Turkana_Demography/02_genotype_calling/genotyped_vcfs/autosomes/gatk.called.raw.autosomes.females.males.vcf.gz
