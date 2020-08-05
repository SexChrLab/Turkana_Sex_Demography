#!/bin/bash
#SBATCH --job-name=neutral  # Job name
#SBATCH --mail-type=END,FAIL           # notifications for job done & fail
#SBATCH --mail-user=tphung3@asu.edu # send-to address
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 48:00:00

snakemake --snakefile select_variants_in_putatively_neutral_regions.snakefile -j 1 --cluster "sbatch --mem=40000 -c 4 -t 24:00:00"
