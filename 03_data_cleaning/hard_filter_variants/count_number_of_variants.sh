#!/bin/bash
#SBATCH --job-name=countvariants  # Job name
#SBATCH --mail-type=END,FAIL           # notifications for job done & fail
#SBATCH --mail-user=tphung3@asu.edu # send-to address
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 24:00:00

# python /home/tphung3/softwares/tanya_repos/vcfhelper/calc_num_site_in_vcf.py --vcf /agavescratch/tphung3/Turkana_Demography/02_genotype_calling/genotyped_vcfs/autosomes/gatk.called.raw.autosomes.females.males.vcf.gz > number_of_variants_raw.txt
# python /home/tphung3/softwares/tanya_repos/vcfhelper/calc_num_site_in_vcf.py --vcf /agavescratch/tphung3/Turkana_Demography/03_data_cleaning/hard_filter_variants/custom_filter_vcf/gatk.called.custom.filter.autosomes.females.males.vcf.gz > number_of_variants_hard_filter.txt

for i in {1..22}
do

python /home/tphung3/softwares/tanya_repos/vcfhelper/calc_num_site_in_vcf.py --vcf /agavescratch/tphung3/Turkana_Demography/04_find_neutral_regions/results/chr${i}/chr${i}.gatk.called.custom.filter.autosomes.females.males.neutral.112719.vcf.gz > /agavescratch/tphung3/Turkana_Demography/04_find_neutral_regions/results/chr${i}/chr${i}_num_variants_neutral.txt

done
