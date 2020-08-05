#!/bin/bash

# conda activate py2

# chr8 females and males
python /home/tphung3/softwares/easySFS/easySFS.py -i /scratch/tphung3/Turkana_Demography/04_find_neutral_regions/results/chr8/chr8.gatk.called.custom.filter.autosomes.females.males.neutral.112719.vcf.gz -p pop_for_easySFS_females_males_chr8.txt -a -v --proj 36 -f -o /scratch/tphung3/Turkana_Demography/05_generate_sfs/results/chr8_projectdown_females_males

# chrX females
python /home/tphung3/softwares/easySFS/easySFS.py -i /scratch/tphung3/Turkana_Demography/04_find_neutral_regions/results/chrX/chrX.gatk.called.custom.filter.females.neutral.112719.vcf.gz -p pop_for_easySFS_females_chrX.txt -a -v --proj 36 -f -o /scratch/tphung3/Turkana_Demography/05_generate_sfs/results/chrX_projectdown_females
