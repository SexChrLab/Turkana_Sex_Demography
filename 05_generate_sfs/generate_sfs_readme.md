1. Use the script `generate_sfs.snakefile` to generate the SFS for each of the chromosome
2. Use the script `tabulate_foldedSFS.py` to combine the SFS from chr1 to chr22 for the autosomes
  ```
  python tabulate_foldedSFS.py --num_bin 82 --directory ~/scratch/Turkana_Demography/05_generate_sfs/results/ --out_filename ~/scratch/Turkana_Demography/05_generate_sfs/results/autosomes_females_males_sfs_112719.txt
  ```
3. Project down SFS
```
python /home/tphung3/softwares/easySFS/easySFS.py -i /agavescratch/tphung3/Turkana_Demography/04_find_neutral_regions/results/chrX/chrX.gatk.called.custom.filter.females.neutral.112719.vcf.gz -p pop_for_easySFS.txt -a -v --proj 36 -f -o test
```
