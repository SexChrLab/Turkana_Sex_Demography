## This file documents the demographic inference

* Working directory: `/home/tphung3/scratch/Turkana_Demography/06_demographic_inference`
1. Convert SFS to dadi file format
  ```
  grep -v af_bin ../05_generate_sfs/results/chrX/chrX.females.sfs.txt > ../05_generate_sfs/results/chrX/chrX.females.sfs.clean.txt

  python convert_sfs_to_dadi_format.py --num_bin 93 --folded_or_unfolded folded --population_name Turkana --sfs_filename ../05_generate_sfs/results/chrX/chrX.females.sfs.clean.txt --num_individuals 46 --out_filename chrX_females_sfs_02082020_dadi_format.sfs
  ```

  ```
  grep -v af_bin ../05_generate_sfs/results/chrX/chrX.males.sfs.txt > ../05_generate_sfs/results/chrX/chrX.males.sfs.clean.txt

  python convert_sfs_to_dadi_format.py --num_bin 73 --folded_or_unfolded folded --population_name Turkana --sfs_filename ../05_generate_sfs/results/chrX/chrX.males.sfs.clean.txt --num_individuals 36 --out_filename chrX_males_sfs_03052020_dadi_format.sfs
  ```

2. Run dadi
  1. chrX females
  ```
  for i in {1..50}; do python 1D.2Epoch.dadi.py --runNum ${i} --pop Turkana --mu 1.5e-8 --L 10526249 --sfs chrX_females_sfs_02082020_dadi_format.sfs --outdir 2Epoch_chrX_females_run_${i}; done;
  ```

  2. chrX males
  ```

  ```

  ```
  python merge_out.py --directory . --run_base_directory 2Epoch_chrX_males --file_basename Turkana.dadi.inference.1D.2Epoch.runNum --out_filename Turkana.dadi.inference.1D.2Epoch.chrX.males.50.runs.csv --out_filename_sorted Turkana.dadi.inference.1D.2Epoch.chrX.males.50.runs.sorted.csv
  ```

  ```
  python parse_dadi_expsfs.py --dadi_expsfs 2Epoch_run_3/Turkana.dadi.inference.1D.2Epoch.runNum.3.expSFS --num_individuals 18 --theta 9028.43395222795 --out_filename 2Epoch_chrX_males_run_3_sfs_normalized_by_theta.txt
  ```

  3. Project down to 36 chromosomes in females
  ```

  ```
