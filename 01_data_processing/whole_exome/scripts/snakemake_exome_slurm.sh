#!/bin/bash
#SBATCH --job-name=kenya_exome # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=END,FAIL           # notifications for job done & fail
#SBATCH --mail-user=amtarave@asu.edu # send-to address
#SBATCH -n 1
#SBATCH -t 96:00:00
cd /home/amtarave/projects/Kenya_sequencing/exome/
source activate kenya
snakemake -j 10 --latency-wait 60 --cluster-config cluster.json --cluster "sbatch -n {cluster.n} -t {cluster.t}:00:00 --mail-type=END,FAIL --mail-user=amtarave@asu.edu "

