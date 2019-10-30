# Angela Taravella
# 2017-11-07
# Notes for Kenya data prep using snakemake

##############
# DIRECTORIES
##############
# Data on cluster
/mnt/storage/CNTRLDRAW/WES-11-29-2017-YALE/ycga-ba/ba_sequencers2/scratch/ccc7/UserDataLinks/Sayres/WES-11-13-2017/
/mnt/storage/CNTRLDRAW/WGS-11-29-2017-YALE/ycga-ba/ba_sequencers2/scratch/ccc7/UserDataLinks/Sayres/WGS-11-9-2017/
# new location
/mnt/storage/CNTRLDPROCESSED/WES-11-29-2017-YALE/ycga-ba/ba_sequencers2/scratch/ccc7/UserDataLinks/Sayres/WES-11-13-2017/
/mnt/storage/CNTRLDPROCESSED/WGS-11-29-2017-YALE/ba_sequencers2/scratch/ccc7/UserDataLinks/Sayres/WGS-11-9-2017/
#### THIS IS OLD #########################################################################
########
# Step 1
########
# I will make a separate environment for each step in data prep.
# Evironment 1: QC
# Set up conda environment for QC (quality control)
cd /home/amtarave/Kenya_sequencing/QC/
conda env create --name kenya_pipeline_QC --file kenya_pipeline_QC.yaml

#
# To activate this environment, use:
# > source activate kenya_pipeline_QC
#
# To deactivate this environment, use:
# > source deactivate kenya_pipeline_QC
#
# I am going to change the environment name to kenya
conda create --name kenya --clone kenya_pipeline_QC
conda remove --name kenya_pipeline_QC --all
#
#
# To activate this environment, use:
# > source activate kenya
#
# To deactivate this environment, use:
# > source deactivate kenya
#

##########################################################################################

# 2018-02-28
# I will create 2 snakefiles, one for the whole genome, one for the exome but all of the
#	steps will be in the snakefile (as opposed to breaking it up by step)

########
# Step 1
########
# Snakemake #

# Step 1.01: make json file for snakemake
# For each unique pair of fastq files, I want to add paticular information into the json
#	file so that we can call it later in the snakefile.
# This is what I want it to look like in the json file:
	"A10_H5G7VDMXX_L002":{
		"fq1": "A10_BH5G7VDMXX_L002_R1_001.fastq.gz",
		"fq2": "A10_BH5G7VDMXX_L002_R2_001.fastq.gz",
		"ID": "A10_H5G7VDMXX_L002",
		"SM": "A10_exome",
		"LB": "A10_exome",
		"PU": "H5G7VDMXX.2",
		"PL": "Illumina"
	},
# I will try to use scripting to do this, to avoid human (copy and pasting) error.
for i in *R1*.gz; do echo $i; zcat $i | head -1; done

# Make an bash script to get a list of all R1/2 pairs of files for all samples. Example for A10:
A10_AH575TDMXX_L001_R1_001.fastq.gz
@A00127:24:H575TDMXX:1:1101:1253:1000 1:N:0:ATTACTCG+CAGGACGT
A10_AH575TDMXX_L002_R1_001.fastq.gz
@A00127:24:H575TDMXX:2:1101:1253:1000 1:N:0:ATTACTCG+CAGGACGT
A10_BH5G7VDMXX_L001_R1_001.fastq.gz
@A00127:26:H5G7VDMXX:1:1101:2537:1063 1:N:0:ATTACTCG+NAGGACGT
A10_BH5G7VDMXX_L002_R1_001.fastq.gz
@A00127:26:H5G7VDMXX:2:1101:1416:1063 1:N:0:ATTACTCG+NAGGACGT

######### SIDE NOTE ########################
#***NOTE: info on header info from fastq file
https://en.wikipedia.org/wiki/FASTQ_format#Illumina_sequence_identifiers
# File name example
A10_AH575TDMXX_L001_R1_001.fastq.gz

# Header
@A00127:24:H575TDMXX:1:1101:1253:1000 2:N:0:ATTACTCG+CAGGACGT

@A00127: machine
24: Run ID
H575TDMXX: Flow cell
1: flow cell lane
1101: tile number within flow cell lane
1253: 'x'-coordinate of the cluster within the tile
1000: 'y'-coordinate of the cluster within the tile
1: the member of a pair, 1 or 2 (paired-end or mate-pair reads only)
N: Y if the read is filtered, N otherwise
0: 0 when none of the control bits are on, otherwise it is an even numbe
ATTACTCG+CAGGACGT: index sequence

# Read group info
https://software.broadinstitute.org/gatk/documentation/article?id=6472
#############################################


## Then I will make a python script that parses this file and generates proper info for
#	json file.
[amtarave@cg1-2:~/projects/Kenya_sequencing/scripts]$ bash get_files_sample_info.sh > exome_info.txt
[amtarave@cg1-2:~/projects/Kenya_sequencing/scripts]$ bash get_files_sample_info.sh > whole_genome_info.txt

# I want the output for one file set to be on one line not two so like:
A10_AH575TDMXX_L001_R1_001.fastq.gz 	@A00127:24:H575TDMXX:1:1101:1253:1000 1:N:0:ATTACTCG+CAGGACGT
# So I made two separate lists (one with fastq file name and one with the fastq info) and I will paste them

[amtarave@cg1-2:~/projects/Kenya_sequencing]$ paste exome_files.txt exome_list.txt > exome_info_one_line.txt
[amtarave@cg1-2:~/projects/Kenya_sequencing]$ head exome_info_one_line.txt
A10_BHMG7GBBXX_L005_R1_001.fastq.gz     @K00175:103:HMG7GBBXX:5:1101:1367:1156 1:N:0:NATCAG+NCAAGT
A100_BHMCGJBBXX_L008_R1_001.fastq.gz    @K00162:217:HMCGJBBXX:8:1101:1083:1156 1:N:0:NGATGT+NTAGCC
A11_BHMG7GBBXX_L005_R1_001.fastq.gz     @K00175:103:HMG7GBBXX:5:1101:1326:1156 1:N:0:NAGCTT+NTGATC

[amtarave@cg1-2:~/projects/Kenya_sequencing]$ paste whole_genome_files.txt whole_genome_list.txt > whole_genome_info_one_line.txt
[amtarave@cg1-2:~/projects/Kenya_sequencing]$ head whole_genome_info_one_line.txt
A10_AH575TDMXX_L001_R1_001.fastq.gz     @K00175:103:HMG7GBBXX:5:1101:1367:1156 1:N:0:NATCAG+NCAAGT
A10_AH575TDMXX_L002_R1_001.fastq.gz     @K00162:217:HMCGJBBXX:8:1101:1083:1156 1:N:0:NGATGT+NTAGCC
A10_BH5G7VDMXX_L001_R1_001.fastq.gz     @K00175:103:HMG7GBBXX:5:1101:1326:1156 1:N:0:NAGCTT+NTGATC

## Run python script (once for exomes once for whole genome) to get a text file with info
#	for each pair of fq files for all individuals. File names are hard coded.
python ~/projects/Kenya_sequencing/scripts/get_json_info.py
python ~/projects/Kenya_sequencing/scripts/get_json_info.py

# Then just copy and paste results into ~/projects/Kenya_sequencing/whole_genome/kenya_config_whole_genome.json
#	and ~/projects/Kenya_sequencing/exome/kenya_config_exome.json


## I also need to make a "unique_identifiers" section in each json file this will have every
#	unique ID for all samples. For example: A10_H575TDMXX_L001
# I edited ~/projects/Kenya_sequencing/scripts/get_json_info.py to do this.



kenya_pipeline_QC) [amtarave@cg2-11:~/projects/Kenya_sequencing/exome/fastqc]$ ll
total 532
-rw-rw-r-- 1 amtarave amtarave 233053 Mar  2 14:00 A10_BHMG7GBBXX_L005_R1_001_fastqc.html
-rw-rw-r-- 1 amtarave amtarave 308582 Mar  2 14:00 A10_BHMG7GBBXX_L005_R1_001_fastqc.zip



## Make symbolic links for fastq files
# See exome/Snakefile


# Step 1.02: Run exome snakemake
# I figured out the syntax for the exome snakefile and I ran a successful dry run.
# Currently, the snakefile makes symbolic links for all the exmome fastq files, runs
# fastqc and multiqc.

# I have to edit the cluster json file so that I can alocate the accurate options per
# rule (for example the correct number of cpus).

[amtarave@agave1:~/projects/Kenya_sequencing/exome]$ sbatch snakemake_exome_slurm.sh


 # Step 1.03: Run whole genome snakemake
 # Make sure to format the json file like the exome json file. I did some things
 #	by hand, but this time I will edit the get_json_info.py to print out the
 # correct information.

 # run get_json_info.py
python get_json_info.py

# dry run Snakefile for whole genome
(kenya_pipeline_QC) [amtarave@cg17-4:~/projects/Kenya_sequencing/whole_genome]$ snakemake -np
Job counts:
        count   jobs
        1       all
        436     fastqc_analysis
        218     make_symbolic_link_for_fastqs
        1       multiqc_analysis
        656

########
# Step 2
########
# Sankemake: trimming and re-run fastqc and multiqc
##########################################################################################
# OLD #
# I will use trimmomatic for trimmming but i am unsure which .fa file to use when trimming
#	I found online that you can download bbmap and use their adapters.fa file which
#	should have all the illumina adapters.
# https://www.biostars.org/p/227726/

# Directory to adapter file.
/home/amtarave/packages/bbmap/resources/adapters.fa

# example trimming command in paired end mode (PE), sliding window with quality score 30 and above, and trimming adapters
trimmomatic PE SRR519926_1.fastq SRR519926_2.fastq trimmed_1.fq unpaired_1.fq trimmed_2.fq unpaired_2.fq SLIDINGWINDOW:4:30 ILLUMINACLIP:/home/amtarave/packages/bbmap/resources/adapter.fa:2:30:5

# Step 2.01: Install trimmomatic into kenya environment
source activate kenya
conda install -c bioconda trimmomatic
##########################################################################################
# I will actually use bbduk to trim the data and then re-run fastqc and multiqc on the
#	trimmed data.

# Step 2.01: Install bbduk into kenya environment. It looks like bbmap has bbduk sh file...
conda install -c bioconda bbmap

# Step 2.02: incorporate trimmomatic command into snakefile.
# Tim provided some code to do this. See below:
rule trim_adapters_paired_bbduk:
	input:
		fq1 = lambda wildcards: os.path.join(
			fastq_directory, config[wildcards.sample]["fq1"]),
		fq2 = lambda wildcards: os.path.join(
			fastq_directory, config[wildcards.sample]["fq2"])
	output:
		out_fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		out_fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
	params:
		bbduksh = bbduksh_path
	shell:
		"{params.bbduksh} -Xmx3g in1={input.fq1} in2={input.fq2} "
		"out1={output.out_fq1} out2={output.out_fq2} "
		"ref=misc/adapter_sequence.fa ktrim=r k=21 mink=11 hdist=2 tbo tpe "
		"qtrim=rl trimq=15 minlen=50 maq=20"

# I should probably change minlen to 75 because my read are longer.
# hdist is conservative can change to 1.

# I think I need to move this to my scratch directory because the fastq files might be too
# 	big for my home directory. in the trimming rule, I have the fastq files outputting to
# 	my scratch directory.

# Step 2.03: Run snakemake dry runs
# WGS
Job counts:
        count   jobs
        1       all
        218     fastqc_analysis_trimmed
        1       multiqc_analysis_trimmed
        218     trim_adapters_paired_bbduk
        438

########
# Step 3
########
# Bam generation step #

# Sankemake: Mapping, adding readgroups, fixmates, sort, indexing, merging bam files
# small issue: I ran trimming on ocotillo (in my scratch directory) for wes and that data
#	is not on agave. So do I re-run it on agave or do I here after run WES on ocotillo?

# WGS
(kenya) [amtarave@cg8-9:~/projects/Kenya_sequencing/whole_genome]$ snakemake -np
Job counts:
        count   jobs
        1       all
        218     index_bam
        218     map_and_process_trimmed_reads
        437

# I will be mapping to sex specific reference files.
# Here I will make the index files using bwa index
cd REF/
mkdir grch38
cd grch38
pwd
/home/amtarave/REF/grch38
ln -s /mnt/storage/SAYRES/XY_Trim_Ref/references/gencode.GRCh38.p7_wholeGenome/GRCh38_wholeGenome_reference.fa GRCh38_wholeGenome_ref.fa
samtools faidx GRCh38_wholeGenome_ref.fa
samtools dict -o GRCh38_wholeGenome_ref.dict GRCh38_wholeGenome_ref.fa
bwa index GRCh38_wholeGenome_ref.fa

### There also seems to be an issue with samtools in general

# this will have to change in config file
"Ref_GRCh38": "/mnt/storage/SAYRES/XY_Trim_Ref/references/gencode.GRCh38.p7_wholeGenome/GRCh38_wholeGenome_reference.fa",
"Ref_GRCh38_Y_HardMasked": "/mnt/storage/SAYRES/XY_Trim_Ref/references/gencode.GRCh38.p7_Ymasked/GRCh38_Ymasked_reference.fa",
"Ref_GRCh38_Y_PARsMasked": "/mnt/storage/SAYRES/XY_Trim_Ref/references/gencode.GRCh38.p7_minusYPARs/GRCh38_minusYPARs_reference.fa",

bwa index <ref.fa>

# practice with this command until it works. I will use the regular ref for testing
bwa mem -t 4 -R '@RG\tID:A74_H3HFLDMXX_L001\tSM:A74_whole_genome\tLB:A74_whole_genome\tPU:H3HFLDMXX.1\tPL:Illumina' /home/amtarave/REF/grch38/GRCh38_minusYPARs_ref.fa /scratch/amtarave/Kenya_agave/whole_genome/trimmed_fastqs/A74_H3HFLDMXX_L001_trimmed_R1.fastq.gz /scratch/amtarave/Kenya_agave/whole_genome/trimmed_fastqs/A74_H3HFLDMXX_L001_trimmed_R2.fastq.gz | samtools fixmate -O bam - - | samtools sort -O bam -o /scratch/amtarave/Kenya_agave/whole_genome/processed_bams/A74_H3HFLDMXX_L001.GRCh38_TEST.sorted.bam


#NOTE: There is no merging bams for the exome since each sample will only have
#	one bam file


#### Testing XYalign ####
python /home/amtarave/packages/XYalign-1.0.0/xyalign_runner.py --ANALYZE_BAM --ref /home/amtarave/projects/Kenya_sequencing/exome/refs/Ref_GRCh38_Y_HardMasked.fa --bam /scratch/amtarave/Kenya/processed_bams/B3_HMCGJBBXX_L005.GRCh38_Ymasked.sorted.bam --output_dir /scratch/amtarave/Kenya/test_xyalign/ --sample_id B3_HMCGJBBXX_L005 --cpus 4 --target_bed  /home/amtarave/projects/Kenya_sequencing/exome/xGen_capture_kit_bed_files/xgen-exome-research-panel-targets.bed --chromosomes chr19 chrX  --x_chromosome chrX
# To activate this environment, use:
# > source activate xyalign_env
#
# To deactivate this environment, use:
# > source deactivate xyalign_env
#

# EXOME #
# Test with including chr Y
# A2 and B3 are females
python /home/amtarave/packages/XYalign-1.0.0/xyalign_runner.py --ANALYZE_BAM --ref /home/amtarave/projects/Kenya_sequencing/exome/refs/Ref_GRCh38_Y_HardMasked.fa --bam /scratch/amtarave/Kenya/processed_bams/B3_HMCGJBBXX_L005.GRCh38_Ymasked.sorted.bam --output_dir /scratch/amtarave/Kenya/test_xyalign/B3_withchrY/ --sample_id B3_HMCGJBBXX_L005 --cpus 4 --target_bed  /home/amtarave/projects/Kenya_sequencing/exome/xGen_capture_kit_bed_files/xgen-exome-research-panel-targets.bed --chromosomes chr19 chrX chrY --x_chromosome chrX --y_chromosome chrY
python /home/amtarave/packages/XYalign-1.0.0/xyalign_runner.py --CHARACTERIZE_SEX_CHROMS --ref /home/amtarave/projects/Kenya_sequencing/exome/refs/Ref_GRCh38_Y_HardMasked.fa --bam /scratch/amtarave/Kenya/processed_bams/B3_HMCGJBBXX_L005.GRCh38_Ymasked.sorted.bam --output_dir /scratch/amtarave/Kenya/test_xyalign/B3_withchrY/ --sample_id B3_HMCGJBBXX_L005 --cpus 4 --target_bed  /home/amtarave/projects/Kenya_sequencing/exome/xGen_capture_kit_bed_files/xgen-exome-research-panel-targets.bed --chromosomes chr19 chrX chrY --x_chromosome chrX --y_chromosome chrY

python /home/amtarave/packages/XYalign-1.0.0/xyalign_runner.py --CHARACTERIZE_SEX_CHROMS --ref /home/amtarave/projects/Kenya_sequencing/exome/refs/Ref_GRCh38_Y_PARsMasked.fa --bam /scratch/amtarave/Kenya/processed_bams/A2_HMG7GBBXX_L004.GRCh38_Ymasked.sorted.bam --output_dir /scratch/amtarave/Kenya/test_xyalign/A2_withchrY/ --sample_id A2_HMG7GBBXX_L004 --cpus 4 --target_bed  /home/amtarave/projects/Kenya_sequencing/exome/xGen_capture_kit_bed_files/xgen-exome-research-panel-targets.bed --chromosomes chr19 chrX chrY --x_chromosome chrX --y_chromosome chrY

python /home/amtarave/packages/XYalign-1.0.0/xyalign_runner.py --CHROM_STATS --use_counts --bam /scratch/amtarave/Kenya/processed_bams/B3_HMCGJBBXX_L005.GRCh38_Ymasked.sorted.bam --ref null --output_dir  /scratch/amtarave/Kenya/test_xyalign/B3_withchrY/ --sample_id B3_HMCGJBBXX_L005 --chromosomes chr19 chrX chrY
python /home/amtarave/packages/XYalign-1.0.0/xyalign_runner.py --CHROM_STATS --use_counts --bam /scratch/amtarave/Kenya/processed_bams/A2_HMG7GBBXX_L004.GRCh38_Ymasked.sorted.bam --ref null --output_dir  /scratch/amtarave/Kenya/test_xyalign/A2_withchrY/ --sample_id A2_HMG7GBBXX_L004 --chromosomes chr19 chrX chrY

# A10 is a male
python /home/amtarave/packages/XYalign-1.0.0/xyalign_runner.py --CHROM_STATS --use_counts --bam /scratch/amtarave/Kenya/processed_bams/A10_HMG7GBBXX_L005.GRCh38_minusYPARs.sorted.bam --ref null --output_dir  /scratch/amtarave/Kenya/test_xyalign/A10_withchrY/ --sample_id A10_HMG7GBBXX_L005 --chromosomes chr19 chrX chrY
python /home/amtarave/packages/XYalign-1.0.0/xyalign_runner.py --CHARACTERIZE_SEX_CHROMS --ref /home/amtarave/projects/Kenya_sequencing/exome/refs/Ref_GRCh38_Y_PARsMasked.fa --bam /scratch/amtarave/Kenya/processed_bams/A10_HMG7GBBXX_L005.GRCh38_minusYPARs.sorted.bam --output_dir /scratch/amtarave/Kenya/test_xyalign/A10_withchrY/ --sample_id A10_HMG7GBBXX_L005 --cpus 4 --target_bed  /home/amtarave/projects/Kenya_sequencing/exome/xGen_capture_kit_bed_files/xgen-exome-research-panel-targets.bed --chromosomes chr19 chrX chrY --x_chromosome chrX --y_chromosome chrY

# using small window instead of target bed
python /home/amtarave/packages/XYalign-1.0.0/xyalign_runner.py --CHARACTERIZE_SEX_CHROMS --ref /home/amtarave/projects/Kenya_sequencing/exome/refs/Ref_GRCh38_Y_PARsMasked.fa --bam /scratch/amtarave/Kenya/processed_bams/A2_HMG7GBBXX_L004.GRCh38_Ymasked.sorted.bam --output_dir /scratch/amtarave/Kenya/test_xyalign/A2_withchrY_5000window/ --sample_id A2_HMG7GBBXX_L004 --cpus 4 --window_size 5000 --chromosomes chr19 chrX chrY --x_chromosome chrX --y_chromosome chrY
python /home/amtarave/packages/XYalign-1.0.0/xyalign_runner.py --CHARACTERIZE_SEX_CHROMS --ref /home/amtarave/projects/Kenya_sequencing/exome/refs/Ref_GRCh38_Y_PARsMasked.fa --bam /scratch/amtarave/Kenya/processed_bams/A10_HMG7GBBXX_L005.GRCh38_minusYPARs.sorted.bam --output_dir /scratch/amtarave/Kenya/test_xyalign/A10_withchrY_5000window/ --sample_id A10_HMG7GBBXX_L005 --cpus 4 --window_size 5000 --chromosomes chr19 chrX chrY --x_chromosome chrX --y_chromosome chrY


# WHOLE GENOME #
python /home/amtarave/packages/XYalign-1.0.0/xyalign_runner.py --CHARACTERIZE_SEX_CHROMS --ref /home/amtarave/projects/Kenya_sequencing/exome/refs/Ref_GRCh38_Y_PARsMasked.fa --bam /scratch/amtarave/Kenya_agave/whole_genome/processed_bams/A2.GRCh38_Ymasked.sorted.merged.bam --output_dir /scratch/amtarave/Kenya_agave/whole_genome/test_xyalign/A2_withchrY/ --sample_id A2 --cpus 4 --window_size 10000 --chromosomes chr19 chrX chrY --x_chromosome chrX --y_chromosome chrY
python /home/amtarave/packages/XYalign-1.0.0/xyalign_runner.py --CHARACTERIZE_SEX_CHROMS --ref /home/amtarave/projects/Kenya_sequencing/exome/refs/Ref_GRCh38_Y_PARsMasked.fa --bam /scratch/amtarave/Kenya_agave/whole_genome/processed_bams/A10.GRCh38_minusYPARs.sorted.merged.bam --output_dir /scratch/amtarave/Kenya_agave/whole_genome/test_xyalign/A10_withchrY/ --sample_id A10 --cpus 4 --window_size 10000 --chromosomes chr19 chrX chrY --x_chromosome chrX --y_chromosome chrY

########
# Step 4
########
# Fix snakemake error
# I received an OUT_OF_MEMORY error when running snakemake. Specifically the bwa mapping step.

# The following samples received this error:
A34_H3HFLDMXX_L002.GRCh38_minusYPARs.sorted.bam
A23_H57HCDMXX_L001.GRCh38_minusYPARs.sorted.bam
A67_H3HCVDMXX_L002.GRCh38_minusYPARs.sorted.bam
A63_H3HCNDMXX_L001.GRCh38_minusYPARs.sorted.bam
A69_H5G7VDMXX_L001.GRCh38_minusYPARs.sorted.bam
A69_H5G7VDMXX_L002.GRCh38_minusYPARs.sorted.bam
A24_H575TDMXX_L001.GRCh38_minusYPARs.sorted.bam
A32_H3GWKDMXX_L001.GRCh38_minusYPARs.sorted.bam
A67_H3HCVDMXX_L001.GRCh38_minusYPARs.sorted.bam
A63_H3HCNDMXX_L002.GRCh38_minusYPARs.sorted.bam
A23_H57HCDMXX_L002.GRCh38_minusYPARs.sorted.bam
A32_H3GWKDMXX_L002.GRCh38_minusYPARs.sorted.bam
A40_H3HCNDMXX_L002.GRCh38_minusYPARs.sorted.bam
A40_H3HCNDMXX_L001.GRCh38_minusYPARs.sorted.bam
A21_H3GWKDMXX_L001.GRCh38_minusYPARs.sorted.bam
A72_H5G7VDMXX_L002.GRCh38_minusYPARs.sorted.bam
A21_H3GWKDMXX_L002.GRCh38_minusYPARs.sorted.bam
A60_H5G7VDMXX_L002.GRCh38_minusYPARs.sorted.bam
A18_H3HFLDMXX_L001.GRCh38_minusYPARs.sorted.bam
A95_H57HCDMXX_L001.GRCh38_minusYPARs.sorted.bam
A48_H3HFLDMXX_L001.GRCh38_minusYPARs.sorted.bam
A18_H3HFLDMXX_L002.GRCh38_minusYPARs.sorted.bam
A93_H3HCVDMXX_L001.GRCh38_minusYPARs.sorted.bam
A48_H3HFLDMXX_L002.GRCh38_minusYPARs.sorted.bam
A52_H3HCNDMXX_L001.GRCh38_minusYPARs.sorted.bam
A93_H3HCVDMXX_L002.GRCh38_minusYPARs.sorted.bam
A22_H3HFLDMXX_L002.GRCh38_minusYPARs.sorted.bam
A52_H3HCNDMXX_L002.GRCh38_minusYPARs.sorted.bam
A81_H3HCVDMXX_L002.GRCh38_minusYPARs.sorted.bam
A36_H5G7VDMXX_L001.GRCh38_minusYPARs.sorted.bam
A98_H57HCDMXX_L002.GRCh38_minusYPARs.sorted.bam

# Remove these files that were generate yet received the out of memory error:
rm A34_H3HFLDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A23_H57HCDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A67_H3HCVDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A63_H3HCNDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A69_H5G7VDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A69_H5G7VDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A24_H575TDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A32_H3GWKDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A67_H3HCVDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A63_H3HCNDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A23_H57HCDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A32_H3GWKDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A40_H3HCNDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A40_H3HCNDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A21_H3GWKDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A72_H5G7VDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A21_H3GWKDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A60_H5G7VDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A18_H3HFLDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A95_H57HCDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A48_H3HFLDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A18_H3HFLDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A93_H3HCVDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A48_H3HFLDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A52_H3HCNDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A93_H3HCVDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A22_H3HFLDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A52_H3HCNDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A81_H3HCVDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A36_H5G7VDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A98_H57HCDMXX_L002.GRCh38_minusYPARs.sorted.bam


## Frist try request more cores (n) in cluster config
# If this doesn't work, try adding a memory option.
# change n = 4 to n = 8 for mapping rules
[amtarave@agave1:~/projects/Kenya_sequencing/whole_genome]$ sbatch snakemake_whole_genome_slurm.sh
Submitted batch job 145559

# This did not work so I will try requesting more memory per cpu (--mem-per-cpu=16000)


## Try requesting more memory per cpu
# Make a list of samples that failed because they ran out of memory. I want to delete these
#	files because they can be truncated.
rm A24_H575TDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A95_H57HCDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A93_H3HCVDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A93_H3HCVDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A18_H3HFLDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A18_H3HFLDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A98_H57HCDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A48_H3HFLDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A23_H57HCDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A52_H3HCNDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A23_H57HCDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A48_H3HFLDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A52_H3HCNDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A22_H3HFLDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A32_H3GWKDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A81_H3HCVDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A21_H3GWKDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A36_H5G7VDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A40_H3HCNDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A40_H3HCNDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A63_H3HCNDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A21_H3GWKDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A67_H3HCVDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A60_H5G7VDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A32_H3GWKDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A63_H3HCNDMXX_L001.GRCh38_minusYPARs.sorted.bam
rm A67_H3HCVDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A72_H5G7VDMXX_L002.GRCh38_minusYPARs.sorted.bam
rm A34_H3HFLDMXX_L002.GRCh38_minusYPARs.sorted.bam

# Added mem-per-cpu flag to the sbatch and config files
# Do a dry run
(kenya) [amtarave@cg17-7:~/projects/Kenya_sequencing/whole_genome]$ snakemake -np
Job counts:
        count   jobs
        1       all
        124     index_bam_females
        96      index_bam_males
        49      index_merged_bam_females
        36      index_merged_bam_males
        124     map_to_sex_specific_refs_females
        49      map_to_sex_specific_refs_males
        49      merge_bams_females
        36      merge_bams_males
        564
# Re-run snakemake
[amtarave@agave1:~/projects/Kenya_sequencing/whole_genome]$ sbatch snakemake_whole_genome_slurm.sh
Submitted batch job 147112

# Two jobs have failed so far. But the error is that there are missing files
[amtarave@agave1:~/projects/Kenya_sequencing/whole_genome]$ cat slurm-147133.out
Missing files after 5 seconds:
/scratch/amtarave/Kenya_agave/whole_genome/trimmed_fastqs/A63_H3HCNDMXX_L001_trimmed_R1.fastq.gz
/scratch/amtarave/Kenya_agave/whole_genome/trimmed_fastqs/A63_H3HCNDMXX_L001_trimmed_R2.fastq.gz
[amtarave@agave1:~/projects/Kenya_sequencing/whole_genome]$ cat slurm-147113.out
Missing files after 5 seconds:
/scratch/amtarave/Kenya_agave/whole_genome/trimmed_fastqs/A40_H3HCNDMXX_L002_trimmed_R1.fastq.gz
/scratch/amtarave/Kenya_agave/whole_genome/trimmed_fastqs/A40_H3HCNDMXX_L002_trimmed_R2.fastq.gz

# The weird thing is that these files are not missing.
[amtarave@agave1:~/projects/Kenya_sequencing/whole_genome]$ ls /scratch/amtarave/Kenya_agave/whole_genome/trimmed_fastqs/A63_H3HCNDMXX_L001_trimmed_*
/scratch/amtarave/Kenya_agave/whole_genome/trimmed_fastqs/A63_H3HCNDMXX_L001_trimmed_R1.fastq.gz
/scratch/amtarave/Kenya_agave/whole_genome/trimmed_fastqs/A63_H3HCNDMXX_L001_trimmed_R2.fastq.gz

[amtarave@agave1:~/projects/Kenya_sequencing/whole_genome]$ ls /scratch/amtarave/Kenya_agave/whole_genome/trimmed_fastqs/A40_H3HCNDMXX_L002_trimmed_*
/scratch/amtarave/Kenya_agave/whole_genome/trimmed_fastqs/A40_H3HCNDMXX_L002_trimmed_R1.fastq.gz
/scratch/amtarave/Kenya_agave/whole_genome/trimmed_fastqs/A40_H3HCNDMXX_L002_trimmed_R2.fastq.gz

# To address this issue, I will add a latency wait flag and re-run snakemake.
#  --latency-wait 60

[amtarave@agave1:~/projects/Kenya_sequencing/whole_genome]$ sbatch snakemake_whole_genome_slurm.sh
Submitted batch job 157939

## This started to work but I got the out of memory error again
A82_H3GWKDMXX_L002
A86_H57HCDMXX_L002
A86_H57HCDMXX_L001
A77_H57HCDMXX_L002
A82_H3GWKDMXX_L001
A71_H5G7VDMXX_L001
A76_H5G7VDMXX_L001
A77_H57HCDMXX_L001
A29_H3HCNDMXX_L001
A29_H3HCNDMXX_L002
A50_H5G7VDMXX_L001
A31_H3HCVDMXX_L001
A85_H3HFLDMXX_L002
A65_H3HCVDMXX_L001
A76_H5G7VDMXX_L002
A65_H3HCVDMXX_L002
A31_H3HCVDMXX_L002
A87_H57HCDMXX_L001
A56_H3GWKDMXX_L001
A87_H57HCDMXX_L002
A7_H3GWKDMXX_L001
A17_H3HCVDMXX_L002
A56_H3GWKDMXX_L002
A17_H3HCVDMXX_L001
A30_H3GWKDMXX_L001
A7_H3GWKDMXX_L002
...

# to make things easier i just deleted all of the generated bam files and i am re-running
# the snake file with --nodes=1 in the sbatch
[amtarave@agave1:~/projects/Kenya_sequencing/whole_genome]$ sbatch snakemake_whole_genome_slurm.sh
Submitted batch job 194118

# The following ran out of memory
A60_H5G7VDMXX_L002
A21_H3GWKDMXX_L002
A79_H3HCNDMXX_L002
A93_H3HCVDMXX_L001
A95_H57HCDMXX_L002
A95_H57HCDMXX_L001
A72_H5G7VDMXX_L002
A21_H3GWKDMXX_L001
A40_H3HCNDMXX_L002
A32_H3GWKDMXX_L002
A67_H3HCVDMXX_L002
A69_H5G7VDMXX_L001
A40_H3HCNDMXX_L001
A63_H3HCNDMXX_L001
A24_H575TDMXX_L002
A69_H5G7VDMXX_L002
A63_H3HCNDMXX_L002
A25_H57HCDMXX_L001
A25_H57HCDMXX_L002
A24_H575TDMXX_L001
A23_H57HCDMXX_L002
A61_H3HFLDMXX_L002
A67_H3HCVDMXX_L001
A34_H3HFLDMXX_L001
A32_H3GWKDMXX_L001
A34_H3HFLDMXX_L002
A61_H3HFLDMXX_L001
A38_H57HCDMXX_L002
A38_H57HCDMXX_L001
A23_H57HCDMXX_L001
A68_H3HFLDMXX_L001
A51_H575TDMXX_L001
A83_H3GWKDMXX_L001
A89_H3HCNDMXX_L001
A75_H3HFLDMXX_L001
A97_H57HCDMXX_L001
A97_H57HCDMXX_L002
A90_H3HCVDMXX_L002
A13_H3GWKDMXX_L002
A53_H3GWKDMXX_L002


# I am going to run an interactive session and map one bam file that failed to see if it
#	will complete this way. Also I can then see if these files are the same or if one is
#	corrupted
[amtarave@agave1:~]$ interactive -n 8 -N 1 --mem-per-cpu=16000
[amtarave@cg15-13:/scratch/amtarave/Kenya_agave/whole_genome/test_bam]$ pwd
/scratch/amtarave/Kenya_agave/whole_genome/test_bam

[amtarave@cg15-13:/scratch/amtarave/Kenya_agave/whole_genome/test_bam]$ source activate kenya
bwa mem -t 4 -R '@RG\tID:A21_H3GWKDMXX_L001\tSM:A21_whole_genome\tLB:A21_whole_genome\tPU:H3GWKDMXX.1\tPL:Illumina' /home/amtarave/projects/Kenya_sequencing/whole_genome/refs/Ref_GRCh38_Y_PARsMasked.fa /scratch/amtarave/Kenya_agave/whole_genome/trimmed_fastqs/A21_H3GWKDMXX_L001_trimmed_R1.fastq.gz /scratch/amtarave/Kenya_agave/whole_genome/trimmed_fastqs/A21_H3GWKDMXX_L001_trimmed_R2.fastq.gz | samtools fixmate -O bam - - | samtools sort -O bam -o A21_H3GWKDMXX_L001_test.bam

# This finished running and the file is not the same size as the bam file ran with snakemake

### I am just running the mapping for males for now
(kenya) [amtarave@cg5-14:~/projects/Kenya_sequencing/whole_genome]$ snakemake -np
Job counts:
        count   jobs
        1       all
        37      map_to_sex_specific_refs_males
        38

[amtarave@agave1:~/projects/Kenya_sequencing/whole_genome]$ sbatch snakemake_whole_genome_slurm.sh
Submitted batch job 203451

# got the OUT_OF_MEMORY error again but the bam files are good (mapped a sample
#	that "failed" and it is exactly the same as mapping in a regular sh script)..
#	apparently it is an Agave error not an error in my script.

# Add time option in config.json and add to sh script that runs snakemake.
# Dry run snakemake (I added female mapping and male indexing)
(kenya) [amtarave@cg15-3:~/projects/Kenya_sequencing/whole_genome]$ snakemake -np
Job counts:
        count   jobs
        1       all
        94      index_bam_males
        108     map_to_sex_specific_refs_females
        203


### Index female bams and merge male bams ###
# Dry run
(kenya) [amtarave@cg17-13:~/projects/Kenya_sequencing/whole_genome]$ snakemake -np

Job counts:
        count   jobs
        1       all
        124     index_bam_females
        36      merge_bams_males
        161

### Merge female bams index all merged bams ###
(kenya) [amtarave@cg15-3:~/projects/Kenya_sequencing/whole_genome]$ snakemake -np

Job counts:
        count   jobs
        1       all
        49      index_merged_bam_females
        35      index_merged_bam_males
        49      merge_bams_females
        134

#######
#Step 5
#######
# get stats on bam files #
# add to snakefile and run
(kenya) [amtarave@cg17-13:~/projects/Kenya_sequencing/whole_genome]$ snakemake -np
Job counts:
        count   jobs
        1       all
        85      bam_stats
        86

[amtarave@agave1:~/projects/Kenya_sequencing/whole_genome]$ sbatch snakemake_whole_genome_slurm.sh
Submitted batch job 207528


# Step 6

# Use gatk haplotypecaller to call snps for individuals
# This will require me to install gatk to my conda environment
conda install -c bioconda gatk
# this is gatk 3.8

# Then I need to incorporate the code into my snakefiles. I will start with the exome
# data and also call the monomorphic sites.
# Before doing the above, I need to mark duplicates in my bam files .
# Below will mark dups, index, call snps, and joint genotype for 10 females.
(kenya) [amtarave@cg15-3:~/projects/Kenya_sequencing/whole_genome]$ snakemake -np

Job counts:
        count   jobs
        1       all
        10      gatk_gvcf
        1       gatk_joint_genotype
        10      index_mkdup_bam
        10      mkdup_bam_stats
        10      picard_mkdups
        42

# joint genotype failed because I did not give java enough memory.
# Fix command and re-run (re-ran on 6/3/18 at 10am)
(kenya) [amtarave@fs1-1:~/projects/Kenya_sequencing/exome]$ snakemake -np
Job counts:
        count   jobs
        1       all
        1       gatk_joint_genotype
        2
[amtarave@ocotillo1:~/projects/Kenya_sequencing/exome]$ sbatch snakemake_exome_slurm.sh
Submitted batch job 1315012

# I found a forum on GATKs web site where someone wanted to figure out how to select SNP
# and reference monomorphic sites from a large call set.
# https://gatkforums.broadinstitute.org/gatk/discussion/8966/select-snp-and-reference-monomorphic-sites-from-a-large-callset
# This may be something we can use to get the monomorphic sites directly from the all sites
# VCF.
# their commands:
-T SelectVariants -R ref.fa \
--selectTypeToExclude INDEL --selectTypeToExclude MIXED --selectTypeToExclude MNP --selectTypeToExclude SYMBOLIC \
--selectTypeToInclude NO_VARIATION --maxNOCALLfraction 0 --excludeFiltered  -V $VCF  -o $OUT.Monomorphic.vcf.gz

-T SelectVariants  -R ref.fa \
--selectTypeToExclude INDEL --selectTypeToExclude MIXED --selectTypeToExclude MNP --selectTypeToExclude SYMBOLIC \
--selectTypeToInclude SNP  --maxNOCALLfraction 0 --excludeFiltered  --restrictAllelesTo BIALLELIC \
-V $VCF -o $OUT.SNP.vcf.gz

-T CombineVariants -R ref.fa -V $OUT.Monomorphic.vcf.gz -V $OUT.SNP.vcf.gz  -o $OUT.Clean.vcf.gz --assumeIdenticalSamples

# I think I want to have a monomorphic vcf and a snps vcf separately
#-T SelectVariants --selectTypeToInclude NO_VARIATION

# Dry run
(kenya) [amtarave@fs1-15:~/projects/Kenya_sequencing/exome]$ snakemake -np
Job counts:
        count   jobs
        1       all
        1       combine_sites
        1       select_SNPs
        3

# Submit snakemake job to get snps and combine with monomorphic sites
[amtarave@ocotillo2:~/projects/Kenya_sequencing/exome]$ sbatch snakemake_exome_slurm.sh
Submitted batch job 1340585


# to get number of monomorphic sites:
# 0/0 for every individual at a site is a monomorphic site. Also 1/1 at a site is
# a monomorphic site. But *.monomorphic.vcf is only 0/0 instances. *SNPs.VCF
# has 1/1 sites. To get number of monomorphic sites that are 1/1 use vcftools
VCFtools - v0.1.12b
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
        --vcf pilot_samples.raw.SNPs.vcf
        --non-ref-af 1
        --out test.vcf
        --recode

After filtering, kept 10 out of 10 Individuals
Outputting VCF file...
After filtering, kept 9564 out of a possible 621265 Sites
Run Time = 6.00 seconds

[amtarave@fs1-4:/scratch/amtarave/Kenya/exome/vcfs]$ vcftools --vcf pilot_samples.raw.monomorphic.vcf --non-ref-af 1

VCFtools - v0.1.12b
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
        --vcf pilot_samples.raw.monomorphic.vcf
        --non-ref-af 1

After filtering, kept 10 out of 10 Individuals
After filtering, kept 130650853 out of a possible 130650853 Sites
Run Time = 998.00 seconds


# Re-run to get vcf with no filters.

[amtarave@ocotillo2:~/projects/Kenya_sequencing/exome]$ sbatch snakemake_exome_slurm.sh 
Submitted batch job 1381578

# Run pilot analysis on WGS data:
[amtarave@agave1:~/projects/Kenya_sequencing/whole_genome]$ sbatch snakemake_whole_genome_slurm.sh
Submitted batch job 515412

# Re-run pilot WGS with more memory for gvcf step
[amtarave@cg8-8:~/projects/Kenya_sequencing/whole_genome]$ sbatch snakemake_whole_genome_slurm.sh
Submitted batch job 538230