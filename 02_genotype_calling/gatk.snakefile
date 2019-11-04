# This is a snakefile for GATK for Kimberly:
# - Joint call across 30 female placentas and 30 male placentas for autosomes
# - Joint call across 30 female placentas for chrX
# - Joint call across 30 male placentas for chrX
# - Joint call across 30 male placentas for chrY

import os

configfile: "gatk_config.json"

# Tool paths:
gatk_path = "/home/tphung3/softwares/gatk-4.1.0.0/gatk"
gatk3_path = "/home/tphung3/softwares/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar"

rule all:
    input: #Generate gvcf files
        expand("gvcfs/{chrm}/{chrm}.{sample}.female.diploid.g.vcf.gz", chrm=config["autosomes"], sample=config["females"]),
        expand("gvcfs/{chrm}/{chrm}.{sample}.male.diploid.g.vcf.gz", chrm=config["autosomes"], sample=config["males"]),
        expand("gvcfs/chrX/sexchr.chrX.{sample}.female.diploid.g.vcf.gz", sample=config["females"]),
        expand("gvcfs/{chrm}/sexchr.{chrm}.{sample}.male.haploid.g.vcf.gz", chrm=config["sex_chr_haploid"], sample=config["males"])


# --------------------------------
# Create symlink for the bam files
# --------------------------------
rule mk_symlink_for_bams_females:
    input:
        bams = os.path.join(config["bams_directory"], "{sample}.GRCh38_Ymasked.sorted.merged.mkdup.bam"),
        bais = os.path.join(config["bams_directory"], "{sample}.GRCh38_Ymasked.sorted.merged.mkdup.bam.bai")
    output:
        bams_symlink = "processed_bams/{sample}.GRCh38_Ymasked.sorted.merged.mkdup.bam",
        bais_symlink = "processed_bams/{sample}.GRCh38_Ymasked.sorted.merged.mkdup.bam.bai"
    shell:
        """
        ln -s {input.bams} {output.bams_symlink};
        ln -s {input.bais} {output.bais_symlink}
        """

rule mk_symlink_for_bams_males:
    input:
        bams = os.path.join(config["bams_directory"], "{sample}.GRCh38_minusYPARs.sorted.merged.mkdup.bam"),
        bais = os.path.join(config["bams_directory"], "{sample}.GRCh38_minusYPARs.sorted.merged.mkdup.bam.bai")
    output:
        bams_symlink = "processed_bams/{sample}.GRCh38_minusYPARs.sorted.merged.mkdup.bam",
        bais_symlink = "processed_bams/{sample}.GRCh38_minusYPARs.sorted.merged.mkdup.bam.bai"
    shell:
        """
        ln -s {input.bams} {output.bams_symlink};
        ln -s {input.bais} {output.bais_symlink}
        """

# -------------------------------------------
# Generate the gVCF files for each individual
# Autosomes (females + males)
# X chrosome diploid (females)
# X and Y chromosomes haploid (males)
# -------------------------------------------
rule gatk_gvcf_females_autosomes:
    input:
        ref = config["Ref_GRCh38_Y_HardMasked"],
        bam = "processed_bams/{sample}.GRCh38_Ymasked.sorted.merged.mkdup.bam",
        bai = "processed_bams/{sample}.GRCh38_Ymasked.sorted.merged.mkdup.bam.bai"
    output:
        "gvcfs/{chrm}/{chrm}.{sample}.female.diploid.g.vcf.gz"
    params:
        gatk = gatk_path,
        chrm_n = "{chrm}"
    shell:
        "{params.gatk} "
        "HaplotypeCaller -R {input.ref} -I {input.bam} -L {params.chrm_n} "
        "--emit-ref-confidence GVCF --output {output}"

rule gatk_gvcf_males_autosomes:
    input:
        ref = config["Ref_GRCh38_Y_HardMasked"],
        bam = "processed_bams/{sample}.GRCh38_minusYPARs.sorted.merged.mkdup.bam",
        bai = "processed_bams/{sample}.GRCh38_minusYPARs.sorted.merged.mkdup.bam.bai"
    output:
        "gvcfs/{chrm}/{chrm}.{sample}.male.diploid.g.vcf.gz"
    params:
        gatk = gatk_path,
        chrm_n = "{chrm}"
    shell:
        "{params.gatk} "
        "HaplotypeCaller -R {input.ref} -I {input.bam} -L {params.chrm_n} "
        "--emit-ref-confidence GVCF --output {output}"

rule gatk_gvcf_females_sexchr:
    input:
        ref = config["Ref_GRCh38_Y_HardMasked"],
        bam = "processed_bams/{sample}.GRCh38_Ymasked.sorted.merged.mkdup.bam",
        bai = "processed_bams/{sample}.GRCh38_Ymasked.sorted.merged.mkdup.bam.bai"
    output:
        "gvcfs/chrX/sexchr.chrX.{sample}.female.diploid.g.vcf.gz"
    params:
        gatk = gatk_path,
        chrm_n = "chrX"
    shell:
        "{params.gatk} "
        "HaplotypeCaller -R {input.ref} -I {input.bam} -L {params.chrm_n} "
        "--emit-ref-confidence GVCF --output {output}"

rule gatk_gvcf_males_sexchr:
    input:
        ref = config["Ref_GRCh38_Y_PARsMasked"],
        bam = "processed_bams/{sample}.GRCh38_minusYPARs.sorted.merged.mkdup.bam",
        bai = "processed_bams/{sample}.GRCh38_minusYPARs.sorted.merged.mkdup.bam.bai"
    output:
        "gvcfs/{chrm}/sexchr.{chrm}.{sample}.male.haploid.g.vcf.gz"
    params:
        gatk = gatk_path,
        chrm_n = "{chrm}"
    shell:
        "{params.gatk} "
        "HaplotypeCaller -R {input.ref} -I {input.bam} -L {params.chrm_n} "
        "--emit-ref-confidence GVCF -ploidy 1 --output {output}"
