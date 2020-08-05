import os

configfile: "gatk_config.json"

# Tool paths:
gatk_path = "/home/tphung3/softwares/gatk-4.1.0.0/gatk"
gatk3_path = "/home/tphung3/softwares/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar"

rule all:
    # -------------------------------------------------------------------------
    # Start rules for processing sex chromosomes in males
    input:
        expand("genotyped_vcfs/{chrm}/{chrm}.gatk.called.raw.sexchr.males.vcf.gz", chrm=config["sex_chr_haploid"])
    input:
        expand("gvcfs/{chrm}/sexchr.{chrm}.{sample}.male.haploid.g.vcf.gz", chrm=config["sex_chr_haploid"], sample=config["males"])
    # End rules for processing sex chromosomes in males
    # -------------------------------------------------------------------------
    # Start rules for processing X chromosome in females
    input:
        "genotyped_vcfs/chrX/chrX.gatk.called.raw.sexchr.females.vcf.gz"
    input:
        expand("gvcfs/chrX/sexchr.chrX.{sample}.female.diploid.g.vcf.gz", sample=config["females"])
    # End rule for processing X chromosome in females
    # -------------------------------------------------------------------------
    # Start rules for processing autosomes
    input: #combine chr1 to chr22
        "genotyped_vcfs/autosomes/gatk.called.raw.autosomes.females.males.vcf.gz"
    input: #genotype for Autosomes
        expand("genotyped_vcfs/{chrm}/{chrm}.gatk.called.raw.autosomes.females.males.all.sites.vcf.gz", chrm=config["autosomes"])
    input: #Generate gvcf files for males autotomes diploid
        expand("gvcfs/{chrm}/{chrm}.{sample}.male.diploid.g.vcf.gz", chrm=config["autosomes"], sample=config["males"])
    input: #Generate gvcf files for females autotomes diploid
        expand("gvcfs/{chrm}/{chrm}.{sample}.female.diploid.g.vcf.gz", chrm=config["autosomes"], sample=config["females"])
    # End rules for processing autosomes
    # -------------------------------------------------------------------------



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

# --------------------------------------------------------------------
# Joint call across all samples (both females and males) for autosomes
# --------------------------------------------------------------------
rule gatk_combinegvcfs_autosomes_females_males:
    input:
        ref = config["Ref_GRCh38_Y_HardMasked"],
        female_gvcfs = lambda wildcards: expand(
			"gvcfs/{chrm}/{chrm}.{sample}.female.diploid.g.vcf.gz",
			sample=config["females"],
			chrm=config["autosomes"]),
        male_gvcfs = lambda wildcards: expand(
			"gvcfs/{chrm}/{chrm}.{sample}.male.diploid.g.vcf.gz",
			sample=config["males"],
			chrm=config["autosomes"])

    params:
        gatk = gatk_path,
        chrm_n = "{chrm}"

    output:
        "combined_gvcfs/{chrm}/{chrm}.gatk.combinegvcf.autosomes.females.males.g.vcf.gz"

    run:
        variant_files = []
        for i in input.female_gvcfs:
	           variant_files.append("--variant " + i)
        for i in input.male_gvcfs:
	           variant_files.append("--variant " + i)
        variant_files = " ".join(variant_files)
        print(
	       """{params.gatk} --java-options "-Xmx10g" """
	          """CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chrm_n} -O {output}"""
        )
        shell(
	       """{params.gatk} --java-options "-Xmx10g" """
	          """CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chrm_n} -O {output}""")

rule gatk_genotypegvcf_autosomes_females_males:
    input:
        ref = config["Ref_GRCh38_Y_HardMasked"],
        gvcf = "combined_gvcfs/{chrm}/{chrm}.gatk.combinegvcf.autosomes.females.males.g.vcf.gz"
    output:
        "genotyped_vcfs/{chrm}/{chrm}.gatk.called.raw.autosomes.females.males.all.sites.vcf.gz"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} --java-options "-Xmx10g" """
        """GenotypeGVCFs -R {input.ref} -V {input.gvcf} --include-non-variant-sites -O {output} """

rule combine_autosomes: #after joint genotypying, combine chr1 to chr22
    input:
        ref = config["Ref_GRCh38_Y_HardMasked"],
        vcfs = lambda wildcards: expand(
			"genotyped_vcfs/{chrm}/{chrm}.gatk.called.raw.autosomes.females.males.vcf.gz",
			chrm=config["autosomes"])
    output:
        "genotyped_vcfs/autosomes/gatk.called.raw.autosomes.females.males.vcf.gz"
    params:
        gatk = gatk3_path
    run:
        variant_files = []
        for i in input.vcfs:
	           variant_files.append("--variant " + i)
        variant_files = " ".join(variant_files)
        shell(
        """java -cp {params.gatk} org.broadinstitute.gatk.tools.CatVariants """
        """-R {input.ref} {variant_files} -out {output} -assumeSorted""" )

# Females
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

rule gatk_combinegvcfs_females_sexchr:
    input:
        ref = config["Ref_GRCh38_Y_HardMasked"],
        gvcfs = lambda wildcards: expand(
			"gvcfs/chrX/sexchr.chrX.{sample}.female.diploid.g.vcf.gz",
			sample=config["females"])
    params:
        gatk = gatk_path,
        chrm_n = "chrX"

    output:
        "combined_gvcfs/chrX/chrX.gatk.combinegvcf.sexchr.females.g.vcf.gz"

    run:
        variant_files = []
        for i in input.gvcfs:
	           variant_files.append("--variant " + i)
        variant_files = " ".join(variant_files)
        print(
	       """{params.gatk} --java-options "-Xmx10g" """
	          """CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chrm_n} -O {output}"""
        )
        shell(
	       """{params.gatk} --java-options "-Xmx10g" """
	          """CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chrm_n} -O {output}""")

rule gatk_genotypegvcf_females_sexchr:
    input:
        ref = config["Ref_GRCh38_Y_HardMasked"],
        gvcf = "combined_gvcfs/chrX/chrX.gatk.combinegvcf.sexchr.females.g.vcf.gz"
    output:
        "genotyped_vcfs/chrX/chrX.gatk.called.raw.sexchr.females.vcf.gz"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} --java-options "-Xmx10g" """
        """GenotypeGVCFs -R {input.ref} -V {input.gvcf} -O {output} """

# Males
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

rule gatk_combinegvcfs_sexchr_males:
    input:
        ref = config["Ref_GRCh38_Y_PARsMasked"],
        gvcfs = lambda wildcards: expand(
			"gvcfs/{chrm}/sexchr.{chrm}.{sample}.male.haploid.g.vcf.gz",
			sample=config["males"],
			chrm=config["sex_chr_haploid"])
    params:
        gatk = gatk_path,
        chrm_n = "{chrm}"

    output:
        "combined_gvcfs/{chrm}/{chrm}.gatk.combinegvcf.sexchr.males.g.vcf.gz"

    run:
        variant_files = []
        for i in input.gvcfs:
	           variant_files.append("--variant " + i)
        variant_files = " ".join(variant_files)
        print(
	       """{params.gatk} --java-options "-Xmx10g" """
	          """CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chrm_n} -O {output}"""
        )
        shell(
	       """{params.gatk} --java-options "-Xmx10g" """
	          """CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chrm_n} -O {output}""")

rule gatk_genotypegvcf_sexchr_males:
    input:
        ref = config["Ref_GRCh38_Y_PARsMasked"],
        gvcf = "combined_gvcfs/{chrm}/{chrm}.gatk.combinegvcf.sexchr.males.g.vcf.gz"
    output:
        "genotyped_vcfs/{chrm}/{chrm}.gatk.called.raw.sexchr.males.vcf.gz"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} --java-options "-Xmx10g" """
        """GenotypeGVCFs -R {input.ref} -V {input.gvcf} -O {output} """
