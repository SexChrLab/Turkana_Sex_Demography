# Environment: CustomFilterVCF
import os

configfile: "hard_filter_variants_config.json"
rule all:
    # Extra stats
    input:
        "stats/autosomes_females_males_annotations_extra.txt",
        "stats/chrX_females_annotations_extra.txt",
        "stats/chrX_males_annotations_extra.txt"
    # chrX females
    input: #apply hard filter
        config["chrX_females_vcf_custom_filter"]
    input: #extract annotations
        "stats/chrX_females_annotations.txt",
        "stats/chrX_females_annotations_AN_plot.png",
        "stats/chrX_females_annotations_QD_plot.png",
        "stats/chrX_females_annotations_DP_plot.png",
        "stats/chrX_females_annotations_MQ_plot.png",
        "stats/chrX_females_annotations_DP_statistics.txt"
    # End chrX females

    # chrX males
    input: #apply hard filter
        config["chrX_males_vcf_custom_filter"]
    input: #extract annotations
        "stats/chrX_males_annotations.txt",
        "stats/chrX_males_annotations_AN_plot.png",
        "stats/chrX_males_annotations_QD_plot.png",
        "stats/chrX_males_annotations_DP_plot.png",
        "stats/chrX_males_annotations_MQ_plot.png",
        "stats/chrX_males_annotations_DP_statistics.txt"
    # End chrX males

    # Autosomes
    input: #apply hard filter on each chromosome of the autosomes
        expand(os.path.join(config["custom_filter_vcf_directory"], "{chr}", "{chr}.gatk.called.custom.filter.autosomes.females.males.vcf.gz"), chr=config["autosomes"])
    input: #apply hard filter for autosomes
        config["autosomes_vcf_custom_filter"]
    input: #extract annotations
        "stats/autosomes_females_males_annotations.txt",
        "stats/autosomes_females_males_annotations_AN_plot.png",
        "stats/autosomes_females_males_annotations_QD_plot.png",
        "stats/autosomes_females_males_annotations_DP_plot.png",
        "stats/autosomes_females_males_annotations_MQ_plot.png",
        "stats/autosomes_females_males_annotations_DP_statistics.txt"
    # End autosomes

# -------------------
# rules for autosomes
# -------------------
rule extract_stats_autosomes_females_males:
    input:
        vcf = config["autosomes_vcf"]
    output:
        stats = "stats/autosomes_females_males_annotations.txt"
    params:
        script = config["extract_stats_from_vcf_script"],
        AN = "AN",
        QD = "QD",
        MQ = "MQ",
        DP = "DP"
    shell:
        """
        python {params.script} {params.AN} {params.QD} {params.MQ} {params.DP} --vcf {input.vcf} --outfile {output.stats}
        """

rule plot_stats:
    input:
        "stats/autosomes_females_males_annotations.txt"
    output:
        AN_plot = "stats/autosomes_females_males_annotations_AN_plot.png",
        QD_plot = "stats/autosomes_females_males_annotations_QD_plot.png",
        DP_plot = "stats/autosomes_females_males_annotations_DP_plot.png",
        MQ_plot = "stats/autosomes_females_males_annotations_MQ_plot.png",
        DP_statistics = "stats/autosomes_females_males_annotations_DP_statistics.txt"
    params:
        script = config["plot_stats_script"]
    shell:
        """
        Rscript {params.script} {input} {output.AN_plot} {output.QD_plot} {output.DP_plot} {output.MQ_plot} {output.DP_statistics}
        """

rule custom_filter_autosomes_females_males:
    input:
        ref = config["Ref_GRCh38_Y_HardMasked"],
        vcf = os.path.join(config["raw_vcf_directory"], "{chr}", "{chr}.gatk.called.raw.autosomes.females.males.vcf.gz")
    output:
        os.path.join(config["custom_filter_vcf_directory"], "{chr}", "{chr}.gatk.called.custom.filter.autosomes.females.males.vcf.gz")
    params:
        gatk = config["gatk3_path"],
        AN = 150,
        MQ = 40.0,
        QD = 5.0,
        DP = 50000.0
    shell:
        """java -jar {params.gatk} """
        """-R {input.ref} """
        """-T SelectVariants """
        """-V {input.vcf} """
        """-o {output} """
        """-select "AN >= {params.AN} && MQ > {params.MQ} && QD > {params.QD} && DP <= {params.DP}" """

rule extract_stats_autosomes_females_males_extra:
    input:
        vcf = config["autosomes_vcf"]
    output:
        stats = "stats/autosomes_females_males_annotations_extra.txt"
    params:
        script = config["extract_stats_from_vcf_script"],
        FS = "FS",
        ReadPosRankSum = "ReadPosRankSum",
        SOR = "SOR"
    shell:
        """
        python {params.script} {params.FS} {params.ReadPosRankSum} {params.SOR} --vcf {input.vcf} --outfile {output.stats}
        """

# ----------------------
# rules for chrX females
# ----------------------
rule extract_stats_chrX_females:
    input:
        vcf = config["chrX_females_vcf"]
    output:
        stats = "stats/chrX_females_annotations.txt"
    params:
        script = config["extract_stats_from_vcf_script"],
        AN = "AN",
        QD = "QD",
        MQ = "MQ",
        DP = "DP"
    shell:
        """
        python {params.script} {params.AN} {params.QD} {params.MQ} {params.DP} --vcf {input.vcf} --outfile {output.stats}
        """

rule plot_stats_chrX_females:
    input:
        "stats/chrX_females_annotations.txt"
    output:
        AN_plot = "stats/chrX_females_annotations_AN_plot.png",
        QD_plot = "stats/chrX_females_annotations_QD_plot.png",
        DP_plot = "stats/chrX_females_annotations_DP_plot.png",
        MQ_plot = "stats/chrX_females_annotations_MQ_plot.png",
        DP_statistics = "stats/chrX_females_annotations_DP_statistics.txt"
    params:
        script = config["plot_stats_script"]
    shell:
        """
        Rscript {params.script} {input} {output.AN_plot} {output.QD_plot} {output.DP_plot} {output.MQ_plot} {output.DP_statistics}
        """

rule custom_filter_chrX_females:
    input:
        ref = config["Ref_GRCh38_Y_HardMasked"],
        vcf = config["chrX_females_vcf"]
    output:
        config["chrX_females_vcf_custom_filter"]
    params:
        gatk = config["gatk3_path"],
        AN = 75,
        MQ = 40.0,
        QD = 5.0,
        DP = 5000.0
    shell:
        """java -jar {params.gatk} """
        """-R {input.ref} """
        """-T SelectVariants """
        """-V {input.vcf} """
        """-o {output} """
        """-select "AN >= {params.AN} && MQ > {params.MQ} && QD > {params.QD} && DP <= {params.DP}" """

rule extract_stats_chrX_females_extra:
    input:
        vcf = config["chrX_females_vcf"]
    output:
        stats = "stats/chrX_females_annotations_extra.txt"
    params:
        script = config["extract_stats_from_vcf_script"],
        FS = "FS",
        ReadPosRankSum = "ReadPosRankSum",
        SOR = "SOR"
    shell:
        """
        python {params.script} {params.FS} {params.ReadPosRankSum} {params.SOR} --vcf {input.vcf} --outfile {output.stats}
        """

# --------------------
# rules for chrX males
# --------------------
rule extract_stats_chrX_males:
    input:
        vcf = config["chrX_males_vcf"]
    output:
        stats = "stats/chrX_males_annotations.txt"
    params:
        script = config["extract_stats_from_vcf_script"],
        AN = "AN",
        QD = "QD",
        MQ = "MQ",
        DP = "DP"
    shell:
        """
        python {params.script} {params.AN} {params.QD} {params.MQ} {params.DP} --vcf {input.vcf} --outfile {output.stats}
        """

rule plot_stats_chrX_males:
    input:
        "stats/chrX_males_annotations.txt"
    output:
        AN_plot = "stats/chrX_males_annotations_AN_plot.png",
        QD_plot = "stats/chrX_males_annotations_QD_plot.png",
        DP_plot = "stats/chrX_males_annotations_DP_plot.png",
        MQ_plot = "stats/chrX_males_annotations_MQ_plot.png",
        DP_statistics = "stats/chrX_males_annotations_DP_statistics.txt"
    params:
        script = config["plot_stats_script"]
    shell:
        """
        Rscript {params.script} {input} {output.AN_plot} {output.QD_plot} {output.DP_plot} {output.MQ_plot} {output.DP_statistics}
        """

rule custom_filter_chrX_males:
    input:
        ref = config["Ref_GRCh38_Y_PARsMasked"],
        vcf = config["chrX_males_vcf"]
    output:
        config["chrX_males_vcf_custom_filter"]
    params:
        gatk = config["gatk3_path"],
        AN = 30,
        MQ = 50.0,
        QD = 15.0,
        DP = 2000.0
    shell:
        """java -jar {params.gatk} """
        """-R {input.ref} """
        """-T SelectVariants """
        """-V {input.vcf} """
        """-o {output} """
        """-select "AN >= {params.AN} && MQ > {params.MQ} && QD > {params.QD} && DP <= {params.DP}" """

rule extract_stats_chrX_males_extra:
    input:
        vcf = config["chrX_males_vcf"]
    output:
        stats = "stats/chrX_males_annotations_extra.txt"
    params:
        script = config["extract_stats_from_vcf_script"],
        FS = "FS",
        ReadPosRankSum = "ReadPosRankSum",
        SOR = "SOR"
    shell:
        """
        python {params.script} {params.FS} {params.ReadPosRankSum} {params.SOR} --vcf {input.vcf} --outfile {output.stats}
        """
