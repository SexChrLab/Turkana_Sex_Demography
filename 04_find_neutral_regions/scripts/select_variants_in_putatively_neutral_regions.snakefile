# Environment: CustomFilterVCF
import os

configfile: "select_variants_in_putatively_neutral_regions_config.json"

rule all:
    input: #chrX females
        os.path.join(config["puatively_neutral_directory"], "chrX", "chrX.gatk.called.custom.filter.females.neutral.112719.vcf.gz")
    input: #chrX males
        os.path.join(config["puatively_neutral_directory"], "chrX", "chrX.gatk.called.custom.filter.males.neutral.112719.vcf.gz")
    input:
        expand(os.path.join(config["puatively_neutral_directory"], "{chr}", "{chr}.gatk.called.custom.filter.autosomes.females.males.neutral.112719.vcf.gz"), chr=config["autosomes"])

rule select_variants_in_neutral:
    input:
        ref = config["Ref_GRCh38_Y_HardMasked"],
        vcf = os.path.join(config["custom_filter_vcf_directory"], "{chr}", "{chr}.gatk.called.custom.filter.autosomes.females.males.vcf.gz"),
        bed = os.path.join(config["puatively_neutral_directory"], "{chr}", "{chr}_putatively_neutral_regions_112719.bed")
    output:
        os.path.join(config["puatively_neutral_directory"], "{chr}", "{chr}.gatk.called.custom.filter.autosomes.females.males.neutral.112719.vcf.gz")
    params:
        gatk = config["gatk3_path"]
    shell:
        """java -jar {params.gatk} """
        """-R {input.ref} """
        """-T SelectVariants """
        """-V {input.vcf} """
        """-L {input.bed} """
        """--restrictAllelesTo BIALLELIC """
        """-o {output} """

rule select_variants_in_neutral_chrX_females:
    input:
        ref = config["Ref_GRCh38_Y_HardMasked"],
        vcf = os.path.join(config["custom_filter_vcf_directory"], "gatk.called.custom.filter.chrX.females.vcf.gz"),
        bed = os.path.join(config["puatively_neutral_directory"], "chrX", "chrX_putatively_neutral_regions_112719.bed")
    output:
        os.path.join(config["puatively_neutral_directory"], "chrX", "chrX.gatk.called.custom.filter.females.neutral.112719.vcf.gz")
    params:
        gatk = config["gatk3_path"]
    shell:
        """java -jar {params.gatk} """
        """-R {input.ref} """
        """-T SelectVariants """
        """-V {input.vcf} """
        """-L {input.bed} """
        """--restrictAllelesTo BIALLELIC """
        """-o {output} """

rule select_variants_in_neutral_chrX_males:
    input:
        ref = config["Ref_GRCh38_Y_PARsMasked"],
        vcf = os.path.join(config["custom_filter_vcf_directory"], "gatk.called.custom.filter.chrX.males.vcf.gz"),
        bed = os.path.join(config["puatively_neutral_directory"], "chrX", "chrX_putatively_neutral_regions_112719.bed")
    output:
        os.path.join(config["puatively_neutral_directory"], "chrX", "chrX.gatk.called.custom.filter.males.neutral.112719.vcf.gz")
    params:
        gatk = config["gatk3_path"]
    shell:
        """java -jar {params.gatk} """
        """-R {input.ref} """
        """-T SelectVariants """
        """-V {input.vcf} """
        """-L {input.bed} """
        """--restrictAllelesTo BIALLELIC """
        """-o {output} """
