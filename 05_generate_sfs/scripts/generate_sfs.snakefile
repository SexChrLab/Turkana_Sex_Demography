import os

configfile: "generate_sfs_config.json"

rule all:
    input:
        expand(os.path.join(config["sfs_directory"], "{chr}", "{chr}.males.sfs.txt"), chr=config["sex_chr"])
    input:
        expand(os.path.join(config["sfs_directory"], "{chr}", "{chr}.females.sfs.txt"), chr=config["sex_chr"])
    input:
        expand(os.path.join(config["sfs_directory"], "{chr}", "{chr}.sfs.txt"), chr=config["autosomes"])
rule generate_sfs_autosomes_females_males:
    input:
        os.path.join(config["vcf_directory"], "{chr}", "{chr}.gatk.called.custom.filter.autosomes.females.males.neutral.112719.vcf.gz")
    output:
        os.path.join(config["sfs_directory"], "{chr}", "{chr}.sfs.txt")
    params:
        script = config["generate_sfs_script"]
    shell:
        """
        python {params.script} --vcf_file {input} --sfs_all --sfs_all_out {output}
        """

rule generate_sfs_chrX_females:
    input:
        os.path.join(config["vcf_directory"], "{chr}", "{chr}.gatk.called.custom.filter.females.neutral.112719.vcf.gz")
    output:
        os.path.join(config["sfs_directory"], "{chr}", "{chr}.females.sfs.txt")
    params:
        script = config["generate_sfs_script"]
    shell:
        """
        python {params.script} --vcf_file {input} --sfs_all --sfs_all_out {output}
        """

rule generate_sfs_chrX_males:
    input:
        os.path.join(config["vcf_directory"], "{chr}", "{chr}.gatk.called.custom.filter.males.neutral.112719.vcf.gz")
    output:
        os.path.join(config["sfs_directory"], "{chr}", "{chr}.males.sfs.txt")
    params:
        script = config["generate_sfs_script"]
    shell:
        """
        python {params.script} --vcf_file {input} --sfs_all --sfs_all_out {output} --ploidy haploid
        """
