if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("SNPRelate")

library(gdsfmt)
library(SNPRelate)

setwd("/scratch/tphung3/Turkana_Demography/02_genotype_calling/genotyped_vcfs/autosomes")

pop_code = rep("turkana", 82)

vcf.fn = "gatk.called.raw.autosomes.females.males.vcf"
snpgdsVCF2GDS(vcf.fn, "autosomes_all_variants_WGS.gds", method="biallelic.only")
# Start file conversion from VCF to SNP GDS ...
# Method: exacting biallelic SNPs
# Number of samples: 82
# Parsing "gatk.called.raw.autosomes.females.males.vcf" ...
# import 23846790 variants.
# + genotype   { Bit2 82x23846790, 466.2M } *
#   Optimize the access efficiency ...
# Clean up the fragments of GDS file:
#   open the file 'autosomes_all_variants_WGS.gds' (579.3M)
# # of fragments: 2751
# save to 'autosomes_all_variants_WGS.gds.tmp'
# rename 'autosomes_all_variants_WGS.gds.tmp' (579.3M, reduced: 32.0K)
# # of fragments: 20
snpgdsSummary("autosomes_all_variants_WGS.gds")
# The file name: /scratch/tphung3/Turkana_Demography/02_genotype_calling/genotyped_vcfs/autosomes/autosomes_all_variants_WGS.gds 
# The total number of samples: 82 
# The total number of SNPs: 23846790 
# SNP genotypes are stored in SNP-major mode (Sample X SNP).
genofile <- snpgdsOpen("autosomes_all_variants_WGS.gds")
pca = snpgdsPCA(genofile, autosome.only=TRUE)
# Principal Component Analysis (PCA) on genotypes:
#   Excluding 0 SNP on non-autosomes
# Excluding 106,687 SNPs (monomorphic: TRUE, MAF: NaN, missing rate: NaN)
# Working space: 82 samples, 23,740,103 SNPs
# using 1 (CPU) core
# PCA:    the sum of all selected genotypes (0,1,2) = 3364049911
# CPU capabilities: Double-Precision SSE2
# Mon Nov 18 17:55:26 2019    (internal increment: 55492)
# [==================================================] 100%, completed, 59s 
# Mon Nov 18 17:56:25 2019    Begin (eigenvalues and eigenvectors)
# Mon Nov 18 17:56:25 2019    Done.
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
# [1] 1.58 1.54 1.45 1.42 1.41 1.39

# Make a data frame
sample.id = pca$sample.id
autosomes_df <- data.frame(sample.id = pca$sample.id,
                       pop = factor(pop_code)[match(pca$sample.id, sample.id)],
                       EV1 = pca$eigenvect[,1],    # the first eigenvector
                       EV2 = pca$eigenvect[,2],    # the second eigenvector
                       EV3 = pca$eigenvect[,3],
                       EV4 = pca$eigenvect[,4],
                       EV5 = pca$eigenvect[,5],
                       stringsAsFactors = FALSE)

pop = factor(pop_code)[match(pca$sample.id, sample.id)]

lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=pop, labels=lbls)

snpgdsClose(genofile)