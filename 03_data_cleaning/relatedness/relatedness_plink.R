genofile <- snpgdsOpen("/scratch/tphung3/Turkana_Demography/02_genotype_calling/genotyped_vcfs/autosomes/autosomes_all_variants_WGS.gds")

# Estimate IBD coefficients
ibd <- snpgdsIBDMoM(genofile, maf=0.05, missing.rate=0.05, num.thread=2)

# Make a data.frame
ibd.coeff <- snpgdsIBDSelection(ibd)
head(ibd.coeff)

plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
     xlab="k0", ylab="k1")
lines(c(0,1), c(1,0), col="red", lty=2)