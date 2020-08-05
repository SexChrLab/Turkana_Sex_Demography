library(ggplot2)
library(ggpubr)

autosomes_sfs = read.table("c://Users/tuyen/Documents/postdoc_asu/projects/Turkana_Demography/05_generate_sfs/results/autosomes_females_males_sfs_112719.txt")

colnames(autosomes_sfs) = c("bins", "count")

autosomes_sfs$bins = as.character(autosomes_sfs$bins)
autosomes_sfs$bins = factor(autosomes_sfs$bins, levels = unique(autosomes_sfs$bins))

png("c://Users/tuyen/Documents/postdoc_asu/projects/Turkana_Demography/05_generate_sfs/results/autosomes_females_males_sfs.png", width = 13, height = 6, units = "in", res = 300)
ggplot(data=autosomes_sfs, aes(x=bins, y=count)) + 
  geom_bar(stat = "identity", position=position_dodge(), colour="black", fill="gray40") + 
  theme_bw() +
  theme(axis.text.x = element_text(size=8, colour="black", angle=50, hjust=1), axis.text.y = element_text(size=12, colour = "black"), axis.title.y = element_text(size=18), axis.line=element_line(colour = "black"),plot.title = element_text(hjust = 0.5, size=18), legend.title = element_blank(), panel.background = element_blank(), legend.text = element_text(size=14), legend.spacing.x = unit(0.15, 'cm')) +
  labs(y="Counts", x="Allele frequency bin", title="Autosomes")
dev.off()

# ------------
# chrX females
# ------------

chrX_females_sfs = read.table("/scratch/tphung3/Turkana_Demography/05_generate_sfs/results/chrX/chrX.females.sfs.txt", header=T)

colnames(chrX_females_sfs) = c("bins", "count")

chrX_females_sfs$bins = as.character(chrX_females_sfs$bins)
chrX_females_sfs$bins = factor(chrX_females_sfs$bins, levels = unique(chrX_females_sfs$bins))

p1 = ggplot(data=chrX_females_sfs, aes(x=bins, y=count)) + 
  geom_bar(stat = "identity", position=position_dodge(), colour="black", fill="gray40") + 
  theme_bw() +
  theme(axis.text.x = element_text(size=8, colour="black", angle=50, hjust=1), axis.text.y = element_text(size=12, colour = "black"), axis.title.y = element_text(size=18), axis.line=element_line(colour = "black"),plot.title = element_text(hjust = 0.5, size=18), legend.title = element_blank(), panel.background = element_blank(), legend.text = element_text(size=14), legend.spacing.x = unit(0.15, 'cm')) +
  labs(y="Counts", x="Allele frequency bin", title="chrX females")

# ------------
# chrX males
# ------------

chrX_males_sfs = read.table("/scratch/tphung3/Turkana_Demography/05_generate_sfs/results/chrX/chrX.males.sfs.txt", header=T)

colnames(chrX_males_sfs) = c("bins", "count")

chrX_males_sfs$bins = as.character(chrX_males_sfs$bins)
chrX_males_sfs$bins = factor(chrX_males_sfs$bins, levels = unique(chrX_males_sfs$bins))

p2 = ggplot(data=chrX_males_sfs, aes(x=bins, y=count)) + 
  geom_bar(stat = "identity", position=position_dodge(), colour="black", fill="gray40") + 
  theme_bw() +
  theme(axis.text.x = element_text(size=8, colour="black", angle=50, hjust=1), axis.text.y = element_text(size=12, colour = "black"), axis.title.y = element_text(size=18), axis.line=element_line(colour = "black"),plot.title = element_text(hjust = 0.5, size=18), legend.title = element_blank(), panel.background = element_blank(), legend.text = element_text(size=14), legend.spacing.x = unit(0.15, 'cm')) +
  labs(y="Counts", x="Allele frequency bin", title="chrX males")

ggarrange(p1, p2, ncol = 2, nrow = 1)