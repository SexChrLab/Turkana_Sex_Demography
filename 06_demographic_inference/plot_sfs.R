library(ggplot2)
library(ggpubr)

empirical = read.csv("c://Users/tuyen/Documents/postdoc_asu/projects/Turkana_Demography/05_generate_sfs/results/autosomes_females_males_sfs_112719.txt", header = FALSE, sep="\t")
model_1epoch = read.csv("c://Users/tuyen/Documents/postdoc_asu/projects/Turkana_Demography/06_demographic_inference/1Epoch_run_1_sfs_normalized_by_theta.txt", header = FALSE, sep="\t")
model_2epoch = read.csv("c://Users/tuyen/Documents/postdoc_asu/projects/Turkana_Demography/06_demographic_inference/2Epoch_run_10_sfs_normalized_by_theta.txt", header = FALSE, sep="\t")

empirical_labels = rep("Empirical", 82)
oneepoch_labels = rep("1Epoch", 82)
twoepoch_labels = rep("2Epoch", 82)

dat = data.frame(Bins = empirical[,1], Counts=c(empirical[,2], model_1epoch[,2], model_2epoch[,2]), category = c(empirical_labels, oneepoch_labels, twoepoch_labels))

dat$category = as.character(dat$category)
dat$category = factor(dat$category, levels = unique(dat$category))

ggplot(data=dat, aes(x=Bins, y=Counts, fill=category)) + 
  theme_bw() +
  geom_bar(stat = "identity", position=position_dodge(), colour="black") + 
  scale_fill_manual(values=c("gray40", "lightblue", "gold")) + 
  theme(axis.text.x = element_text(size=14, colour="black"), axis.text.y = element_text(size=14, colour = "black"), axis.title.y = element_text(size=18), axis.line=element_line(colour = "black"),plot.title = element_text(hjust = 0.5, size=18), legend.title = element_blank(), panel.background = element_blank(), legend.text = element_text(size=14), legend.spacing.x = unit(0.15, 'cm')) +
  labs(y="Counts", x="") +
  coord_cartesian(xlim=c(1,19.5))

# chrX females
empirical = read.csv("/agavescratch/tphung3/Turkana_Demography/05_generate_sfs/results/chrX/chrX.females.sfs.clean.txt", header = FALSE, sep="\t")
model_2epoch = read.csv("/agavescratch/tphung3/Turkana_Demography/06_demographic_inference/2Epoch_chrX_females_run_31_sfs_normalized_by_theta.txt", header = FALSE, sep="\t")

empirical_labels = rep("Empirical", 46)
twoepoch_labels = rep("2Epoch", 46)

dat = data.frame(Bins = empirical[,1], Counts=c(empirical[,2], model_2epoch[,2]), category = c(empirical_labels, twoepoch_labels))

dat$category = as.character(dat$category)
dat$category = factor(dat$category, levels = unique(dat$category))

ggplot(data=dat, aes(x=Bins, y=Counts, fill=category)) + 
  theme_bw() +
  geom_bar(stat = "identity", position=position_dodge(), colour="black") + 
  scale_fill_manual(values=c("gray40", "lightblue")) + 
  theme(axis.text.x = element_text(size=14, colour="black"), axis.text.y = element_text(size=14, colour = "black"), axis.title.y = element_text(size=18), axis.line=element_line(colour = "black"),plot.title = element_text(hjust = 0.5, size=18), legend.title = element_blank(), panel.background = element_blank(), legend.text = element_text(size=14), legend.spacing.x = unit(0.15, 'cm')) +
  labs(y="Counts", x="") +
  coord_cartesian(xlim=c(1,19.5))