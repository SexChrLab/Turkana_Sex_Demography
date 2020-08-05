library(ggplot2)
library(ggpubr)

chr8_labels = rep("chr8", 18)
female_labels = rep("chrX Females", 18)
male_labels = rep("chrX Males", 18)

dat = data.frame(Bins = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18), Counts=c(21446.69727614205, 9349.686898578557, 5969.48493816625, 4439.873384940814, 3569.451327457855, 3009.17179504345, 2613.271066995873, 2313.468072772218, 2077.685677808519, 1889.867532698215, 1740.471754197876, 1622.647233457915, 1530.865075212979, 1460.674089417163, 1408.831525414207, 1373.207706455185, 1352.421601312946, 672.7995862642485, 15651.5517327103, 6533.260849007543, 4000.922433337113, 2951.149067875618, 2392.32366773181, 2040.31124206564, 1791.4062921717, 1606.302125757627, 1465.20311011423, 1355.647476072807, 1268.557423511563, 1196.552866391428, 1135.501093976001, 1085.642650679327, 1049.310105122403, 1027.183559119234, 1016.679221436713, 506.9062637696446, 13581, 5562, 3451, 2608, 2021, 1875, 1676, 1447, 1346, 1325, 1242, 1180, 1046, 1082, 974, 915, 960, 446), category = c(chr8_labels, female_labels, male_labels))

dat$category = as.character(dat$category)
dat$category = factor(dat$category, levels = unique(dat$category))
png("/scratch/tphung3/Turkana_Demography/05_generate_sfs/plots/compare_empirical_chr8_chrX.females.males.png", width = 9, height = 5, units = 'in', res = 300)
ggplot(data=dat, aes(x=Bins, y=Counts, fill=category)) + 
  theme_bw() +
  geom_bar(stat = "identity", position=position_dodge(), colour="black") + 
  scale_fill_manual(values=c("gray40", "lightblue", "gold")) + 
  theme(axis.text.x = element_text(size=14, colour="black"), axis.text.y = element_text(size=14, colour = "black"), axis.title.y = element_text(size=18), axis.line=element_line(colour = "black"),plot.title = element_text(hjust = 0.5, size=18), legend.title = element_blank(), panel.background = element_blank(), legend.text = element_text(size=14), legend.spacing.x = unit(0.15, 'cm')) +
  labs(y="Counts", x="")
dev.off()