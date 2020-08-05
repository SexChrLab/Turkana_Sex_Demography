library(ggplot2)
library(ggpubr)

autosomes_stats = read.csv('/scratch/tphung3/Turkana_Demography/03_data_cleaning/hard_filter_variants/stats/autosomes_females_males_annotations.txt', sep = '\t')
chrX_females_stats = read.csv('/scratch/tphung3/Turkana_Demography/03_data_cleaning/hard_filter_variants/stats/chrX_females_annotations.txt', sep = '\t')
chrX_males_stats = read.csv('/scratch/tphung3/Turkana_Demography/03_data_cleaning/hard_filter_variants/stats/chrX_males_annotations.txt', sep = '\t')

# ---
# QD
# ---
png('/scratch/tphung3/Turkana_Demography/03_data_cleaning/hard_filter_variants/stats/QD_autosomes_chrX.png', width=13, height=5, units = "in", res = 300)
p1 = ggplot(autosomes_stats, aes(x=QD)) +
  geom_density() +
  theme_bw() +
  theme(axis.text.x = element_text(size=16, colour="black"), axis.text.y = element_text(size=16, colour = "black"), axis.title= element_text(size=16), axis.line=element_line(colour = "black"), axis.ticks = element_line(colour = "black", size=1.5), axis.ticks.length = unit(0.3, 'cm'), plot.title = element_text(hjust = 0.5, size=18)) +
  geom_vline(xintercept = 2, linetype=2, size=1.5, color="darkgray") +
  labs(title="Autosomes \n Females & Males") + 
  coord_cartesian(ylim=c(0, 0.1))

p2 = ggplot(chrX_females_stats, aes(x=QD)) +
  geom_density() +
  theme_bw() +
  theme(axis.text.x = element_text(size=16, colour="black"), axis.text.y = element_text(size=16, colour = "black"), axis.title= element_text(size=16), axis.line=element_line(colour = "black"), axis.ticks = element_line(colour = "black", size=1.5), axis.ticks.length = unit(0.3, 'cm'), plot.title = element_text(hjust = 0.5, size=18)) +
  geom_vline(xintercept = 2, linetype=2, size=1.5, color="darkgray") +
  labs(title="X chromosome \n Females") + 
  coord_cartesian(ylim=c(0, 0.1))

p3 = ggplot(chrX_males_stats, aes(x=QD)) +
  geom_density() +
  theme_bw() +
  theme(axis.text.x = element_text(size=16, colour="black"), axis.text.y = element_text(size=16, colour = "black"), axis.title= element_text(size=16), axis.line=element_line(colour = "black"), axis.ticks = element_line(colour = "black", size=1.5), axis.ticks.length = unit(0.3, 'cm'), plot.title = element_text(hjust = 0.5, size=18)) +
  geom_vline(xintercept = 2, linetype=2, size=1.5, color="darkgray") +
  labs(title="X chromosome \n Males") + 
  coord_cartesian(ylim=c(0, 0.1))

ggarrange(p1, p2, p3, ncol = 3, nrow = 1)
dev.off()

# ---
# DP
# ---
png('/scratch/tphung3/Turkana_Demography/03_data_cleaning/hard_filter_variants/stats/DP_autosomes_chrX.png', width=13, height=5, units = "in", res = 300)
p1 = ggplot(autosomes_stats, aes(x=DP)) +
  geom_density() +
  theme_bw() +
  theme(axis.text.x = element_text(size=16, colour="black"), axis.text.y = element_text(size=16, colour = "black"), axis.title= element_text(size=16), axis.line=element_line(colour = "black"), axis.ticks = element_line(colour = "black", size=1.5), axis.ticks.length = unit(0.3, 'cm'), plot.title = element_text(hjust = 0.5, size=18)) +
  labs(title="Autosomes \n Females & Males") + 
  coord_cartesian(ylim=c(0, 0.05))

p2 = ggplot(chrX_females_stats, aes(x=DP)) +
  geom_density() +
  theme_bw() +
  theme(axis.text.x = element_text(size=16, colour="black"), axis.text.y = element_text(size=16, colour = "black"), axis.title= element_text(size=16), axis.line=element_line(colour = "black"), axis.ticks = element_line(colour = "black", size=1.5), axis.ticks.length = unit(0.3, 'cm'), plot.title = element_text(hjust = 0.5, size=18)) +
  labs(title="X chromosome \n Females") + 
  coord_cartesian(ylim=c(0, 0.05))

p3 = ggplot(chrX_males_stats, aes(x=DP)) +
  geom_density() +
  theme_bw() +
  theme(axis.text.x = element_text(size=16, colour="black"), axis.text.y = element_text(size=16, colour = "black"), axis.title= element_text(size=16), axis.line=element_line(colour = "black"), axis.ticks = element_line(colour = "black", size=1.5), axis.ticks.length = unit(0.3, 'cm'), plot.title = element_text(hjust = 0.5, size=18)) +
  labs(title="X chromosome \n Males") + 
  coord_cartesian(ylim=c(0, 0.05))

ggarrange(p1, p2, p3, ncol = 3, nrow = 1)
dev.off()

# -----------
# Extra stats
# -----------
autosomes_stats = read.csv('/scratch/tphung3/Turkana_Demography/03_data_cleaning/hard_filter_variants/stats/autosomes_females_males_annotations_extra.txt', sep = '\t')
chrX_females_stats = read.csv('/scratch/tphung3/Turkana_Demography/03_data_cleaning/hard_filter_variants/stats/chrX_females_annotations_extra.txt', sep = '\t')
chrX_males_stats = read.csv('/scratch/tphung3/Turkana_Demography/03_data_cleaning/hard_filter_variants/stats/chrX_males_annotations_extra.txt', sep = '\t')

# ---
# FS
# ---
png('/scratch/tphung3/Turkana_Demography/03_data_cleaning/hard_filter_variants/stats/FS_autosomes_chrX.png', width=13, height=5, units = "in", res = 300)
p1 = ggplot(autosomes_stats, aes(x=FS)) +
  geom_density() +
  theme_bw() +
  theme(axis.text.x = element_text(size=16, colour="black"), axis.text.y = element_text(size=16, colour = "black"), axis.title= element_text(size=16), axis.line=element_line(colour = "black"), axis.ticks = element_line(colour = "black", size=1.5), axis.ticks.length = unit(0.3, 'cm'), plot.title = element_text(hjust = 0.5, size=18)) +
  labs(title="Autosomes \n Females & Males") +
  coord_cartesian(ylim=c(0, 1.2))

p2 = ggplot(chrX_females_stats, aes(x=FS)) +
  geom_density() +
  theme_bw() +
  theme(axis.text.x = element_text(size=16, colour="black"), axis.text.y = element_text(size=16, colour = "black"), axis.title= element_text(size=16), axis.line=element_line(colour = "black"), axis.ticks = element_line(colour = "black", size=1.5), axis.ticks.length = unit(0.3, 'cm'), plot.title = element_text(hjust = 0.5, size=18)) +
  labs(title="X chromosome \n Females") +
  coord_cartesian(ylim=c(0, 1.2))

p3 = ggplot(chrX_males_stats, aes(x=FS)) +
  geom_density() +
  theme_bw() +
  theme(axis.text.x = element_text(size=16, colour="black"), axis.text.y = element_text(size=16, colour = "black"), axis.title= element_text(size=16), axis.line=element_line(colour = "black"), axis.ticks = element_line(colour = "black", size=1.5), axis.ticks.length = unit(0.3, 'cm'), plot.title = element_text(hjust = 0.5, size=18)) +
  labs(title="X chromosome \n Males") +
  coord_cartesian(ylim=c(0, 1.2))

ggarrange(p1, p2, p3, ncol = 3, nrow = 1)
dev.off()
