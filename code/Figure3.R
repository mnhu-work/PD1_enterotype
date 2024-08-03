############Bacteria alpha diversity from microeco
alpha <- read.delim("../data/Bacteria_alpha_div_microeco_output.tsv",row.names=1)
dim(alpha)
meta <- read.csv("../data/clustering_results.csv")
identical(rownames(meta),rownames(alpha))
df <- cbind(meta,alpha)
source("/data/Humuni/source_functions/FN_colsets.R")
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(ggsignif) # Significance Brackets for 'ggplot2'
library(gghalves) 
library(ggpubr)
##
ggplot(df, aes(group, Shannon, fill = group)) +
  
  geom_half_violin(position = position_nudge(x = 0.3), side = "r", width = 0.8, color = NA, alpha = 0.6) +
  geom_jitter(aes(fill = group), shape = 21, size = 2.5, width = 0.2, alpha = 0.6) +
  geom_boxplot(width = 0.5, size = 0.75, outlier.color = NA, alpha = 0.6) +
  stat_compare_means(comparisons = list(c("favor_type", "unfavor_type")),
                     method = "wilcox.test", label = "p.format", size = 4) +
  scale_y_continuous() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 1),
        axis.text.x = element_text(color = "black", size = 13),
        axis.text.y = element_text(color = "black", size = 13),
        legend.position = "none",
        axis.ticks = element_line(color = "black", linewidth = 1)) +
  labs(x = "Group", y = "Alpha diversity (Shannon)") +
  scale_fill_manual(values = colset.d.2)
ggsave(file = "../../rarefying_596samples/Shannon_Bacter_microeco.pdf",height = 4,width = 2.5)


phylum_abu <- read.delim("../data/Bacteria_phylum_abu_microeco_output.tsv",row.names=1)
Phylum <- phylum_abu
#rowSums(Phylum[,c(1:7)])
mat <- Phylum[,c(1:40)]/100
mat <- as.data.frame(mat)
mat <- mat[, order(-colMeans(mat))]
mat$Others <-rowSums(mat[, 14:ncol(mat)])
mat <- mat[,-c(14:40)]
mat$Sample_ID <- rownames(mat)
meta$Sample_ID <- rownames(meta)
identical(mat$Sample_ID,meta$Sample_ID)
dat <- meta
#mat$group <- df$group
mat$Group <- dat$group
library(reshape2)
melted_df <- melt(mat, id.vars = "Group", variable.name = "phylum_name", value.name = "percentage")
melted_df$percentage <- as.numeric(melted_df$percentage)
sum(is.na(melted_df$percentage))
melted_df <- melted_df[!is.na(melted_df$percentage),]
library(dplyr)
results_df <- melted_df %>%
  group_by(Group, phylum_name) %>%
  summarize(average_percentage = mean(percentage))
############plot phylumn
###
col.hm <- c(colset.d.4,colset.d.12)
long_data <- results_df
names(long_data)
colnames(long_data)[2] <- "Phylumn"
library(ggplot2)
#dev.off()
g <- ggplot(long_data,aes(x=Group, y=average_percentage, fill=Phylumn))  +
  geom_bar(stat = "identity", width=0.8, col='black',alpha=0.4)  +
  #  theme_pander() +
  scale_fill_manual(values = col.hm) +
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x  = element_text(hjust=1, angle=45,size = 12),
        panel.background = element_rect(fill=NULL, colour = 'white')) 
facet_grid(. ~ Group)
g
#pdf(file = "../plots/clust_fungi_phylum_barplot_V2.pdf",width = 4,height = 6)
ggsave(filename = "../../rarefying_596samples/Bacter_Phylum_percentage_barplot.pdf",width = 4,height = 6)




identical(rownames(dat),rownames(phylum_abu))
data2 <- cbind(dat,Phylum)
names(Phylum)
table(data2$group)
###########################################################
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(ggsignif) # Significance Brackets for 'ggplot2'
library(gghalves)
library(ggpubr)
names(phylum_abu)
phylum_abu$group <- dat$group
F_B_ratio <- phylum_abu$F_B_ratio
phylum_abu$logF_B_ratio <- log2(F_B_ratio+1)
#pdf(file = "plots/Fungi_Basidiomycota.pdf",height = 4,width = 2.5)
ggplot(phylum_abu, aes(group, logF_B_ratio, fill = group)) +
  geom_half_violin(position = position_nudge(x = 0.25), side = "r", width = 0.8, color = NA, alpha = 0.6) +
  
  geom_jitter(aes(fill = group), shape = 21, size = 2.5, width = 0.2, alpha = 0.6) +
  geom_boxplot(width = 0.45, size = 0.75, outlier.color = NA, alpha = 0.6) +
  
  stat_compare_means(comparisons = list(c("favor_type", "unfavor_type")),
                     method = "wilcox.test", label = "p.format", size = 4) +
  scale_y_continuous() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 1),
        axis.text.x = element_text(color = "black", size = 8, angle = 45, hjust = 1),  # Rotate x-axis labels
        axis.text.y = element_text(color = "black", size = 8),
        legend.position = "none",
        axis.ticks = element_line(color = "black", linewidth = 1)) +
  labs(x = "Group", y = "logF_B_ratio") +
  scale_fill_manual(values = colset.d.2)
ggsave(filename = "../../rarefying_596samples/logF_B_ratio.pdf",height = 4,width = 2.5)
