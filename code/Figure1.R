####load Figure1_data
data <- read.csv("../data/clustering_results.csv")

#### PCOA
#library(doBy)

library(tidyverse)
library(vegan)
library(RColorBrewer)
#' Tests
library(rstatix)
library(ggplot2)
library(patchwork)
library(dplyr)
mytheme <- theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = 'right',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
pca <- data
summary <- pca %>%
  group_by(group) %>%
  summarise(
    mean_pcoa1 = mean(PCoA1, na.rm = TRUE),
    mean_pcoa2 = mean(PCoA2, na.rm = TRUE)
  )
summary
pca_mean <- summary

load("../data/BC_distance_pcoa_results.RData")
source("../data/FN_colsets.R")
eig = summary(eigenvals(pcoa))
head(eig)
axis = paste0("PCoA", 1:ncol(eig))
axis
eig = data.frame(Axis = axis, t(eig)[, -3])
head(eig)
pca <- merge(pca, pca_mean, by = 'group')
head(pca)
pco1 = round(eig[1, 3] * 100, 2)
pco2 = round(eig[2, 3] * 100, 2)
xlab = paste0("PCoA1 (",pco1,"%)")
ylab = paste0("PCoA2 (",pco2,"%)")

pp <- ggplot(data = pca, aes(x = PCoA1,y = PCoA2)) +
  geom_point(size = 0.5,
             aes(col = group, shape = group)) +
  scale_shape_manual(values = c(16,15)) +
  scale_color_manual(values = colset.d.2) +
  labs(x = xlab, y = ylab) +
  mytheme
pp
summary <- pca %>%
  group_by(group) %>%
  summarise(
    mean_pcoa1 = mean(PCoA1, na.rm = TRUE),
    mean_pcoa2 = mean(PCoA2, na.rm = TRUE)
  )

# Visualization for radiating points
p9 <- pp +
  geom_segment(data = pca, aes(x = mean_pcoa1, y = mean_pcoa2,
                               xend = PCoA1, yend = PCoA2, color = group, alpha = 0.6),
               show.legend = FALSE) +
  geom_point(data = pca, aes(x = mean_pcoa1, y = mean_pcoa2, color = group),
             shape = 1, size = 3, stroke = 1.5, alpha = 0.8, show.legend = FALSE) +
  
  # Visualization for center points
  geom_point(data = pca %>% group_by(group) %>% summarize(center_x = mean(mean_pcoa1), center_y = mean(mean_pcoa2)),
             aes(x = center_x, y = center_y),
             shape = 1, size = 5, color = "black", fill = "transparent", alpha = 0.6, show.legend = FALSE)

p9
ggsave(filename="PcoA_bray_Fungi.pdf",height=4.5,width=6)



##########
observed_dat<-data.frame(data$group,data$ICB_response)
colnames(observed_dat) <- c("Cluster","Clin_Response")
observed_dat$Clin_Response <- as.factor(observed_dat$Clin_Response)
class(observed_dat$Cluster)
chisq_result <- chisq.test(observed_dat$Clin_Response,observed_dat$Cluster)
observed_counts <- chisq_result$observed
observed_df <- as.data.frame(observed_counts)
observed_df$Category <- rownames(observed_counts)
observed_df$Group <- observed_df$observed_dat.Cluster
observed_df <- observed_df %>%
  group_by(Group) %>%
  mutate(Percentage = Freq / sum(Freq))
gg <- ggplot(observed_df, aes(x = Group, y = Percentage, fill = Category)) +
  geom_col(alpha = 0.8, position = "fill", col = "black") + 
  labs(title = "",
       x = "Clusters",
       y = "Percentage") +
  theme_minimal() +
  scale_fill_manual(values = c("NR" = '#B2182B', "R" = "#2166AC")) +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 8), 
        plot.title = element_text(size = 8), 
        axis.title = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1))
p_value_label <- sprintf("p = %.3e", chisq_result$p.value)  
gg + annotate("text", x = 0.5, y = 0.5, label = p_value_label, size = 6, hjust = 0)

ggsave(filename = "clustering_ICB_chisq.pdf",height = 6,width = 2.5)





