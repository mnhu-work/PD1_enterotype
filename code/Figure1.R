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
#自定义主题和配色：
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

#设置各维度的名字，从PCoA1开始
axis = paste0("PCoA", 1:ncol(eig))
axis
#各轴解释率
eig = data.frame(Axis = axis, t(eig)[, -3])
head(eig)
pca <- merge(pca, pca_mean, by = 'group')
head(pca)
pco1 = round(eig[1, 3] * 100, 2)
pco2 = round(eig[2, 3] * 100, 2)
# 设置画图时的x轴和y轴的标题，衔接逗号中的内容
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
# 提取卡方检验结果中的观察频数
observed_counts <- chisq_result$observed
# 将观察频数数据转换为数据框
observed_df <- as.data.frame(observed_counts)
observed_df$Category <- rownames(observed_counts)
observed_df$Group <- observed_df$observed_dat.Cluster
# 计算百分比
observed_df <- observed_df %>%
  group_by(Group) %>%
  mutate(Percentage = Freq / sum(Freq))
gg <- ggplot(observed_df, aes(x = Group, y = Percentage, fill = Category)) +
  geom_col(alpha = 0.8, position = "fill", col = "black") +  # 使用geom_col绘制竖直柱状图，并设置颜色透明度为0.4
  labs(title = "",
       x = "Clusters",
       y = "Percentage") +
  theme_minimal() +
  scale_fill_manual(values = c("NR" = '#B2182B', "R" = "#2166AC")) +
  theme(panel.background = element_blank(),  # 去除背景
        panel.grid.major = element_blank(),  # 去除主要网格线
        panel.grid.minor = element_blank(),  # 去除次要网格线
        axis.text = element_text(size = 8),  # 调整坐标轴文字字体大小
        plot.title = element_text(size = 8),  # 调整图标题字体大小
        axis.title = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1))
p_value_label <- sprintf("p = %.3e", chisq_result$p.value)  # 将p值格式化为科学计数法，保留三位有效数字
gg + annotate("text", x = 0.5, y = 0.5, label = p_value_label, size = 6, hjust = 0)

ggsave(filename = "../rarefying_596samples/pam_clustering_ICB_chisq.pdf",height = 6,width = 2.5)





