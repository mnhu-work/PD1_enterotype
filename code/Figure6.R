####load data
meta.R <- read.csv("../data/Unique_RCP_meta_for_chisqtest.csv",row.names = 1)
table(meta.R$donor_origin_clust,meta.R$Clin_Response)
meta.R$donor_origin_clust <- paste0("Donor_origin_clust",meta.R$donor_origin_clust)
observed_dat<-data.frame(meta.R$donor_origin_clust,meta.R$Clin_Response)
colnames(observed_dat) <- c("Cluster","Clin_Response")
observed_dat$Clin_Response <- as.factor(observed_dat$Clin_Response)
class(observed_dat$Cluster)
chisq_result <- chisq.test(observed_dat$Clin_Response,observed_dat$Cluster)
observed_counts <- chisq_result$observed
observed_df <- as.data.frame(observed_counts)
observed_df$Category <- rownames(observed_counts)
observed_df$Group <- observed_df$observed_dat.Cluster

library(tidyverse)
observed_df <- observed_df %>%
  group_by(Group) %>%
  mutate(Percentage = Freq / sum(Freq))
#pdf(file="Cluster_chisq_FMT_output.pdf",height=4,width=6)
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
#ggsave(filename="../../results/Cluster_chisq_FMT_output.pdf",height=4,width=2.5)


###alpha diversity Fungi
alpha_fungi <- read.csv("../data/RCP_paired_alpha_index_fungi.csv",row.names=1)
package.list=c("tidyverse","ggsignif","ggsci","ggprism")

for (package in package.list) {
  if (!require(package,character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}
library(tidyverse)
library(ggpubr)

# Prepare the data frame
df <- alpha_fungi[, c("diversity_shannon", "patient_id", "S1_or_S4", "donor_origin_clust")]
df <- df[order(df$patient_id, df$S1_or_S4),]

# Create a new data frame with additional columns
df <- df %>% mutate(paired = rep(1:(n()/2), each=2), year = factor(S1_or_S4)) 
df1 <- df[, c("year", "diversity_shannon", "donor_origin_clust", "S1_or_S4", "paired")]

# Function to identify outliers
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

# Identify outliers
df1 <- df1 %>%
  group_by(donor_origin_clust) %>%
  mutate(outlier = ifelse(is_outlier(diversity_shannon), TRUE, FALSE)) %>%
  ungroup()

# Remove outliers and their pairs
df1_filtered <- df1 %>%
  filter(!paired %in% paired[outlier]) %>%
  select(-outlier)

# Print the structure of df1_filtered to diagnose the issue
print(str(df1_filtered))
print(summary(df1_filtered))

# Check if df1_filtered is empty and proceed only if it is not
if(nrow(df1_filtered) > 0) {
  # Calculate the y-axis limits based on the range of diversity_shannon
  y_min <- min(df1_filtered$diversity_shannon, na.rm = TRUE)
  y_max <- max(df1_filtered$diversity_shannon, na.rm = TRUE)
  
  # Plot the data
  p <- df1_filtered %>%
    ggplot(aes(year, diversity_shannon)) +
    stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.2), width = 0.1) +
    geom_boxplot(position = position_dodge(width = 0.2), width = 0.4, outlier.shape = NA) + # Hide outliers
    geom_line(aes(group = paired), position = position_dodge(0.2), color = "grey80") +
    geom_point(aes(fill = year, group = paired), pch = 21, position = position_dodge(0.2), size = 4, alpha = 0.6) + 
    scale_size_continuous(range = c(1, 3)) +
    facet_wrap(. ~ donor_origin_clust, nrow = 1) +
    scale_fill_npg() +
    labs(x = "FMT Timepoint", y = "Alpha diversity (Shannon)") +
    theme_prism(base_line_size = 0.5) +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(), # Remove major grid lines
      axis.text.y = element_text(color = "black", size = 10),
      axis.text.x = element_text(margin = margin(t = -5), color = "black", size = 10),
      legend.position = "none"
    ) +
    coord_cartesian(ylim = c(y_min, y_max)) + # Adjust y-axis limits
    stat_compare_means(
      method = "wilcox.test",
      label = "p.format",
      comparisons = list(c("S1", "S4")),
      size = 3,
      position = position_dodge(width = 0.2),
      label.y = y_max * 0.95,  # Adjust the y-coordinate for p-value labels
      vjust = 0.5,
      paired = TRUE  # Specify that the test is paired
    )
  
  print(p)
  
  # Save the plot
  ggsave(filename = "FMT_paired_fungi_alpha_diversity_shannon_wilcox_test.pdf", plot = p, width = 5, height = 4)
} else {
  print("The subset dataframe df1_filtered is empty. Please check the subsetting criteria.")
}



###alpha diversity Bacteria
alpha_fungi <- read.csv("../data/RCP_paired_alpha_indexes_bacteria.csv",row.names=1)
package.list=c("tidyverse","ggsignif","ggsci","ggprism")

for (package in package.list) {
  if (!require(package,character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}
library(tidyverse)
library(ggpubr)

library(tidyverse)
library(ggpubr)

library(tidyverse)
library(ggpubr)

# Prepare the data frame
df <- alpha_fungi[, c("diversity_shannon", "patient_id", "S1_or_S4", "donor_origin_clust")]
df <- df[order(df$patient_id, df$S1_or_S4),]

# Create a new data frame with additional columns
df <- df %>% mutate(paired = rep(1:(n()/2), each=2), year = factor(S1_or_S4)) 
df1 <- df[, c("year", "diversity_shannon", "donor_origin_clust", "S1_or_S4", "paired")]

# Function to identify outliers
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

# Identify outliers
df1 <- df1 %>%
  group_by(donor_origin_clust) %>%
  mutate(outlier = ifelse(is_outlier(diversity_shannon), TRUE, FALSE)) %>%
  ungroup()

# Remove outliers and their pairs
df1_filtered <- df1 %>%
  filter(!paired %in% paired[outlier]) %>%
  select(-outlier)

# Print the structure of df1_filtered to diagnose the issue
print(str(df1_filtered))
print(summary(df1_filtered))

# Check if df1_filtered is empty and proceed only if it is not
if(nrow(df1_filtered) > 0) {
  # Calculate the y-axis limits based on the range of diversity_shannon
  y_min <- min(df1_filtered$diversity_shannon, na.rm = TRUE)
  y_max <- max(df1_filtered$diversity_shannon, na.rm = TRUE)
  
  # Plot the data
  p <- df1_filtered %>%
    ggplot(aes(year, diversity_shannon)) +
    stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.2), width = 0.1) +
    geom_boxplot(position = position_dodge(width = 0.2), width = 0.4, outlier.shape = NA) + # Hide outliers
    geom_line(aes(group = paired), position = position_dodge(0.2), color = "grey80") +
    geom_point(aes(fill = year, group = paired), pch = 21, position = position_dodge(0.2), size = 4, alpha = 0.6) + 
    scale_size_continuous(range = c(1, 3)) +
    facet_wrap(. ~ donor_origin_clust, nrow = 1) +
    scale_fill_npg() +
    labs(x = "FMT Timepoint", y = "Alpha diversity (Shannon)") +
    theme_prism(base_line_size = 0.5) +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(), # Remove major grid lines
      axis.text.y = element_text(color = "black", size = 10),
      axis.text.x = element_text(margin = margin(t = -5), color = "black", size = 10),
      legend.position = "none"
    ) +
    coord_cartesian(ylim = c(y_min, y_max)) + # Adjust y-axis limits
    stat_compare_means(
      method = "wilcox.test",
      label = "p.format",
      comparisons = list(c("S1", "S4")),
      size = 3,
      position = position_dodge(width = 0.2),
      label.y = y_max * 0.95,  # Adjust the y-coordinate for p-value labels
      vjust = 0.5,
      paired = TRUE  # Specify that the test is paired
    )
  
  print(p)
  
  # Save the plot
  ggsave(filename = "FMT_paired_bacteria_alpha_diversity_shannon_wilcox_test.pdf", plot = p, width = 5, height = 4)
} else {
  print("The subset dataframe df1_filtered is empty. Please check the subsetting criteria.")
}