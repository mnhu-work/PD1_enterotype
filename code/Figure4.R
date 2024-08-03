######### Figrue 4

SparCC_favor <- read.csv("../data/SparCC_favor_correlations.csv",row.names = 1)
SparCC_unfavor <- read.csv("../data/SparCC_unfavor_correlations.csv",row.names = 1)
table(SparCC_favor$new_column)
table(SparCC_unfavor$new_column)


data <- data.frame(
  new_column = c("Bacter_Bacter", "Bacter_fungi", "fungi_fungi"),
  favor_type = c(382, 757, 339),
  unfavor_type = c(387, 821, 375)
)

# 加载必要的包
library(tidyr)
library(dplyr)
library(ggplot2)

# 原始数据
data <- data.frame(
  new_column = c("Bacter_Bacter", "Bacter_fungi", "fungi_fungi"),
  favor_type = c(382, 757, 339),
  unfavor_type = c(387, 821, 375)
)

# 数据变形为长格式
data_long <- data %>%
  pivot_longer(cols = c(favor_type, unfavor_type), names_to = "category", values_to = "count")

# 计算百分比
data_long <- data_long %>%
  group_by(category) %>%
  mutate(percent = count / sum(count) * 100) %>%
  ungroup()

# 绘制圆饼图
ggplot(data_long, aes(x="", y=percent, fill=new_column)) + 
  geom_bar(stat="identity", width=1, alpha=0.6) +
  coord_polar("y") +
  facet_wrap(~category, scales="free") +
  scale_fill_manual(values= colset.d.12) +
  geom_text(aes(label=sprintf("%.1f%%", percent)), 
            position=position_stack(vjust=0.5), size=3) +
  labs(title="Comparison of Groups in Favor and Unfavor",
       fill="Group") +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank())

ggsave(filename="SparCC_DE_correlations_compar.pdf",width=5,height=3.5)


favor_fungi_bacter <- subset(unique_favor,unique_favor$new_column=="Bacter_fungi")
unfavor_fungi_bacter <- subset(unique_unfavor,unique_unfavor$new_column=="Bacter_fungi")
table(favor_fungi_bacter$direction)
table(unfavor_fungi_bacter$direction)

data1 <- data.frame(
  new_column = c("positive", "negtive"),
  favor_type = c(148,609),
  unfavor_type = c(43,778)
)

# 加载必要的包
library(tidyr)
library(dplyr)
library(ggplot2)

# 原始数据
data1 <- data.frame(
  new_column = c("positive", "negative"),
  favor_type = c(148, 609),
  unfavor_type = c(43, 778)
)

# 数据变形为长格式
data1_long <- data1 %>%
  pivot_longer(cols = c(favor_type, unfavor_type), names_to = "category", values_to = "count")

# 计算百分比
data1_long <- data1_long %>%
  group_by(category) %>%
  mutate(percent = count / sum(count) * 100) %>%
  ungroup()

# 自定义颜色
my_colors <- c("positive" = "#b15928", "negative" = "#b2df8a")

# 绘制堆叠百分比柱形图
p <- ggplot(data1_long, aes(fill = new_column, y = percent, x = category)) + 
  geom_bar(position = "stack", stat = "identity", alpha = 0.6) +
  scale_fill_manual(values = my_colors) +
  labs(title = "Comparison of Groups in Favor and Unfavor with negative correlations",
       x = "Category",
       y = "Percentage",
       fill = "Group") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# 显示图表
print(p)

# 保存图表为 PDF 文件
ggsave(filename="Bacter_fungi_correlation_compar.pdf", plot = p, width=4.5, height=6)


contingency_table <- data1 %>%
  select(-new_column) %>%
  as.matrix()
chisq.test(contingency_table)