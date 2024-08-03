###########KO
Mebo_KO <- read.csv("../data/KO_differential_results.csv",row.names = 1)

##Sankey plot

library(ggplot2)
library(ggalluvial)
library(networkD3)
library(tidyverse)
#rm(list = ls())
plot_dat <- subset(Mebo_KO,abs(Mebo_KO$cohen.d) >= 0.65 & Mebo_KO$p.adj < 10e-05)
names(plot_dat)
library(ggalluvial)
df <- plot_dat[,c("Level2","Level3","KO_id","cohen.d","p.adj")]
df <- df[order(df$cohen.d), ] 
df$Level3 <- factor(df$Level3,levels = df$Level3 %>% unique())
df$KO_id <- factor(df$KO_id ,levels = df$KO_id  %>% unique())
df$Level2 <- factor(df$Level2 ,levels = df$Level2  %>% unique())
#df$id<-factor(df$id,levels = df$id)
df$Level3 <- gsub("\\[.*", "", df$Level3)
library(dplyr)
table(df$Level2)
names(df)
df1 <- df[,c(1:3)]
#df1 <- df1[order(df1$KO_id),]
#df1$Level2 <- as.factor(df1$Level2)
#df1$Level3 <- as.factor(df1$Level3)
#df1$KO_id <- as.factor(df1$KO_id)
df1$KO_id <- factor(df1$KO_id, levels = unique(df1$KO_id))  # Convert Cohort to factor with unique levels
# Sort the data frame by Cohort in descending order
df1 <- arrange(df1, desc(KO_id))
df1$KO_id

#df1$KO_id <- factor(df1$KO_id,levels = rev(df2$KO_id) %>% unique())
UCB_lodes <- to_lodes_form(df1[,1:ncol(df1)],
                           axes = 1:ncol(df1),
                           id = "Cohort")
dim(UCB_lodes)
head(UCB_lodes)
tail(UCB_lodes)
pdf(file="outFile_V2.pdf",width=7,height=6)
source("FN_colsets.R")
mycol <- rep(c( "#FFD1DC","#6699CC", "#F0F0F0", "#CCE1F2", "#E6E6FA",  "#99C2EB", "#FFFACD", "#F7C6C6", "#FFF8DC"), 15)
print(mycol)
#"#A9A9A9", "#FFD700", "#FF69B4", "#ADFF2F", "#B0C4DE", "#FF4500", "#FFFFE0", "#B22222"
ggplot(UCB_lodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
  scale_x_discrete(expand = c(0, 0)) +  
  geom_flow(width = 2/10,aes.flow = "forward") + 
  geom_stratum(alpha = .9,width = 2/10) +
  scale_fill_manual(values = mycol) +
  geom_text(stat = "stratum", size = 2,color="black") +
  xlab("") + ylab("") + theme_bw() + 
  theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank()) + #ȥ????????
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_blank()) + 
  ggtitle("") + guides(fill = FALSE) 
dev.off()
######################################################

# buuble
names(df)

df2 <- df
df2$KO_id <- factor(df2$KO_id,levels = df2$KO_id %>% unique())
df2 <- arrange(df2, desc(KO_id))
pdf(file="outFile_buble_V2.pdf",width=4,height=6)
pl<-ggplot(df2,aes(cohen.d,KO_id,colour=-log10(p.adj)))+   
  geom_point(size=3,alpha=0.8)+
  #scale_color_gradientn(colours=c("#CBD588","#599861"))+
  scale_color_gradientn(colours=c("#CCE1F2","#2166AC"))+
  #scale_fill_brewer(palette = "Blues",direction = -1)+
  #scale_size_continuous(range = c(5,6))
  #scale_x_continuous(limits = c(1.5,5))+
  theme_bw(15)+
  xlab("Cohen_Distance")+
  ylab("")+
  theme(legend.position=c(1,0),legend.justification = c(1,0))+
  theme(legend.background = element_blank())+
  theme(legend.key = element_blank()) 

#coord_flip()+
#theme(axis.text.x = element_text(angle = 45, hjust = 1))
pl
dev.off()


#################MetaCyc
ID <- c("CENTFERM-PWY: pyruvate fermentation to butanoate","PWY-5676: acetyl-CoA fermentation to butanoate II",
        "PWY-6590: superpathway of Clostridium acetobutylicum acidogenic fermentation","GALACTUROCAT-PWY: D-galacturonate degradation I",
        "PWY66-422: D-galactose degradation V (Leloir pathway)","PWY-6737: starch degradation V")

results <- read.csv("../data/MetaCyc_results",row.names = 1)
##################################
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
#rm(list = ls())
df <- subset(results,results$p.adj!=0 & results$path_ID!="UNMAPPED") 
## 添加一列logCPM
df$logFDR = -log10(df$p.adj)
#确定是上调还是下调，用于给图中点上色
df$threshold = factor(ifelse(df$p.adj  < 0.05 & abs(df$fc) >0, ifelse(df$fc >= 0 ,'Enriched in R_type','Enriched in NR_type'),'NoSignifi'),levels=c('Enriched in R_type','Enriched in NR_type','NoSignifi'))

table(df$threshold)
df <- as.data.frame(df)
ggplot(df,aes(x=fc,y= -log10(p.adj),size = logFDR,fill = threshold,alpha=0.6))+
  geom_point(colour = "black", shape = 21, stroke = 0.5)+
  #scale_size(limits  = c(2, 30))+ #
  scale_fill_manual(values=c("#2166AC","#B2182B","#bdbdbd"))+
  geom_text_repel(
    data = df[df$p.adj<0.05&abs(df$cohen.d)>0.75,],
    aes(label = path_ID),
    size = 4.5,
    color = "black",
    segment.color = "black", show.legend = FALSE )+
  ylab('-log10 (Pvalue)')+
  xlab('Generalized Fold Change')+
  geom_vline(xintercept=c(-0.1,0.1),lty=2,col="black",lwd=0.5)
  geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.5)
  theme_classic(
    base_line_size = 1
  )+
  guides(fill=guide_legend(override.aes = list(size=5)))+
  theme(axis.title.x = element_text(size = 10, 
                                    color = "black",
                                    face = "bold"),
        axis.title.y = element_text(size = 10,
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 7, 
                                   face = "bold")
  )



ggsave("metacyc_path_volcano.pdf", height = 5, width = 6)





###### Yachidas metabolism validation
dat <- read.csv("../data/rf_pred_Yachidas_results.csv",row.names = 1)
names(dat)
table(dat$rf_pred)
source("../data/FN_colsets.R")
names(dat)
dat$group <- ifelse(dat$rf_pred=="1","favor-type","unfavor-type")
ggplot(dat, aes(group, log2Butyrate_abu, fill = group)) +
  
  geom_violin(width = 1.0, color = NA, alpha = 0.6) +
  geom_jitter(aes(fill = group), shape = 21, size = 2.5, width = 0.2, alpha = 0.6) +
  geom_boxplot(width = 0.6, size = 0.75, outlier.color = NA, alpha = 0.6) +
  stat_compare_means(comparisons = list(c("favor-type", "unfavor-type")),
                     method = "wilcox.test", label = "p.format", size = 4) +
  scale_y_continuous() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        #panel.border = element_rect(size = 1),
        axis.text.x = element_text(color = "black", size = 6),
        axis.text.y = element_text(color = "black", size = 13),
        legend.position = "none",
        axis.ticks = element_line(color = "black", linewidth = 0.5)) +
  labs(x = "Group", y = "Normalized abudance of Butyrate") +
  scale_fill_manual(values = colset.d.2)
ggsave(filename = "Yachidas_Butyrate_rf_pred_compr.pdf",height = 4,width = 2.5)


#############
#Faecalibacterium prausnitzii/Faecalibacterium duncaniae/Lachnospira eligens
ggscatter(data, x = "Lachnospira_eligens", y = "log2Butyrate_abu",
          color = "black", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "#92c5de", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          xlab=c('Relative Abundance of Lachnospira eligens'),
          ylab='Normalized abundance of C00246_Butyrate',
          cor.coeff.args = list(method = "spearman", label.sep = "\n")
)
ggsave(filename = "Yachidas_Lachnospira eligens_abucor.pdf",height = 4,width = 3.5)