
# 设置环境
rm(list = ls()) 
options(stringsAsFactors = F) 

library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(stringr)
library(data.table)
library(pheatmap)
library(ggheatmap)
library(RColorBrewer)
library(ggthemes)
library(patchwork)
library(scRNAtoolVis)
library(ggpubr)
library(ggsignif)
library(rstatix)

# 设置数据
load(file = "./4-celltype/3-sce.all.filt.int.celltype.Rdata")
load(file = "./1-datapepare/mycol.Rdata")
table(sce.all.filt.int$celltype)

# 细胞数目
meta = sce.all.filt.int@meta.data
meta = meta[!duplicated(meta$orig.ident),]
hcid = meta$orig.ident[meta$diagnosis == "HC"]


cellcount = as.data.frame(prop.table(table(sce.all.filt.int$celltype, sce.all.filt.int$orig.ident), 2))
colnames(cellcount) = c("Celltype", "Patient", "Proportion")
cellcount$group = ifelse(cellcount$Patient %in% hcid, "HC", "AD")

df1 = subset(cellcount, group == "HC")
df2 = subset(cellcount, group == "AD")


p1 = ggplot(df1, aes(Patient, Proportion, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_classic() +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = col_vec)+
  geom_bar(stat = "identity") +
  ggtitle("HC")+
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

p2 = ggplot(df2, aes(Patient, Proportion, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_classic() +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = col_vec)+
  geom_bar(stat = "identity") +
  ggtitle("AD")+
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())



p3all = p1 + p2 + plot_layout(guides = "collect", widths = c(3,4))
ggsave(p3all, filename = "./4-celltype/5-patient-proportion.pdf", width = 12, height = 2.5)



cellcount = as.data.frame(prop.table(table(sce.all.filt.int$celltype, sce.all.filt.int$diagnosis), 2))
colnames(cellcount) = c("Celltype", "Group", "Proportion")

p1 = ggplot(cellcount, aes(Group, Proportion, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_classic() +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = col_vec)+
  geom_bar(stat = "identity") +
  ggtitle("HC")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave(filename = "./3-celltype/6-bar-group.pdf", width = 13, height = 8)



# 细胞比例
id1 = meta$orig.ident[meta$gse_id == "GSE157827"]
id2 = meta$orig.ident[meta$gse_id == "GSE167494"]
id3 = meta$orig.ident[meta$gse_id == "GSE174367"]

cellcount = as.data.frame(prop.table(table(sce.all.filt.int$celltype, sce.all.filt.int$orig.ident), 2))
colnames(cellcount) = c("Celltype", "Patient", "Proportion")
cellcount$group = ifelse(cellcount$Patient %in% hcid, "HC", "AD")
cellcount$group = factor(cellcount$group, levels = c("HC", "AD"))
cellcount$gse = ifelse(cellcount$Patient %in% id1, "GSE157827",
                       ifelse(cellcount$Patient %in% id2, "GSE167494", "GSE174367"))

p1 = ggplot(cellcount, aes(x = Celltype, y = Proportion, fill = group)) +
  stat_summary(fun = mean, geom = 'col', 
               position = position_dodge2(),
               width = 0.9) +
  geom_point(aes(color = factor(group)),  ## 添加抖动点
             show.legend = F, ## 不显示图例
             position = position_jitterdodge(),
             shape = 21,
             size = 1.5,
             fill = 'white')+
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', ## 添加误差线
               color = 'black',
               position = position_dodge(width = 0.9),
               width = 0.3, size = 0.4) +
  theme_classic()+
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  scale_fill_manual(values = col_vec)+
  scale_color_manual(values = col_vec)+
  theme(axis.text.x = element_text(colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black"))

p1
ggsave(p1, filename = "./4-celltype/6-dot-group.pdf", width = 6, height = 3.5)

stat_test <- cellcount %>%
  group_by(Celltype) %>%
  t_test(Proportion ~ group) %>%
  adjust_pvalue(method = "BH") %>% 
  add_significance("p") 

write.csv(stat_test, file = "./3-celltype/Celltype_signaficant.csv")

