
# XY biasd info
data1 <- read.table("mid_data/allele_specific/HJ/260124/XYbiased_info_geneUnit.v3.rmdup.table", header = T, sep = "\t")

head(data1)
data2 <- data1 %>% filter(tissue == "Mba_Mt")
table(data2$geneID %>% duplicated())
# testis exp ratio
ratio_df <- read.table("mid_data/allele_specific/HJ/250530/quantileScaled_expRatio.v2.table")
head(ratio_df)
my_cutoff <- quantile(ratio_df$Mbr_TR, probs = 0.9)


# combined XY biased info and testis exp ratio info
data3 <- data2 %>% left_join(y = ratio_df, by = "geneID") %>% na.omit()
head(data3)
table(is.na(data3$Mbr_KR))
a <- data3 %>% filter(adjusted_p_value < 0.05 & abs(logYXratio)>0.5 & Mbr_TR > my_cutoff)
# write.table(x = a, file = "mid_data/allele_specific/HJ/correctBiase_250427/testis_YXratio_testisExpRatio_selected_genes.txt", quote = F, row.names = F, col.names = T, sep = "\t")
a$type %>% table()
data3$type %>% table()

# 作图
library(ggrepel)



# 满足加边框条件的点
highlight_df <- subset(data3, 
                       Mbr_TR > my_cutoff & 
                         abs(logYXratio) > 0.5 & 
                         -log10(adjusted_p_value) > -log10(0.05))

# 不满足条件的点
background_df <- data3[!(rownames(data3) %in% rownames(highlight_df)), ]

data3$myType <- data3$type
data3$myType[data3$geneID %in% highlight_df$geneID] <- "testis and Y bias"
table(data3$myType)

p1 <- ggplot() +
  # 背景点（无边框）
  geom_point(data = data3, 
             aes(x = logYXratio, y = -log10(adjusted_p_value), 
                 size = Mbr_TR, color = myType), 
             alpha = 0.8) +
  # 辅助线
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#4A6FA5") +
  geom_vline(xintercept = c(0.5, -0.5), linetype = "dashed", color = "#4A6FA5") +
  scale_color_manual(values = c("Ybiased" = "#a1caf1",
                               "testis and Y bias" = "#0067a5",
                               "Xbiased" = "#be0032", 
                               "Not biased" = "#bebebe")) +
  scale_size_continuous(range = c(1, 5)) +
  
  # 文字标注
  geom_text_repel(data = highlight_df,
                  aes(x = logYXratio, y = -log10(adjusted_p_value), label = geneSymbol),
                  size = 3, max.overlaps = 50,
                  box.padding = 0.3, point.padding = 0.3,
                  segment.color = "grey50") +
  
  # 标签与主题
  labs(y = expression(-log[10](adjusted~italic(p))), 
       x = "log2(Y/X ratio)", 
       color = "Group", 
       fill = "Group", 
       size = "Mbr_TR") +
  ggthemes::theme_gdocs()
p1
ggsave(filename = "plots/DNA_pipe/hj/hifi/final2/allele_specific/260124/testisRatio_XYbiased.v3.pdf", plot = p1, width = 6, height = 4)
