library(tidyverse)
library(pafr)



test_alignment <- "mid_data/基因组组装/01.HJ/final/final2/synteny/BJ_HJr5md_20.paf"
ali <- read_paf(test_alignment)
ali <- ali %>% filter(nmatch > 300)
df_filtered <- ali %>% filter(qname == "1")
head(df_filtered)
table(df_filtered$tname)
my_color <- c("5" = "#e41a1c","12" = "#377eb8",
              "13" = "#fdb462","15" = "black",
              "1" = "grey","2"="grey","3"="grey","4"="grey","6"="grey","7"="grey","8"="grey",
              "9" = "grey","10"= "grey","11"="grey","14"="grey","19"="grey")
p1 <- ggplot()+  
  geom_rect(data=df_filtered, mapping=aes(xmin=qstart/1e6, xmax=qend/1e6, ymin=0, ymax=5, fill=tname), color=NA, alpha=0.8)+
  theme_classic()+
  scale_fill_manual(values = my_color)+
  theme(legend.position = "bottom") +
  labs(y = "Synteny with Dd")
p1



# male snp number
snpNum <- read.table("mid_data/DNA_pipe/HJ/hifi/final2/sex_snp/all_heterM_homoF_snp.100k.bed") %>% filter(V1=="1")
head(snpNum)
snpNum

p2 <- ggplot() +
  geom_line(snpNum, mapping = aes(x = V2/1e6, y = V4), color = "#008856") +
  theme_classic() +
  labs(y = "Male-specific heterozygote site number(per100kb)",x = NULL)
p2
# change point analysis
library(changepoint)
fit_changepoint = cpt.meanvar(snpNum$V4, method = "BinSeg", Q = 6)

plot(fit_changepoint)
# Return estimates
c(ints = param.est(fit_changepoint)$mean,
  cp = cpts(fit_changepoint))
# change point
my_cpt <- cpts(fit_changepoint)
my_cpt_point <- c(min(snpNum$V3)/1e6,(snpNum[my_cpt,] %>% pull(V2))/1e6, max(snpNum$V3)/1e6)
my_cpt_value <- param.est(fit_changepoint)$mean
p22 <- p2+ geom_rect(
  mapping = aes(xmin=my_cpt_point[-length(my_cpt_point)],
                xmax=my_cpt_point[-1], 
                ymin=my_cpt_value,
                ymax = my_cpt_value +2), color = "black", fill = "black")+
  geom_vline(xintercept = my_cpt_point[2], linetype = "dashed", color = "red")
p22
ggsave(filename = "plots/DNA_pipe/hj/hifi/final2/sex_chr/sigNum_changepoint.v3.pdf", width = 5, height = 1.5, plot = p22)
#significance test between strata
snpNum2 <- snpNum %>%
  mutate(
    midMb = (V2 + V3) / 2 / 1e6,
    type = case_when(
      midMb >= 0     & midMb < 8.9     ~ "PAR",
      midMb >= 8.9   & midMb < 71.4    ~ "S1",
      midMb >= 71.4  & midMb <= 98.9973 ~ "S0",
      TRUE ~ NA_character_
    )
  )
table(snpNum2$type, useNA = "ifany")
pairwise.wilcox.test(
  snpNum2$V4,
  snpNum2$type,
  p.adjust.method = "BH"
)
aggregate(V4 ~ type, snpNum2, median)

# FST plot
fst <- read.table("mid_data/DNA_pipe/HJ/hifi/final2//circos_file/fst.window.value") %>% filter(V1=="1")

head(fst)
fst$V4 <- as.numeric(fst$V4)
fst$V2 <- as.integer(fst$V2)
fst$V3 <- as.integer(fst$V3)
head(fst)
fst[fst$V4<0,"V4"] <- 0
p3 <- ggplot() +
  geom_point(fst, mapping = aes(x = V2/1e6, y = V4), color = "882d17") +
  theme_classic()+
  labs(x=NULL, y = "Fst")

p3
fit_changepoint2 = cpt.meanvar(fst$V4, method = "BinSeg", Q = 2)
# Return estimates
c(ints = param.est(fit_changepoint2)$mean,
  cp = cpts(fit_changepoint2))
my_cpt2 <- cpts(fit_changepoint2)
my_cpt_point2 <- c(min(fst$V3)/1e6,(fst[my_cpt2,] %>% pull(V2))/1e6, max(fst$V3)/1e6)
my_cpt_value2 <- param.est(fit_changepoint2)$mean
plot(fit_changepoint2)
p33 <- p3+ geom_rect(
  mapping = aes(xmin=my_cpt_point2[-length(my_cpt_point)],
                xmax=my_cpt_point2[-1], 
                ymin=my_cpt_value2,
                ymax = my_cpt_value2 +0.01), color = "black", fill = "black")
p33

fst2 <- fst %>%
  mutate(
    midMb = (V2 + V3) / 2 / 1e6,
    type = case_when(
      midMb >= 0     & midMb < 8.45     ~ "PAR",
      midMb >= 8.45   & midMb < 19.4    ~ "R1",
      midMb >= 19.4  & midMb <= 99 ~ "R2",
      TRUE ~ NA_character_
    )
  )
table(fst2$type, useNA = "ifany")
pairwise.wilcox.test(
  fst2$V4,
  fst2$type,
  p.adjust.method = "BH"
)
aggregate(V4 ~ type, fst2, median)



# ---------------
# GC and repeat
GC <- read.table("mid_data/DNA_pipe/HJ/hifi/final2/circos_file/HJ.GCcontent.100k.bed") %>% filter(V1=="1")
head(GC)
p4 <- ggplot() +
  geom_line(GC, mapping = aes(x = V2/1e6, y = V4), color = "#0067a5") + 
  theme_classic()+
  labs(x=NULL, y = "GC%")
p4

head(GC)
GC$Type <- ifelse(test = GC$V2 < 8900000, yes = "PAR", no = "SDR")
ggplot(GC, mapping = aes(x = Type, y = V4)) +
  geom_boxplot()
ggplot(GC, aes(x = Type, y = V4)) +
  geom_boxplot() +
  ggpubr::stat_compare_means(
    method = "wilcox.test",
    label = "p.format"
  )
ggplot(GC, aes(x = Type, y = V4)) +
  geom_boxplot() +
  ggpubr::stat_compare_means(
    method = "wilcox.test",
    label = "p.signif"
  )

rep_data <- read.table("mid_data/DNA_pipe/HJ/hifi/final2/circos_file/HJ.repeatDensity.100k.bed") %>% filter(V1=="1")
p5 <- ggplot() +
  geom_line(rep_data, mapping = aes(x = V2/1e6, y = V4), color = "#dcd300")+
  theme_classic()+
  labs(x="Chr1(Mb)", y = "Repeat%")

p5

centromere_data <- read.table("mid_data/DNA_pipe/HJ/hifi/final2/circos_file/chr16_TR00407.blastOut.100k.bed") %>% filter(V1=="1")
p6 <- ggplot() +
  geom_line(centromere_data, mapping = aes(x = V2/1e6, y = V4), color = "#be0032")+
  theme_classic()+
  labs(x=NULL, y = "Repeat%")

p6

library(grid)
pdf(file = "plots/DNA_pipe/hj/hifi/final2/sex_chr/combined.v2.pdf",width = 5, height = 8)
grid.newpage()
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p5),ggplotGrob(p22),ggplotGrob(p33),ggplotGrob(p4) ,ggplotGrob(p6),size = "last"))
dev.off()



# male testis biased snp allele
X_biasedSNP <- read.table("mid_data/allele_specific/HJ/ASEsignif_Xbiased_Mtestis.table", header = T, sep = "\t")
head(X_biasedSNP)
X_biasedSNP$type3 <- ifelse(X_biasedSNP$Xcount > X_biasedSNP$Ycount, yes = "Xbiased", no = "Ybiased")
X_biasedSNP$type3[X_biasedSNP$adjusted_p_value > 0.05] <- "Not signif"
p7 <- ggplot(X_biasedSNP, mapping = aes(x = position, y = -log(adjusted_p_value+0.001), colour = type3))+
  geom_point()+
  geom_hline(yintercept = -log(0.05), linetype = "dashed")+
  scale_color_manual(values = c("Xbiased" = "#be0032","Ybiased"="#0067a5", "Not signif" = "#bebebe"))+
  theme_classic()+
  theme(legend.position = "none")
p7
ggsave(filename = "plots/DNA_pipe/hj/hifi/final2/sex_chr/XYbiasedSNP_plot.pdf", width = 5, height = 1.5, plot = p7)

p8 <- ggplot(X_biasedSNP, mapping = aes(x = position, colour = type3))+
  geom_density(linewidth =2) +
  scale_color_manual(values = c("Xbiased" = "#be0032","Ybiased"="#0067a5", "Not signif" = "#bebebe"))+
  theme_classic()+
  theme(legend.position = "none")
p8
