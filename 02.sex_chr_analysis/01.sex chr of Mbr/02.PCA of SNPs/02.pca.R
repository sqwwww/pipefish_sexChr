library(tidyverse)

a <- read.table("mid_data/DNA_pipe/HJ/hifi/final2/pca/chr1_heterM_homoF/snp.eigenvec")
b <- read.table("mid_data/DNA_pipe/HJ/hifi/final2/pca/chr1_heterM_homoF/snp.eigenval")
head(a)
head(b)
pc1 <- round(b$V1[1] / sum(b$V1) *100)
pc2 <- round(b$V1[2] / sum(b$V1) *100)
a$type <- a$V1 %>% str_sub(start = 1, end = 1)
table(a$type)
p1 <- ggplot(data = a, mapping = aes(x = V3, y = V4, colour = type))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c("F" = "#be0032", "M" = "#0067a5") ) +
  labs(x= paste0("PC1(explained variance, ",pc1,"%)"), y = paste0("PC2(explained variance, ",pc2,"%)"))
p1
ggsave(filename = "plots/DNA_pipe/hj/hifi/final2/pca/chr1_heterM_homoF_SNP.pdf", width = 4, height = 3, plot = p1)


a <- read.table("mid_data/DNA_pipe/HJ/hifi/final2/pca/all_snp/snp.eigenvec")
b <- read.table("mid_data/DNA_pipe/HJ/hifi/final2/pca/all_snp/snp.eigenval")
head(a)
head(b)
pc1 <- round(b$V1[1] / sum(b$V1) *100)
pc2 <- round(b$V1[2] / sum(b$V1) *100)
a$type <- a$V1 %>% str_sub(start = 1, end = 1)
table(a$type)
p1 <- ggplot(data = a, mapping = aes(x = V3, y = V4, colour = type))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c("F" = "#be0032", "M" = "#0067a5") )+
  labs(x= paste0("PC1(explained variance, ",pc1,"%)"), y = paste0("PC2(explained variance, ",pc2,"%)"))
p1
ggsave(filename = "plots/DNA_pipe/hj/hifi/final2/pca/allSNP.pdf", width = 4, height = 3, plot = p1)

# chr1 significant snp
a <- read.table("mid_data/DNA_pipe/HJ/hifi/final2/pca/chr1_signif/snp.eigenvec")
b <- read.table("mid_data/DNA_pipe/HJ/hifi/final2/pca/chr1_signif/snp.eigenval")
head(a)
head(b)
pc1 <- round(b$V1[1] / sum(b$V1) *100)
pc2 <- round(b$V1[2] / sum(b$V1) *100)
a$type <- a$V1 %>% str_sub(start = 1, end = 1)
table(a$type)
p1 <- ggplot(data = a, mapping = aes(x = V3, y = V4, colour = type))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c("F" = "#be0032", "M" = "#0067a5") )+
  labs(x= paste0("PC1(explained variance, ",pc1,"%)"), y = paste0("PC2(explained variance, ",pc2,"%)"))
p1
ggsave(filename = "plots/DNA_pipe/hj/hifi/final2/pca/chr1_siginifSNP.pdf", width = 4, height = 3, plot = p1)


# chr1 all snp
a <- read.table("mid_data/DNA_pipe/HJ/hifi/final2/pca/chr1_allSNP/snp.eigenvec")
b <- read.table("mid_data/DNA_pipe/HJ/hifi/final2/pca/chr1_allSNP/snp.eigenval")
head(a)
head(b)
pc1 <- round(b$V1[1] / sum(b$V1) *100)
pc2 <- round(b$V1[2] / sum(b$V1) *100)
a$type <- a$V1 %>% str_sub(start = 1, end = 1)
table(a$type)
p1 <- ggplot(data = a, mapping = aes(x = V3, y = V4, colour = type))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c("F" = "#be0032", "M" = "#0067a5") )+
  labs(x= paste0("PC1(explained variance, ",pc1,"%)"), y = paste0("PC2(explained variance, ",pc2,"%)"))
p1
ggsave(filename = "plots/DNA_pipe/hj/hifi/final2/pca/chr1_allSNP.pdf", width = 4, height = 3, plot = p1)




# chr2 all snp
a <- read.table("mid_data/DNA_pipe/HJ/hifi/final2/pca/chr2_allSNP/snp.eigenvec")
b <- read.table("mid_data/DNA_pipe/HJ/hifi/final2/pca/chr2_allSNP/snp.eigenval")
head(a)
head(b)
pc1 <- round(b$V1[1] / sum(b$V1) *100)
pc2 <- round(b$V1[2] / sum(b$V1) *100)
a$type <- a$V1 %>% str_sub(start = 1, end = 1)
table(a$type)
p1 <- ggplot(data = a, mapping = aes(x = V3, y = V4, colour = type))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c("F" = "#be0032", "M" = "#0067a5") )+
  labs(x= paste0("PC1(explained variance, ",pc1,"%)"), y = paste0("PC2(explained variance, ",pc2,"%)"))
p1
ggsave(filename = "plots/DNA_pipe/hj/hifi/final2/pca/chr2_allSNP.pdf", width = 4, height = 3, plot = p1)
