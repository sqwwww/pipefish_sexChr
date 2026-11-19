library(tidyverse)


XYallel <- read.table("mid_data/allele_specific/HJ/XYallele.table", header = T)
head(XYallel)
closetGene <- read.table("mid_data/DNA_pipe/HJ/hifi/final2/sex_snp/all_heterM_homoF_snp.bed.closest")
head(closetGene)
anno_df <- read.delim2("ref/eggnog/03.HJ/out.emapper (1).annotations", sep = "\t", comment.char = "#", header = F)
head(anno_df,1)
anno_df <- anno_df %>% select(c("V1","V9"))
head(anno_df)
anno_df$V3 <- anno_df$V1 %>% str_split(pattern = "\\.") %>% map_chr(1)
anno_df <- anno_df %>% select(c("V3","V9")) %>% distinct(V3, .keep_all = T)
names(anno_df) <- c("geneID","geneSymbol")
head(anno_df)
head(closetGene)
closetGene <- closetGene %>%  select(-c("V4","V5","V6","V9","V11"))
names(closetGene) <- c("chrom","snp_start","snp_end","genestart","geneEnd","geneStrand","geneID","distance")
closetGene <- closetGene %>% left_join(y = anno_df, by = "geneID")
head(closetGene)
table(is.na(closetGene$geneSymbol))
head(XYallel)
XYallel$chrom <- XYallel$variantID %>% str_split(pattern = ":") %>% map_chr(1) %>% as.integer()
XYallel$snp_start <- XYallel$variantID %>% str_split(pattern = ":") %>% map_chr(2)%>% as.integer()
XYallel <- XYallel %>% left_join(y = closetGene, by = c("chrom","snp_start"))
head(XYallel)

# for snps in two genes, only select fist line of the duplicates, I need to check it back in the fianl!!
XYallel <- XYallel %>% distinct(variantID, .keep_all = T)



combined_df <- data.frame()
for (my_tissue in c("Mba_Mt","Mba_Mb","Mba_Mp","Mba_Ms","Mba_Mk")){
  data1 <- read.table(paste0("mid_data/allele_specific/HJ/",my_tissue,".AS.v5.table"), header = T)
  data1 <- data1 %>% filter(contig == 1 & totalCount >0)
  # data1$ASEscore <- abs(data1$refCount/data1$totalCount -0.5) +0.5
  # expected_p <- 0.6
  # data1$binomial_p_value <- sapply(1:nrow(data1), function(i) {
  #   binom.test(x = data1$refCount[i], n = data1$totalCount[i], p = expected_p, alternative = "two.sided")$p.value
  # })
  # data1$adjusted_p_value <- p.adjust(data1$binomial_p_value, method = "BH")
  data2 <- data1
  data2 <- data2 %>% select(-c("variantID"))
  data2 <- data2 %>% left_join(y = XYallel, by = c("contig" = "chrom","position" = "snp_start"))
  data2$Xcount <- ifelse(test = data2$refAllele == data2$Xallele, yes = data2$refCount, no = data2$altCount)
  data2$Ycount <- ifelse(test = data2$refAllele == data2$Yallele, yes = data2$refCount, no = data2$altCount)
  # add strata infomation 8.9, 71.4
  data2$type <- ifelse(test = data2$position < 8.9e6, yes = "PAR", no = ifelse(test = data2$position < 71.4e6, yes = "young strata", no = "ancient strata"))
  #data2$type2 <- ifelse(test = data2$adjusted_p_value < 0.05, yes = "signif",no = "not signif")
  #data2$type3 <- ifelse(data2$Xcount > data2$Ycount, yes = "Xbiased", no = "Ybiased")
  #data2$type3[data2$adjusted_p_value > 0.05] <- "Not signif"
  data2$tissue <-  my_tissue
  combined_df <- plyr::rbind.fill(combined_df, data2)
}

head(combined_df)
table(combined_df$tissue)
combined_df$tissue <- combined_df$tissue %>% factor(levels = c("Mba_Mt","Mba_Mb","Mba_Mp","Mba_Ms","Mba_Mk"))



head(combined_df)

snpInExon <- read.table('mid_data/allele_specific/HJ/correctBiase_250427/ref/all_heterM_homoF_snp.m3_chr1_inExon.txt')
head(snpInExon)
snpInExon$location <- "exon"
combined_df <- combined_df %>% left_join(y = snpInExon, by = c("position" = "V2"))
combined_df$location[is.na(combined_df$location)] <- "others"
table(combined_df$location)

# check snps in exon
combined_df2 <- combined_df %>% filter(location == "exon")
head(combined_df2)
table(combined_df2$distance)

combined_df3 <- combined_df2 %>% group_by(tissue, geneID) %>% summarise(sXcount=sum(Xcount), sYcount=sum(Ycount), type =first(type), geneSymbol = first(geneSymbol), firstPosition=first(position), firstTissue = first(tissue), num = n()) %>% ungroup()
head(combined_df3)
table(combined_df3$tissue)
combined_df3$totalCount <- combined_df3$sXcount + combined_df3$sYcount
# ASE score for a gene
combined_df3$ASEscore <- abs(combined_df3$sXcount/combined_df3$totalCount -0.5) +0.5
expected_p <- 0.6
combined_df3$binomial_p_value <- sapply(1:nrow(combined_df3), function(i) {
  binom.test(x = combined_df3$sXcount[i], n = combined_df3$totalCount[i], p = expected_p, alternative = "two.sided")$p.value
})
combined_df3$adjusted_p_value <- p.adjust(combined_df3$binomial_p_value, method = "BH")
head(combined_df3)
combined_df3$type2 <- ifelse(test = combined_df3$adjusted_p_value < 0.05, yes = "signif",no = "not signif")
combined_df3$type3 <- ifelse(combined_df3$sXcount > combined_df3$sYcount, yes = "Xbiased", no = "Ybiased")
combined_df3$type3[combined_df3$adjusted_p_value > 0.05] <- "Not signif"
table(combined_df3$type3)
table(combined_df3$sXcount == 0)
combined_df3$logYXratio <- log2((combined_df3$sYcount+1)/(combined_df3$sXcount+1))

# set abs(logYXratio)>0.5 & signif as biased, else not biased
combined_df3$type4 <- combined_df3$type3
combined_df3$type4[combined_df3$type4 == "Not signif"] <- "Not biased"
combined_df3$type4[combined_df3$type3 != "Not signif" & abs(combined_df3$logYXratio) < 0.5] <- "Not biased"
table(combined_df3$type4)
table(combined_df3$type3)
combined_df4 <- combined_df3 %>% select(-c("type2", "type3","binomial_p_value", "firstTissue"))
head(combined_df4)
data.table::setnames(combined_df4, old = c("type", "type4"), new = c("location","type"))

write.table(x = combined_df4,file = "mid_data/allele_specific/HJ/260124/XYbiased_info_geneUnit.v3.rmdup.table",quote = F, row.names = F, col.names = T, sep = "\t")


# check with previous version, the result is the same, no need to do any change
# data1 <- read.table("mid_data/allele_specific/HJ/correctBiase_250427/XYbiased_info_geneUnit.table", header = T, sep = "\t")
# head(data1)

combined_df4 <- read.table(file = "mid_data/allele_specific/HJ/correctBiase_250427/XYbiased_info_geneUnit.v3.rmdup.table", header = T, sep = "\t")
head(combined_df4)
table(combined_df4$type)
# plot ratio bar for each tissue 
table(combined_df4$tissue)
combined_df4 <- combined_df4 %>% filter(!tissue %in% c("Mba_Mp","Mba_Ms"))
combined_df4$tissue <- combined_df4$tissue %>% factor(levels = c("Mba_Mt", "Mba_Mb", "Mba_Mk"))
p1 <- ggplot(data = combined_df4, mapping = aes(x = tissue, fill = type))+
  geom_bar(width = 0.8, position = "fill")+
  scale_fill_manual(values = c("Xbiased" ="#be0032" , "Ybiased" = "#0067a5", "Not biased" = "#bebebe"))+ theme_classic()
p1
ggsave(filename = "plots/DNA_pipe/hj/hifi/final2/allele_specific/260124/geneTypeRatio.v3.pdf", plot = p1, width = 5, height = 3)

# ring plot
plot_data <- combined_df4 %>% filter(tissue == "Mba_Mt") %>% group_by(type) %>% summarise(n=n())
head(plot_data)
# Make the plot
hsize <- 4
p1 <- ggplot(plot_data, aes(x = hsize, y = n, fill = type)) +
  geom_col() +
  coord_polar(theta = "y") +
  xlim(c(1.5, hsize+0.5))+
  scale_fill_manual(values = c("Xbiased" ="#be0032" , "Ybiased" = "#0067a5", "Not biased" = "#bebebe"))+ theme_void()
p1
ggsave(filename = "plots/DNA_pipe/hj/hifi/final2/allele_specific/260124/ring_YbiasedGeneInTestis.v3.pdf", plot = p1, width = 5, height = 3)

p2 <- ggplot(data = combined_df4 %>% filter(tissue == "Mba_Mb") %>% group_by(type) %>% summarise(n=n()), aes(x = hsize, y = n, fill = type)) +
  geom_col() +
  coord_polar(theta = "y") +
  xlim(c(1.5, hsize+0.5))+
  scale_fill_manual(values = c("Xbiased" ="#be0032" , "Ybiased" = "#0067a5", "Not biased" = "#bebebe"))+ theme_void()
p2
ggsave(filename = "plots/DNA_pipe/hj/hifi/final2/allele_specific/260124/ring_YbiasedGeneInBrain.v3.pdf", plot = p2, width = 5, height = 3)
p3 <- ggplot(data = combined_df4 %>% filter(tissue == "Mba_Mk") %>% group_by(type) %>% summarise(n=n()), aes(x = hsize, y = n, fill = type)) +
  geom_col() +
  coord_polar(theta = "y") +
  xlim(c(1.5, hsize+0.5))+
  scale_fill_manual(values = c("Xbiased" ="#be0032" , "Ybiased" = "#0067a5", "Not biased" = "#bebebe"))+ theme_void()
p3
ggsave(filename = "plots/DNA_pipe/hj/hifi/final2/allele_specific/260124/ring_YbiasedGeneInKindey.v3.pdf", plot = p3, width = 5, height = 3)
p4 <- ggplot(data = combined_df4 %>% filter(tissue == "Mba_Mp") %>% group_by(type) %>% summarise(n=n()), aes(x = hsize, y = n, fill = type)) +
  geom_col() +
  coord_polar(theta = "y") +
  xlim(c(1.5, hsize+0.5))+
  scale_fill_manual(values = c("Xbiased" ="#be0032" , "Ybiased" = "#0067a5", "Not biased" = "#bebebe"))+ theme_void()
p4
ggsave(filename = "plots/DNA_pipe/hj/hifi/final2/allele_specific/260124/ring_YbiasedGeneInPouch.v3.pdf", plot = p4, width = 5, height = 3)
p5 <- ggplot(data = combined_df4 %>% filter(tissue == "Mba_Ms") %>% group_by(type) %>% summarise(n=n()), aes(x = hsize, y = n, fill = type)) +
  geom_col() +
  coord_polar(theta = "y") +
  xlim(c(1.5, hsize+0.5))+
  scale_fill_manual(values = c("Xbiased" ="#be0032" , "Ybiased" = "#0067a5", "Not biased" = "#bebebe"))+ theme_void()
p5
ggsave(filename = "plots/DNA_pipe/hj/hifi/final2/allele_specific/260124/ring_YbiasedGeneInSkin.v3.pdf", plot = p5, width = 5, height = 3)


# calculate Y baised genes ratio
table(combined_df4$tissue)
table(combined_df4 %>% filter(tissue == "Mba_Mt") %>% pull(type))
prop.table(table(combined_df4 %>% filter(tissue == "Mba_Mt") %>% pull(type)))
prop.table(table(combined_df4 %>% filter(tissue == "Mba_Mb") %>% pull(type)))
prop.table(table(combined_df4 %>% filter(tissue == "Mba_Mk") %>% pull(type)))
prop.table(table(combined_df4 %>% filter(tissue == "Mba_Mp") %>% pull(type)))
prop.table(table(combined_df4 %>% filter(tissue == "Mba_Ms") %>% pull(type)))



# extract gene symbol for metascape
head(combined_df4)
a <- combined_df4 %>% filter(tissue =="Mba_Mt", type=="Ybiased")
head(a)
write.table(x = a$geneSymbol %>% setdiff(y = "-") %>% na.omit(), 
            file = "mid_data/allele_specific/HJ/260124/testisYbiasedGenes.v3.txt", quote = F, col.names = F, row.names = F)
b <- combined_df4 %>% filter(tissue =="Mba_Mt", type=="Not biased")
write.table(x = b$geneSymbol %>% setdiff(y = "-") %>% na.omit(), 
            file = "mid_data/allele_specific/HJ/260124/testisNotbiasedGenes.v3.txt", quote = F, col.names = F, row.names = F)
sink("mid_data/allele_specific/HJ/260124/testisYbiasedGenes_NotbiasedGenes.v3.metascape.txt")
paste0(a$geneSymbol %>% setdiff(.,"-"), collapse = ",")
paste0(b$geneSymbol %>% setdiff(.,"-"), collapse = ",")
sink()

tmp <- read.table("mid_data/allele_specific/HJ/correctBiase_250427/XYbiasedGene/testisYbiasedGenes.v2.txt")
head(tmp)


# Ybiased genes GO enrichment
phenotype_df <- readxl::read_excel("mid_data/allele_specific/HJ/260124/testisYbiasedGenes_NotbiasedGenes.v3.metascape/selectGOterm.xlsx", sheet = 1)
head(phenotype_df)
# phenotype_df <- phenotype_df %>% arrange(LogP)
phenotype_df$Description <- phenotype_df$Description %>% factor(levels =rev(phenotype_df$Description ))
plot_data2 <- phenotype_df %>% select(c("_LogP_Not","_LogP_Ybias","Description"))
plot_data3 <- plot_data2 %>% reshape2::melt(id.vars = "Description", measure.vars = c("_LogP_Not","_LogP_Ybias"))
plot_data3$variable <- plot_data3$variable %>% factor(levels = c("_LogP_Ybias","_LogP_Not"))
head(plot_data3)
write.table(x = plot_data3, file = "mid_data/allele_specific/HJ/260124/plot_data.table", row.names = F, col.names = T, sep = "\t", quote = F)
# manual add testosterone term 
plot_data3 <- read.table("mid_data/allele_specific/HJ/260124/plot_data.table", header = T, sep="\t")
head(plot_data3)
plot_data3$variable <- plot_data3$variable %>% factor(levels = c("_LogP_Ybias","_LogP_Not"))
plot_data3$Description <- plot_data3$Description %>% factor(levels =c("Testosterone",rev(as.character(phenotype_df$Description))))
library(ggthemes)
p1 <- ggplot(data = plot_data3, mapping = aes(x = variable,y = Description))+
  geom_point(mapping = aes(color = -value, size = -value))+
  scale_color_gradient2(low = "#0067a5", mid = "#f3f3f4", high = "#be0032", midpoint = -log10(0.05))+
  theme_gdocs()
p1

ggsave(filename = "plots/DNA_pipe/hj/hifi/final2/allele_specific/260124/Ybias_NotBias_GOenrich.v2.pdf", plot =p1, width = 5, height = 2.5) 

