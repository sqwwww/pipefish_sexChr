library(tidyverse)


file_info <- read.table(file = "mid_data/genomeAssembly/GCspike/files_info.txt")
head(file_info)
ratio_df <- data.frame()
for (i in seq_len(nrow(file_info))){
  current_sp <- file_info[i,"V1"]
  current_pouch_complex <- file_info[i,"V5"]
  current_GCspike_ratio_file <- paste0("mid_data/genomeAssembly/GCspike/",current_sp,".spike.tipRatio2.txt")
  tmp_file <- read.table(current_GCspike_ratio_file, header = 1)
  tmp_file <- tmp_file %>% drop_na()
  tmp_file$type <- current_pouch_complex
  ratio_df <- ratio_df %>% plyr::rbind.fill(tmp_file)
}
table(ratio_df$species)
head(ratio_df)
table(ratio_df$type)

a <- ratio_df %>% filter(type == 2 | type == 3)
b <- ratio_df %>% filter(type == 1)
prop.table(table(a$species, a$Ratio > 1.1), margin = 1)
prop.table(table(b$species, b$Ratio > 1.1), margin = 1)
# file_info %>% arrange(V5) %>% pull(V1)

library(ggpubr)

guppy_mean <- mean(ratio_df %>% filter(species == "Pret") %>% pull(Ratio))
library(ggthemes)
library(ggtext)

my_type_order <- c("Pret","Onil","Ha","He","Ssco","Hhep","Cint","Ptae","Sbia","Ddac","Mbra","Eaeq")
my_type_color <- c("Ha"="#f6a600",
                   "He"="#f6a600",
                   "Ssco" = "#f6a600",
                   "Hhep" = "#f6a600",
                   "Ptae" = "#dcd300",
                   "Pequ" = "#dcd300",
                   "Sbia" = "#dcd300",
                   "Cint" = "#dcd300",
                   "Ddac" = "#0067a5",
                   "Mbra" = "#0067a5",
                   "Eaeq" = "#0067a5",
                   "Pret" = "grey",
                   "Onil" = "grey")
my_comparisons <- combn(my_type_order, 2, simplify = FALSE)
table(ratio_df$type)
ratio_df$type <- ratio_df$type %>% as.factor()

# GC value per genome read
GCperGenome <- readxl::read_xlsx(path = "mid_data/genomeAssembly/GCspike/GC_stats.xlsx", sheet = 1)
head(GCperGenome)
GCperGenome$species <- factor(GCperGenome$species, levels = c('Pret',"Onil","Eaeq",'Mbra','Ddac','Cint','Sbia','Pequ','Ptae','Hhep','Ssco','Ha','He'))
GCperGenome$GCcontent <- GCperGenome$GCcontent/100
GCperGenome <- GCperGenome[order(GCperGenome$species), ]

ratio_df$species <- ratio_df$species %>% 
  factor(levels = c('Pret',"Onil","Eaeq",'Mbra','Ddac','Cint','Sbia','Pequ','Ptae','Hhep','Ssco','Ha','He'))
my_comparisons <- combn(file_info$V1, 2, simplify = FALSE)

p2 <- ggplot(data = ratio_df, mapping = aes(y = species, x = Ratio))+
  geom_boxplot(aes(fill = species, color =species),outlier.shape =NA)+
  labs(y="GC spike Ratio")+
  # geom_vline(xintercept = guppy_mean,linetype="3313",color="black",linewidth=1)+
  labs(y="GC spike Ratio")+
  scale_fill_manual(values=my_type_color)+
  scale_color_manual(values=my_type_color) + 
  theme_classic()
p3 <- ggplot(data = ratio_df, mapping = aes(y = species, x = background))+
  geom_boxplot(aes(fill = species, color =species),outlier.shape =NA)+
  labs(y="Background GC")+
  # geom_vline(xintercept = guppy_mean,linetype="3313",color="black",linewidth=1)+
  scale_fill_manual(values=my_type_color)+
  scale_color_manual(values=my_type_color) + 
  theme_classic()
p3 <- ggplot(data = ratio_df, aes(y = species, x = background)) +
  geom_boxplot(aes(fill = species, color = species), outlier.shape = NA) +
  geom_point(data = GCperGenome, 
             aes(x = GCcontent, y = species),
             color = "black", size = 2) +
  geom_line(data = GCperGenome, 
            aes(x = GCcontent, y = species, group = 1),
            color = "black", linewidth = 0.8) +
  labs(y = "Background GC") +
  scale_fill_manual(values = my_type_color) +
  scale_color_manual(values = my_type_color) +
  theme_classic()
p3 <- ggplot(ratio_df, aes(y = species, x = background)) +
  geom_boxplot(aes(fill = species, color = species), outlier.shape = NA) +
  geom_point(data = GCperGenome,
             aes(x = GCcontent, y = species),
             shape = 5,
             color = "black",
             size = 1.8,
             stroke = 0.7) +
  
  labs(x = "background", y = "Background GC") +
  scale_fill_manual(values = my_type_color) +
  scale_color_manual(values = my_type_color) +
  theme_classic()


p3
ggsave(plot = p2,width =3,height = 3,filename = "plots/基因组组装/GC_spike/final/species.v6.pdf")
ggsave(plot = p3,width =3,height = 3,filename = "plots/基因组组装/GC_spike/final/species.backgroundGC.v2.pdf")
my_comparisons <- combn(file_info$V1, 2, simplify = FALSE)
head(ratio_df)
i = 1
test_out <- data.frame()
for(i in seq_len(length(my_comparisons))){
  sp1 <- my_comparisons[[i]][1]
  sp2 <- my_comparisons[[i]][2]
  sp1_data <- ratio_df %>% filter(species == sp1)
  sp2_data <- ratio_df %>% filter(species == sp2)
  testResult <- wilcox.test(x = sp1_data$Ratio, y = sp2_data$Ratio)
  tmp_sig <- ifelse(testResult$p.value<0.05, yes = "yes","no")
  tmp_df <- data.frame("sp1" = sp1, "sp2" = sp2,
                       "Pvalue" = testResult$p.value, "Sig"=tmp_sig)
  test_out <- plyr::rbind.fill(test_out,tmp_df)
}

# closed pouch vs open pouch
test_out %>% filter((sp1 %in% c("Ha","He","Hhep","Ssco")) | (sp2 %in% c("Ha","He","Hhep","Ssco"))) %>% 
  filter((sp1 %in% c("Mbra","Ddac","Eaeq")) | (sp2 %in% c("Mbra","Ddac","Eaeq")))
# closed pouch vs semi pouch
test_out %>% filter((sp1 %in% c("Ha","He","Hhep","Ssco")) | (sp2 %in% c("Ha","He","Hhep","Ssco"))) %>% 
  filter((sp1 %in% c("Pequ","Ptae","Sbia","Cint")) | (sp2 %in% c("Pequ","Ptae","Sbia","Cint")))



my_comparisons2 <- combn(unique(file_info$V5), 2, simplify = FALSE)
test_out2 <- data.frame()
for(i in seq_len(length(my_comparisons2))){
  sp1 <- my_comparisons2[[i]][1]
  sp2 <- my_comparisons2[[i]][2]
  sp1_data <- ratio_df %>% filter(type == sp1)
  sp2_data <- ratio_df %>% filter(type == sp2)
  testResult <- wilcox.test(x = sp1_data$Ratio, y = sp2_data$Ratio)
  tmp_sig <- ifelse(testResult$p.value<0.05, yes = "yes","no")
  tmp_df <- data.frame("sp1" = sp1, "sp2" = sp2,
                       "Pvalue" = testResult$p.value, "Sig"=tmp_sig)
  test_out2 <- plyr::rbind.fill(test_out2,tmp_df)
}

# closed pouch vs open pouch
test_out2 %>% filter((sp1 %in% c("1")) | (sp2 %in% c("1"))) %>% 
  filter((sp1 %in% c("3")) | (sp2 %in% c("3")))
# closed pouch vs semi pouch
test_out2 %>% filter((sp1 %in% c("1")) | (sp2 %in% c("1"))) %>% 
  filter((sp1 %in% c("2")) | (sp2 %in% c("2")))
test_out2 %>% filter((sp1 %in% c("3")) | (sp2 %in% c("3"))) %>% 
  filter((sp1 %in% c("2")) | (sp2 %in% c("2")))


# for background GC value
my_comparisons2 <- combn(unique(file_info$V5), 2, simplify = FALSE)
test_out3 <- data.frame()
for(i in seq_len(length(my_comparisons2))){
  sp1 <- my_comparisons2[[i]][1]
  sp2 <- my_comparisons2[[i]][2]
  sp1_data <- ratio_df %>% filter(type == sp1)
  sp2_data <- ratio_df %>% filter(type == sp2)
  testResult <- wilcox.test(x = sp1_data$background, y = sp2_data$background)
  tmp_sig <- ifelse(testResult$p.value<0.05, yes = "yes","no")
  tmp_df <- data.frame("sp1" = sp1, "sp2" = sp2,
                       "Pvalue" = testResult$p.value, "Sig"=tmp_sig)
  test_out3 <- plyr::rbind.fill(test_out3,tmp_df)
}

# closed pouch vs open pouch
test_out3 %>% filter((sp1 %in% c("1")) | (sp2 %in% c("1"))) %>% 
  filter((sp1 %in% c("3")) | (sp2 %in% c("3")))
# open pouch vs semi pouch
test_out3 %>% filter((sp1 %in% c("1")) | (sp2 %in% c("1"))) %>% 
  filter((sp1 %in% c("2")) | (sp2 %in% c("2")))
# closed pouch vs semi pouch
test_out3 %>% filter((sp1 %in% c("3")) | (sp2 %in% c("3"))) %>% 
  filter((sp1 %in% c("2")) | (sp2 %in% c("2")))



# plot the telomere infomation
files <- Sys.glob("mid_data/genomeAssembly/telomere/*stats2.txt") 
telo_df <- data.frame()
for (i in files){
  current_sp <- basename(i) %>% str_remove(".stats2.txt")
  tmp_file <- read.table(i, header = 1)
  twoNo <- sum(tmp_file$head_telomere == "No" & tmp_file$tail_telomere == "No")
  twoYes <- sum(tmp_file$head_telomere == "Yes" & tmp_file$tail_telomere == "Yes")
  oneNo <- nrow(tmp_file) - twoNo - twoYes
  tmp_row <- data.frame(c(twoYes,oneNo, twoNo, current_sp)) %>% t() %>% as.data.frame()
  telo_df <- telo_df %>% plyr::rbind.fill(tmp_row)
}
telo_df

names(telo_df) <- c("twoYes","oneNo", "twoNo", "Species")
telo_df$Species <- telo_df$Species %>% factor(levels = c("Her","Ha","Ssc","HF","Pta","Peq","Sbi","papa","BJ","HJ","Eae","Onil","guppy") %>% rev())
telo_df2 <- telo_df %>% reshape2::melt( measure.vars = c("twoYes","oneNo", "twoNo"), variable.name = "type")
telo_df2$value <- telo_df2$value %>% as.numeric()
p4 <- ggplot(data = telo_df2, mapping = aes(x = value, y = Species, fill= type))+
  geom_bar(stat = "identity",position = "fill")+
  labs(x="telomoere assembly completeness")+
  scale_fill_manual(values=c("twoYes"="#c2b280","oneNo"="#848482", "twoNo"="lightgrey"))+
  theme_minimal()
p4 
p4 <- ggplot(telo_df2, aes(x = value, y = Species, fill = type)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(
    data = subset(telo_df2, value > 0),
    aes(label = value),
    position = position_fill(vjust = 0.5),
    size = 3
  ) +
  labs(x = "telomere assembly completeness") +
  scale_fill_manual(values = c("twoYes" = "#c2b280",
                               "oneNo" = "#848482",
                               "twoNo"  = "lightgrey")) +
  theme_minimal()
ggsave(filename = "plots/基因组组装/telomere/telomereCompleteness.pdf", plot = p4, width = 5, height = 6)
  


