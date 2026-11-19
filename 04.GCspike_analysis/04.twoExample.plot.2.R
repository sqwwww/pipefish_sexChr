library(tidyverse)
library(changepoint)
library(plyr)


# for the Mb_chr14 GC spike plot
data <- read.table("mid_data/基因组组装/GCspike/BJ.GCcontent.100k.bed")
head(data)
data$chrom %>% table()
names(data) <- c("chrom","start","end","value")
data$value %>% fivenum()

chrom_info <- read.table("mid_data/基因组组装/GCspike/ringedPipefish.chromosome.fasta.masked.fai")
head(chrom_info)
names(chrom_info) <- c("chrom","length","N1","N2","N3")
chrom_info <- chrom_info %>% select(c("chrom","length"))

data[data$value>0.6,"value"] <- 0.6

current_chr <- 19
cp_out_df <- data.frame()
myJudge <- c()
tipRatio <- c()
for (current_chr in seq_len(nrow(chrom_info))){
  current_data <- data %>% filter(chrom == current_chr)
  fit_changepoint = cpt.meanvar(current_data$value, method = "BinSeg",Q = 3)
  cp_coord <- c(0,cpts(fit_changepoint)/10,chrom_info[current_chr,"length"]/1e6)
  cp_values<- param.est(fit_changepoint)$mean
  if(which.max(cp_values) %in% c(1,length(cp_values))){
    print(paste(current_chr,"have tip peak."))
    tip_region_length <- cp_coord[(which.max(cp_values)+1)]-cp_coord[which.max(cp_values)]
    tipJudge <- "Y"
    if ((tip_region_length * 3) < (chrom_info[current_chr,"length"]/1e6)){
      print(paste(current_chr,"s tip region have length less than 1/5 chr length"));
      print(tip_region_length)
      lengthJudge <- "Y";
    }else{lengthJudge <- "N"}
  }else{
    tipJudge <- "N"
  }
  
  my_max_GC <- max(cp_values,na.rm = T)
  my_longest_GC <- cp_values[which.max(diff(cp_coord))]
  
  if((tipJudge == "Y") & (lengthJudge == "Y") & my_max_GC > 0.4195){
    myJudge[current_chr] <- "Y"
    tipRatio[current_chr] <- my_max_GC/my_longest_GC
  }else{
    myJudge[current_chr] <- "N"
  }
  
  cp_df <- data.frame(
    X1=cp_coord[1:length(cp_coord)-1],
    XE=cp_coord[2:length(cp_coord)],
    Y1=cp_values,
    YE=cp_values,
    chrom=current_chr,
    tipJudge=tipJudge,
    lengthJudge = lengthJudge)
  cp_out_df <- rbind.fill(cp_out_df,cp_df)
  
}
myJudge
tipRatio
tipRatio_df <- reshape2::melt(tipRatio, value.name = "Ratio") %>% 
  mutate(chrom = rownames(.)) %>%
  mutate(species = "Mbra") %>% 
  select(c(species,chrom,Ratio))
head(data)
data1 <- data
cp_out_df1 <- cp_out_df



# for Hh_chr22 ---
data <- read.table("mid_data/基因组组装/GCspike/HF.GCcontent.100k.bed")
head(data)
names(data) <- c("chrom","start","end","value")
data$value %>% fivenum()
my_min_value = min(data$value)

chrom_info <- read.table("mid_data/基因组组装/GCspike/HF.genome.chr.fasta.hardMasked.fai")
head(chrom_info)
names(chrom_info) <- c("chrom","length","N1","N2","N3")
chrom_info <- chrom_info %>% select(c("chrom","length"))




current_chr <- 19
cp_out_df <- data.frame()
myJudge <- c()
tipRatio <- c()
for (current_chr in seq_len(nrow(chrom_info))){
  current_data <- data %>% filter(chrom == current_chr)
  fit_changepoint = cpt.meanvar(current_data$value, method = "BinSeg",Q = 3)
  cp_coord <- c(0,cpts(fit_changepoint)/10,chrom_info[current_chr,"length"]/1e6)
  cp_values<- param.est(fit_changepoint)$mean
  if(which.max(cp_values) %in% c(1,length(cp_values))){
    print(paste(current_chr,"have tip peak."))
    tip_region_length <- cp_coord[(which.max(cp_values)+1)]-cp_coord[which.max(cp_values)]
    tipJudge <- "Y"
    if ((tip_region_length * 3) < (chrom_info[current_chr,"length"]/1e6)){
      print(paste(current_chr,"s tip region have length less than 1/3 chr length"));
      print(tip_region_length)
      lengthJudge <- "Y";
    }else{lengthJudge <- "N"}
  }else{
    tipJudge <- "N"
  }
  
  my_max_GC <- max(cp_values,na.rm = T)
  my_longest_GC <- cp_values[which.max(diff(cp_coord))]
  
  if((tipJudge == "Y") & (lengthJudge == "Y") & my_max_GC > 0.39){
    myJudge[current_chr] <- "Y"
    tipRatio[current_chr] <- my_max_GC/my_longest_GC
  }else{
    myJudge[current_chr] <- "N"
  }
  
  cp_df <- data.frame(
    X1=cp_coord[1:length(cp_coord)-1],
    XE=cp_coord[2:length(cp_coord)],
    Y1=cp_values,
    YE=cp_values,
    chrom=current_chr,
    tipJudge=tipJudge,
    lengthJudge = lengthJudge)
  cp_out_df <- rbind.fill(cp_out_df,cp_df)
  
}
myJudge
tipRatio
tipRatio_df <- reshape2::melt(tipRatio, value.name = "Ratio") %>% 
  mutate(chrom = rownames(.)) %>%
  mutate(species = "Hhep") %>% 
  select(c(species,chrom,Ratio))
data2 <- data
cp_out_df2 <- cp_out_df


# for papa

data <- read.table("mid_data/基因组组装/GCspike/Cin.GCcontent.100k.bed")
head(data)
names(data) <- c("chrom","start","end","value")
data$value %>% fivenum()
my_min_value = min(data$value)

chrom_info <- read.table("mid_data/基因组组装/GCspike/Cin.numChr.fasta.masked.fai")
head(chrom_info)
names(chrom_info) <- c("chrom","length","N1","N2","N3")
chrom_info <- chrom_info %>% select(c("chrom","length"))




current_chr <- 19
cp_out_df <- data.frame()
myJudge <- c()
tipRatio <- c()
for (current_chr in seq_len(nrow(chrom_info))){
  current_data <- data %>% filter(chrom == current_chr)
  fit_changepoint = cpt.meanvar(current_data$value, method = "BinSeg",Q = 3)
  cp_coord <- c(0,cpts(fit_changepoint)/10,chrom_info[current_chr,"length"]/1e6)
  cp_values<- param.est(fit_changepoint)$mean
  if(which.max(cp_values) %in% c(1,length(cp_values))){
    print(paste(current_chr,"have tip peak."))
    tip_region_length <- cp_coord[(which.max(cp_values)+1)]-cp_coord[which.max(cp_values)]
    tipJudge <- "Y"
    if ((tip_region_length * 3) < (chrom_info[current_chr,"length"]/1e6)){
      print(paste(current_chr,"s tip region have length less than 1/3 chr length"));
      print(tip_region_length)
      lengthJudge <- "Y";
    }else{lengthJudge <- "N"}
  }else{
    tipJudge <- "N"
  }
  
  my_max_GC <- max(cp_values,na.rm = T)
  my_longest_GC <- cp_values[which.max(diff(cp_coord))]
  
  if((tipJudge == "Y") & (lengthJudge == "Y") & my_max_GC > 0.39){
    myJudge[current_chr] <- "Y"
    tipRatio[current_chr] <- my_max_GC/my_longest_GC
  }else{
    myJudge[current_chr] <- "N"
  }
  
  cp_df <- data.frame(
    X1=cp_coord[1:length(cp_coord)-1],
    XE=cp_coord[2:length(cp_coord)],
    Y1=cp_values,
    YE=cp_values,
    chrom=current_chr,
    tipJudge=tipJudge,
    lengthJudge = lengthJudge)
  cp_out_df <- rbind.fill(cp_out_df,cp_df)
  
}
myJudge
tipRatio
tipRatio_df <- reshape2::melt(tipRatio, value.name = "Ratio") %>% 
  mutate(chrom = rownames(.)) %>%
  mutate(species = "Hhep") %>% 
  select(c(species,chrom,Ratio))
data3 <- data
cp_out_df3 <- cp_out_df




# for HJ

data <- read.table("mid_data/genomeAssembly//GCspike/HJ.GCcontent.100k.bed")
head(data)
names(data) <- c("chrom","start","end","value")
data$value %>% fivenum()
my_min_value = min(data$value)

chrom_info <- read.table("mid_data/genomeAssembly/GCspike/HJ.genome.chr.fasta.hardMasked.fai")
head(chrom_info)
names(chrom_info) <- c("chrom","length","N1","N2","N3")
chrom_info <- chrom_info %>% select(c("chrom","length"))




current_chr <- 19
cp_out_df <- data.frame()
myJudge <- c()
tipRatio <- c()
for (current_chr in seq_len(nrow(chrom_info))){
  current_data <- data %>% filter(chrom == current_chr)
  fit_changepoint = cpt.meanvar(current_data$value, method = "BinSeg",Q = 3)
  cp_coord <- c(0,cpts(fit_changepoint)/10,chrom_info[current_chr,"length"]/1e6)
  cp_values<- param.est(fit_changepoint)$mean
  if(which.max(cp_values) %in% c(1,length(cp_values))){
    print(paste(current_chr,"have tip peak."))
    tip_region_length <- cp_coord[(which.max(cp_values)+1)]-cp_coord[which.max(cp_values)]
    tipJudge <- "Y"
    if ((tip_region_length * 3) < (chrom_info[current_chr,"length"]/1e6)){
      print(paste(current_chr,"s tip region have length less than 1/3 chr length"));
      print(tip_region_length)
      lengthJudge <- "Y";
    }else{lengthJudge <- "N"}
  }else{
    tipJudge <- "N"
  }
  
  my_max_GC <- max(cp_values,na.rm = T)
  my_longest_GC <- cp_values[which.max(diff(cp_coord))]
  
  if((tipJudge == "Y") & (lengthJudge == "Y") & my_max_GC > 0.39){
    myJudge[current_chr] <- "Y"
    tipRatio[current_chr] <- my_max_GC/my_longest_GC
  }else{
    myJudge[current_chr] <- "N"
  }
  
  cp_df <- data.frame(
    X1=cp_coord[1:length(cp_coord)-1],
    XE=cp_coord[2:length(cp_coord)],
    Y1=cp_values,
    YE=cp_values,
    chrom=current_chr,
    tipJudge=tipJudge,
    lengthJudge = lengthJudge)
  cp_out_df <- rbind.fill(cp_out_df,cp_df)
  
}
myJudge
tipRatio
tipRatio_df <- reshape2::melt(tipRatio, value.name = "Ratio") %>% 
  mutate(chrom = rownames(.)) %>%
  mutate(species = "Mbr") %>% 
  select(c(species,chrom,Ratio))
data4 <- data
cp_out_df4 <- cp_out_df



p1 <- ggplot(data = data1 %>% filter(chrom == "18"),
             mapping = aes(x = start/1e6, y = value))+
  geom_ribbon(aes(ymin = 0.3, ymax = value), fill = "#0067a5")+
  geom_segment(data = cp_out_df1 %>% filter(chrom ==18),aes(x = X1, xend = XE, 
                                                            y = Y1, yend = YE),
               colour = "black", linewidth = 1)+
  facet_wrap(~chrom, scales = "free") +
  lims(y = c(0.3,0.6)) +
  labs(x= "Position(Mb)",y="GC content")+
  theme_classic()
p1
p2 <- ggplot(data = data2 %>% filter(chrom == 22),
             mapping = aes(x = start/1e6, y = value))+
  geom_ribbon(aes(ymin = 0.3, ymax = value), fill = "#0067a5")+
  geom_segment(data = cp_out_df2 %>% filter(chrom == 22),aes(x = X1, xend = XE, 
                                    y = Y1, yend = YE),
               colour = "black", linewidth = 1)+
  facet_wrap(~chrom, scales = "free") +
  lims(y = c(0.3,0.6)) +
  labs(x= "Position(Mb)",y="GC content") +
  theme_classic()

p2

p3 <- ggplot(data = data3 %>% filter(chrom == 20),
             mapping = aes(x = start/1e6, y = value))+
  geom_ribbon(aes(ymin = 0.3, ymax = value), fill = "#0067a5")+
  geom_segment(data = cp_out_df3 %>% filter(chrom == 20),aes(x = X1, xend = XE, 
                                                             y = Y1, yend = YE),
               colour = "black", linewidth = 1)+
  facet_wrap(~chrom, scales = "free") +
  lims(y = c(0.3,0.6)) +
  labs(x= "Position(Mb)",y="GC content") +
  theme_classic()

p3

p4 <- ggplot(data = data4 %>% filter(chrom == 19),
             mapping = aes(x = start/1e6, y = value))+
  geom_ribbon(aes(ymin = 0.3, ymax = value), fill = "#0067a5")+
  geom_segment(data = cp_out_df4 %>% filter(chrom == 19),aes(x = X1, xend = XE, 
                                                             y = Y1, yend = YE),
               colour = "black", linewidth = 1)+
  facet_wrap(~chrom, scales = "free") +
  lims(y = c(0.3,0.6)) +
  labs(x= "Position(Mb)",y="GC content") +
  theme_classic()

p4

ggsave(plot = p4,filename = "plots/基因组组装/GC_spike/HJ_examples.pdf",width = 8,height = 2.5)

library(ggpubr)
p <- ggarrange(p2,p3,p1, nrow = 1)

ggsave(plot = p,filename = "plots/基因组组装/GC_spike/three_examples.v2.pdf",width = 8,height = 2.5)

