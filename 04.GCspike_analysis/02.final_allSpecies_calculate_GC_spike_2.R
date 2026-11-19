library(tidyverse)
library(changepoint)
library(plyr)

file_info <- read.table(file = "mid_data/genomeAssembly//GCspike/files_info.txt")
head(file_info)

i <- 1
for (i in seq_len(nrow(file_info))){
  current_GC_file <- file_info[i,"V2"]
  current_fai_file <- file_info[i,"V3"]
  current_sp <- file_info[i,"V1"]
  ave_GC <- file_info[i,"V4"]
  
  data <- read.table(paste0("mid_data/genomeAssembly//GCspike/",current_GC_file))
  head(data)
  names(data) <- c("chrom","start","end","value")

  chrom_info <- read.table(paste0("mid_data/genomeAssembly/GCspike/",current_fai_file))
  head(chrom_info)
  names(chrom_info) <- c("chrom","length","N1","N2","N3")
  chrom_info <- chrom_info %>% select(c("chrom","length"))

  if(current_sp == "Onil"){
    data$chrom <- data$chrom %>% gsub(pattern = "Onil_LG", replacement = "") %>% as.integer()
    chrom_info$chrom <- chrom_info$chrom %>% gsub(pattern = "Onil_LG", replacement = "") %>% as.integer() 
    chrom_info <- arrange(df = chrom_info, chrom)
  }
  
  
# current_chr <- 19
  cp_out_df <- data.frame()
  myJudge <- c()
  tipRatio <- c()
  tmp_longest_GC <- c()
  for (current_chr in seq_len(nrow(chrom_info))){
    if(current_sp == "Onil" & current_chr> 20){current_chr = current_chr+1}
    current_data <- data %>% filter(chrom == current_chr)
    fit_changepoint = cpt.meanvar(current_data$value, method = "BinSeg",Q = 3)
    cp_coord <- c(0,cpts(fit_changepoint)/10,chrom_info[current_chr,"length"]/1e6)
    cp_values<- param.est(fit_changepoint)$mean
    if(which.max(cp_values) %in% c(1,length(cp_values))){
      print(paste(current_chr,"have tip peak."))
      tip_region_length <- cp_coord[(which.max(cp_values)+1)]-cp_coord[which.max(cp_values)]
      tipJudge <- "Y"
      if ((tip_region_length * 2) < (chrom_info[current_chr,"length"]/1e6)){
        print(paste(current_chr,"s tip region have length less than 1/2 chr length"));
        print(tip_region_length)
        lengthJudge <- "Y";
      }else{lengthJudge <- "N"}
    }else{
      tipJudge <- "N";
      lengthJudge <- "N"
    }
    
    my_max_GC <- max(cp_values,na.rm = T)
    my_longest_GC <- cp_values[which.max(diff(cp_coord))]
    
    if((tipJudge == "Y") & (lengthJudge == "Y") & (my_max_GC > my_longest_GC)){
      myJudge[current_chr] <- "Y"
      tipRatio[current_chr] <- my_max_GC/my_longest_GC
      tmp_longest_GC[current_chr] <- my_longest_GC
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
    mutate(species = current_sp) %>% 
    select(c(species,chrom,Ratio))
  tipRatio_df$background <- tmp_longest_GC
  write.table(x = tipRatio_df, file = paste0("mid_data/genomeAssembly//GCspike/",current_sp,".spike.tipRatio2.txt"),quote = F, row.names = F, col.names = T)
  myJudge_df <- reshape2::melt(myJudge)
  #write.table(x = myJudge_df, file = paste0("mid_data/基因组组装/GCspike/",current_sp,".spike.judge.txt"),quote = F, col.names = T)



# data$value %>% fivenum()
# my_min_value = min(data$value)
  
  # p1 <- ggplot(data = data ,
  #              mapping = aes(x = start/1e6, y = value))+
  #   geom_ribbon(aes(ymin = min(value), ymax = value), fill = "#86b2d8")+
  #   geom_segment(data = cp_out_df,aes(x = X1, xend = XE,
  #                                     y = Y1, yend = YE),
  #                colour = "red", linewidth = 1)+
  #   facet_wrap(~chrom, scales = "free") +
  #   labs(x= "Position(Mb)",y="GC%")
  
  #ggsave(plot = p1,width = 12, height = 10,filename = paste0("plots/基因组组装/GC_spike/",current_sp,".GC_spike_min.pdf"))
  #ggsave(plot = p1,width = 12, height = 10,filename = paste0("plots/基因组组装/GC_spike/",current_sp,".GC_spike_min.jpg"))

}


