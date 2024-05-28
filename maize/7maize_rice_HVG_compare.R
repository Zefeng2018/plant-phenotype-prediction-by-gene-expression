library(stringr)
setwd("/home/wuzefeng/MyResearch/Tissue_prediction/maize")
orthos <- read.table("../orthologs/maize2rice.orthologs.txt",header = T)
orthos$subject_id <- str_replace(string = orthos$subject_id,pattern ="t" ,replacement = "g") 


maize_hvg <- read.table("HVG.txt",header = T)
rice_hvg <- read.table("../rice/HVG.txt",header = F)

orthos$maize.hvg <- ifelse(orthos$query_id%in%maize_hvg$Gene,1,0) 
orthos$os.hvg <- ifelse(orthos$subject_id%in%rice_hvg$V1,1,0) 


common_hvgs <- table(orthos[,c("maize.hvg","os.hvg")])[2,2] # 486 

### random simulation
commons <- c()
maize_genes <- read.table("../maize/maize.genes.id.txt")
rice_genes <- read.table("rice.genes.id.txt")
for (m in seq(1:1000)){
  maize_hvg_sim <- sample(maize_genes$V1,2880,replace = F)
  rice_hvg_sim <- sample(rice_genes$V1,3997,replace = F)
  orthos$maize.hvg <- ifelse(orthos$query_id%in%maize_hvg_sim,1,0) 
  orthos$os.hvg <- ifelse(orthos$subject_id%in%rice_hvg_sim,1,0) 
  common_hvgs <- table(orthos[,c("maize.hvg","os.hvg")])[2,2]
  message(m,common_hvgs)
  commons <- c(commons,common_hvgs)
}

library(ggpubr)
library(tidyverse)
library(ggbreak)

commons %>% as.data.frame(nm = "Number") %>% gghistogram(x="Number",bins = 100,col="steelblue")+theme_classic2(base_size = 16)+
  annotate("segment", x = 486, xend = 486, y = 50, 
           yend = 5, colour = "orange", size = 1, arrow = arrow())+
  scale_x_break(c(150,470),#截断位置及范围
                space = 0.2,#间距大小
                scales = 0.4)
