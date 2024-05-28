library(M3Drop)
setwd("../rice/")
exp_data <- read.table("TPM.txt",header = T)
HVG <- BrenneckeGetVariableGenes(expr_mat=exp_data, spikes=NA, suppress.plot=FALSE, fdr=0.1, minBiolDisp=0.5, fitMeanQuantile=0.4) # 3997
write.table(HVG,file = "HVG.txt",sep="\t",quote = F,row.names = F,col.names = T)
#HVG <- read.table("HVG.txt",header = T)
HVG_tpm <- exp_data[HVG$Gene,]
HVG_tpm <- t(HVG_tpm)
write.table(HVG_tpm,file = "HVG.tpm.txt",row.names = T,col.names = T,quote = F,sep = "\t")


###############################################################################
## multiple class （wechat）
###############################################################################


# tao tissue specific calculation
HVG_tpm <- read.table("HVG.tpm.txt") 
HVG_tpm <- as.data.frame(t(HVG_tpm))
tao <- apply(HVG_tpm,1,function(x)sum((1-x/max(x))/(ncol(HVG_tpm)-1)))


exp_data <- read.table("TPM.txt",header = T)
tao_all <- apply(exp_data,1,function(x)sum((1-x/max(x))/(ncol(exp_data)-1)))

tao_nonHVG <- tao_all[!names(tao_all)%in%names(tao)]

d <- rbind(data.frame(gene="HVG",tao=tao),data.frame(gene="All",tao=tao_all),data.frame(gene="Non HVG",tao=tao_nonHVG))

ggplot(d,aes(x=gene,y=tao,fill=gene))+
  geom_boxplot()+
  #geom_jitter()+
  geom_signif(comparisons = list(c("HVG","Non HVG"),
                                 c("HVG","All")),
              test="wilcox.test",
              test.args=list(alternative="greater"),
              step_increase = 0.05,tip_length = 0.01)+
  theme_bw()+
  scale_fill_npg()


library(ggplot2)
library(ggsignif)
library(EnvStats)
library(ggsci)
d$gene<-factor(d$gene,levels = c("HVG","Non HVG","All"))
p<-ggplot(d,aes(x=gene,y=tao,col=gene))+
  geom_violin(trim = FALSE,aes(fill=gene),alpha=0.5,col="white",scale = "width")+
  geom_boxplot(width=0.1,outliers = F)+
  scale_fill_aaas()+
  scale_color_aaas()+
  geom_signif(comparisons = list(c("HVG","Non HVG"),
                                 c("HVG","All")),
              test="wilcox.test", test.args=list(alternative="greater"),
              step_increase = 0.05,tip_length = 0.01)+
  theme_bw(base_size = 20)+
  #scale_x_discrete(labels=c("Positive","Negative"))+
  theme(legend.position="none",axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0.5,size=20,face = "bold"),
        theme(legend.position = "None",text = element_text(size = 8),
              axis.title =element_text(size=10)))+ # title posistion
  ylab("Tao index")+
  stat_n_text()+
  ylim(c(0.8,NA))

ggsave(p,file="result_pictures/rice_tao.pdf", width = 105, height = 105, units = "mm")



### GO enrichment
library(gprofiler2)
go_result <-gost(query = rownames(HVG_tpm),organism = "oindica",significant = T,highlight = T)
go_res <- go_result$result
go_res <- subset(go_res,go_res$highlighted==TRUE)

mytheme <- theme(
  axis.title = element_text(size = 9),
  axis.text = element_text(size = 11),
  plot.title = element_text(size = 14,
                            hjust = 0.5,
                            face = "bold"),
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 11),
  plot.margin = margin(t = 5.5,
                       r = 10,
                       l = 5.5,
                       b = 5.5)
)
mytheme2 <- mytheme + theme(axis.text.y = element_blank())




p1 <- ggplot(data = go_res, aes(y = reorder(term_name,-p_value),x=intersection_size,fill = -log10(p_value))) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  labs(x = "Number of genes", y = "GO terms") +
  geom_text(aes(x = 0.03, #用数值向量控制文本标签起始位置
                label = term_name),
            hjust = 0)+ #hjust = 0,左对齐
  theme_classic() + 
  mytheme2+
  facet_wrap(source~.,scales = "free_y",nrow = 3,strip.position = "left")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),text=element_text(size = 2))
