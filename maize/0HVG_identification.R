####################################### HVG selections#####################################################################################
exp_data <- read.table("TPM.txt",header = T)
HVG <- BrenneckeGetVariableGenes(expr_mat=exp_data, spikes=NA, suppress.plot=FALSE, fdr=0.1, minBiolDisp=0.5, fitMeanQuantile=0.8) # 2880
write.table(HVG,file = "HVG.txt",sep="\t",quote = F,row.names = F,col.names = T)
HVG <- read.table("HVG.txt",header = T)
HVG_tpm <- exp_data[HVG$Gene,]
HVG_tpm <- t(HVG_tpm)
write.table(HVG_tpm,file = "HVG.tpm.txt",row.names = T,col.names = T,quote = F,sep = "\t")


###################################### Gene function enrichment analysis ###################################################################
HVG <- read.table("HVG.txt",header = T)
HVG <- HVG$Gene
## GO ontology
maize_go <-  read.table("940_slimGO.txt")
terms <- read.delim("term.txt",sep="\t",quote = "",header = F)

maize_go$onto <- terms$V2[match(maize_go$V4,terms$V1)]
maize_go$term <- terms$V3[match(maize_go$V4,terms$V1)]

CC <- maize_go[maize_go$onto=="C",]
BP <- maize_go[maize_go$onto=="P",]
MF <- maize_go[maize_go$onto=="F",]

term2gene <- CC[,c(4,3)]
term2name <- CC[,c(4,6)]
# 富集分析
library(clusterProfiler)
CC_enrich <- enricher(HVG,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)

term2gene <- MF[,c(4,3)]
term2name <- MF[,c(4,6)]
# 富集分析
MF_enrich <- enricher(HVG,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)

term2gene <- BP[,c(4,3)]
term2name <- BP[,c(4,6)]
# 富集分析
BP_enrich <- enricher(HVG,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)



library(ggpubr)
BP <- dotplot(BP_enrich, title = "BP",label_format=100)
CC <- dotplot(CC_enrich, title = "CC",label_format=100)
MF <- dotplot(MF_enrich, title = "MF",label_format=100)
ggarrange(BP, CC, MF, ncol = 1, nrow = 3, align = "hv")


## elegant  plot

CC_enrich@result$Category <- "CC"
MF_enrich@result$Category <- "MF"
BP_enrich@result$Category <- "BP"

GO_enrich <- rbind(CC_enrich@result[1:10,],MF_enrich@result[1:10,],BP_enrich@result[1:10,])


## plot
mytheme <- theme(
  axis.title = element_text(size = 13),
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
#常规富集条形图绘图：
p <- ggplot(data = GO_enrich, aes(x = Count, y = Description, fill = -log10(pvalue))) +
  scale_fill_distiller(palette = "RdPu",direction = 1) + #更改配色
  geom_bar(stat = "identity", width = 0.8) + #绘制条形图
  labs(x = "Number of Gene", y = "", title = "KEGG enrichment barplot") + #修改/添加各标题
  theme_bw() + mytheme+ #主题更改+
  facet_wrap(.~Category,nrow = 3,scales = "free")


mytheme2 <- mytheme + theme(axis.text.y = element_blank())


p1 <- ggplot(data = GO_enrich, aes(y = reorder(Description,-p.adjust),x=Count,fill = -log10(p.adjust))) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  labs(x = "Number of genes", y = "GO terms") +
  geom_text(aes(x = 0.03, #用数值向量控制文本标签起始位置
                label = Description),
            hjust = 0)+ #hjust = 0,左对齐
  theme_classic() + 
  mytheme2+
  facet_wrap(Category~.,scales = "free_y",nrow = 3,strip.position = "left")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12))
ggsave(filename = "/DATA4T/Brassica_napus_RNA_Seq/3.5tmp_deg/T2toT0_DEG.go.pdf",plot = p1,width=12, height=10)



################################################ HVG expression specificity analysis
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
        plot.title = element_text(hjust = 0.5,size=20,face = "bold"))+ # title posistion
  ylab("Tao index")+
  stat_n_text()+
  ylim(c(0.8,NA))
