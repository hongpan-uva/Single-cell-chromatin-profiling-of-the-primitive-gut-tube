library(limma)

setwd("/scratch/hz9fq/toronto_project/E185_Cdx2KO/raw_data")

gMSmtx <- read.table("Symbol_logExp_norm_matrix.txt",header=T)

gMS.design<-cbind(Cdx2KO=c(1,1,1,0,0,0),WT=c(0,0,0,1,1,1))
contrast.matrix <- makeContrasts(Cdx2KO - WT, levels = gMS.design)

fit <- lmFit(gMSmtx, design=gMS.design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

filteredDEG<-topTable(fit2, coef = 1 ,adjust.method = "BH" ,p.value = 0.05, lfc = log2(1.5), number = Inf,sort.by="logFC")
upDEG<-filteredDEG[which(filteredDEG$t>0),]
head(upDEG,n=10)
nrow(upDEG)

downDEG<-filteredDEG[which(filteredDEG$t<0),]
head(downDEG,n=10)
nrow(downDEG)

write.table(upDEG,"../limma_DEG/Cdx2KOvsWT_DEG_table.txt",col.names=T,row.names=T,sep="\t",quote=F)
genenames <- rownames(upDEG)
write.table(genenames,"../limma_DEG/Cdx2KOvsWT_DEG_genelist.txt",col.names=F,row.names=F,sep="\t",quote=F)

write.table(downDEG,"../limma_DEG/WTvsCdx2KO_DEG_table.txt",col.names=T,row.names=T,sep="\t",quote=F)
genenames <- rownames(downDEG)
write.table(genenames,"../limma_DEG/WTvsCdx2KO_DEG_genelist.txt",col.names=F,row.names=F,sep="\t",quote=F)





