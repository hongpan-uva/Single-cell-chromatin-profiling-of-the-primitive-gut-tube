################E9.5/E13.5 expression barplot (compare within organ)
setwd("/scratch/hz9fq/toronto_project/E95_scRNA")
t1 <- read.table("E95_FPKM_matrix_4Organs.txt",header=T)

t0<-read.table("/scratch/hz9fq/toronto_project/ALL_RNA_collection/mod_exp/FPKM_collection.txt",header=T)
t2 <- data.frame(intestine=apply(t0[,c(7,8)],1,mean),pancreas=apply(t0[,c(9,10)],1,mean),stomach=apply(t0[,c(11,12)],1,mean),lung=apply(t0[,c(13,14)],1,mean))
rownames(t2) <- t0$Gene.ID

pairs = list(c("pancreas","Mafb"),c("intestine","Hnf4g"),c("stomach","Grhl3"),c("lung","Hopx"))
pal=c("#D3C035","#D51F26","#C06CAB","#89C75F")
ticks = list(seq(0,12,2),seq(0,20,5),seq(0,35,5),seq(0,10,2))

pdf("fig2_genetrack_barplot.pdf")
for(pr in 1:length(pairs)){
    organ=pairs[[pr]][1]
    gene=pairs[[pr]][2]

    plotdat <- c(t1[gene,organ],t2[gene,organ])
    names(plotdat) <- c(paste0("E9.5_",organ),paste0("E13.5_",organ))
    print(max(plotdat))

    barplot(plotdat,space=1,col=pal[pr],axes=F,ylim=c(0,max(ticks[[pr]])))
    axis(2,at=ticks[[pr]],labels=ticks[[pr]])
}
dev.off()

################E9.5/E13.5 expression barplot (compare all organs)
library(magrittr)

setwd("/scratch/hz9fq/toronto_project/E95_scRNA")
t1 <- read.table("E95_FPKM_matrix_4Organs.txt",header=T)

t0<-read.table("/scratch/hz9fq/toronto_project/ALL_RNA_collection/mod_exp/FPKM_collection.txt",header=T)
t2 <- data.frame(intestine=apply(t0[,c(7,8)],1,mean),pancreas=apply(t0[,c(9,10)],1,mean),stomach=apply(t0[,c(11,12)],1,mean),lung=apply(t0[,c(13,14)],1,mean))
rownames(t2) <- t0$Gene.ID

geneVec=c("Mafb","Hnf4g","Grhl3","Hopx")
pal=c("#D3C035","#D51F26","#C06CAB","#89C75F")
pal=c("#D3C035","#D3C035","#D51F26","#D51F26","#C06CAB","#C06CAB","#89C75F","#89C75F")
ticks = list(seq(0,12,2),seq(0,20,5),seq(0,35,5),seq(0,40,5))

#functions to shade bar plot and legend
pdf("fig2_genetrack_barplot_complete.pdf")
for(i in 1:length(geneVec)){
    gene=geneVec[i]

    c(t1[gene,"intestine"],t2[gene,"intestine"],
      t1[gene,"pancreas"],t2[gene,"pancreas"],
      t1[gene,"stomach"],t2[gene,"stomach"],
      t1[gene,"lung"],t2[gene,"lung"]) %>%
    matrix(byrow=F,nrow=2) -> plotmtx

    row.names(plotmtx) <- c("E9.5","E13.5")
    colnames(plotmtx) <- c("intestine","pancreas","stomach","lung")

    barplot(main=paste(gene,"FPKM"), ylim=c(0,max(ticks[[i]])), height = plotmtx, beside = TRUE, legend=TRUE)
    axis(2,at=ticks[[i]],labels=ticks[[i]])
}
dev.off()
