setwd("/scratch/hz9fq/toronto_project/E185_Cdx2KO/limma_DEG")

genelist<- read.table("Cdx2KOvsWT_DEG_genelist_confirmedHGNCsymbol.txt",header=F)
genelist <- genelist$V1
anno <- read.table("/nv/vol190/zanglab/hz9fq/annotation/Mouse/mm10_gene_annotation_geneID_LenOrder_TSS50kb.bed",header=F)
genetable <- anno[which(anno$V5%in%genelist),]

Peak_dir <- c("/scratch/hz9fq/toronto_project/Sox2Cdx2_peak/signal_Cdx2/")

resVec <- c()

for(pk in c("Peak_group1","Peak_group2","Peak_group3")){
    peakset <- read.table(paste0(Peak_dir,pk,".bed"),header=F)

    library(GenomicRanges)
    gr1 <- GRanges(
        seqnames = factor(genetable[,1]),
        ranges = IRanges(start = genetable[,2],end = genetable[,3])
    )
    gr2 <- GRanges(
        seqnames = factor(peakset[,1]),
        ranges = IRanges(start = peakset[,2],end = peakset[,3])
    )
    res<-countOverlaps(gr1, gr2)
    
    resVec <- c(resVec,res)
}

resmtx <- matrix(resVec,ncol=3)
rownames(resmtx) <- genelist
colnames(resmtx) <- c("Peak_group1(gained)","Peak_group2(maintained)","Peak_group3(lost)")

#fisher exact test (gene-centered)
#              target peak   other peak
#target gene       X11           X12
#other gene        X21           X22

othergene <- anno[-which(anno$V5 %in% genetable$V5),]

resVec2 <- c()

for(pk in c("Peak_group1","Peak_group2","Peak_group3")){
    peakset <- read.table(paste0(Peak_dir,pk,".bed"),header=F)

    library(GenomicRanges)
    gr1 <- GRanges(
        seqnames = factor(othergene[,1]),
        ranges = IRanges(start = othergene[,2],end = othergene[,3])
    )
    gr2 <- GRanges(
        seqnames = factor(peakset[,1]),
        ranges = IRanges(start = peakset[,2],end = peakset[,3])
    )
    res<-countOverlaps(gr1, gr2)
    
    resVec2 <- c(resVec2,res)
}

resmtx2 <- matrix(resVec2,ncol=3)
rownames(resmtx2) <- othergene$V5
colnames(resmtx2) <- c("Peak_group1(gained)","Peak_group2(maintained)","Peak_group3(lost)")

X11 <- length(which(resmtx[,1]>0))
X12 <- length(which(resmtx[,1]==0 & (resmtx[,2]>0 | resmtx[,3]>0) ))
X21 <- length(which(resmtx2[,1]>0))
X22 <- length(which(resmtx2[,1]==0 & (resmtx2[,2]>0 | resmtx2[,3]>0) ))
Contimtx <- matrix(c(X11,X21,X12,X22),nrow=2)
fishres <- fisher.test(Contimtx,alternative = "greater")
print(fishres)
print(fishres$p.value)

X11 <- length(which(resmtx[,2]>0))
X12 <- length(which(resmtx[,2]==0 & (resmtx[,1]>0 | resmtx[,3]>0) ))
X21 <- length(which(resmtx2[,2]>0))
X22 <- length(which(resmtx2[,2]==0 & (resmtx2[,1]>0 | resmtx2[,3]>0) ))
Contimtx <- matrix(c(X11,X21,X12,X22),nrow=2)
fishres <- fisher.test(Contimtx,alternative = "two.sided")
print(fishres)
print(fishres$p.value)

X11 <- length(which(resmtx[,3]>0))
X12 <- length(which(resmtx[,3]==0 & (resmtx[,1]>0 | resmtx[,2]>0) ))
X21 <- length(which(resmtx2[,3]>0))
X22 <- length(which(resmtx2[,3]==0 & (resmtx2[,1]>0 | resmtx2[,2]>0) ))
Contimtx <- matrix(c(X11,X21,X12,X22),nrow=2)
fishres <- fisher.test(Contimtx,alternative = "less")
print(fishres)
print(fishres$p.value)


#matrix output
write.table(resmtx,"Cdx2KOvsWT_peak_association_matrix.txt",row.names=T,col.names=T,quote=F,sep="\t")

resmtx <- read.table("Cdx2KOvsWT_peak_association_matrix.txt",header=T)

numericmtx <- resmtx
row.names(numericmtx) <- 1:nrow(numericmtx)

#peak association euler plot
s4 <- list(gained = as.numeric(rownames(numericmtx)[which(numericmtx[,1]>0)]),
        maintained = as.numeric(rownames(numericmtx)[which(numericmtx[,2]>0)]),
        lost = as.numeric(rownames(numericmtx)[which(numericmtx[,3]>0)]),
        no_peak = as.numeric(rownames(numericmtx)[which(numericmtx[,1]+numericmtx[,2]+numericmtx[,3]==0)])
        )

library(eulerr)

pdf("Cdx2KOvsWT_peak_association_euler.pdf")
plot(euler(s4,shape="circle"), quantities = TRUE)
dev.off()

#venn plot
library(VennDiagram)

pdf(paste0("Cdx2KOvsWT_peak_association_venn.pdf"))
draw.triple.venn(area1=length(which(numericmtx[,1]>0)),
                 area2=length(which(numericmtx[,2]>0)),
                 area3=length(which(numericmtx[,3]>0)),
                 n12=length(which(numericmtx[,1]>0 & numericmtx[,2]>0)),
                 n23=length(which(numericmtx[,2]>0 & numericmtx[,3]>0)),
                 n13=length(which(numericmtx[,1]>0 & numericmtx[,3]>0)),
                 n123=length(which(numericmtx[,1]>0 & numericmtx[,2]>0 & numericmtx[,3]>0)),
                 fill=c("#fdd2ce","#d1e1ee","#bbe1cc"),
                 lwd=c(1,1,1),cex=2,cat.cex=2,cat.default.pos="outer",cat.pos=0,rotation.degree=120)

dev.off()


#UpSet diagram
onehotmtx <- numericmtx
for(i in 1:3){
    onehotmtx[which(onehotmtx[,i] > 0),i] <- 1
}

colnames(onehotmtx) <- c("genes_with_gained_peaks","genes_with_maintained_peaks","genes_with_lost_peaks")

library(UpSetR)

pdf("Cdx2KOvsWT_peak_association_UpSet.pdf")
upset(onehotmtx, nsets = 3, nintersects = 20, sets = c("genes_with_lost_peaks","genes_with_maintained_peaks","genes_with_gained_peaks"), 
      keep.order=TRUE , mb.ratio = c(0.5, 0.5),
      order.by = c("degree"), decreasing = c(FALSE),
      main.bar.color = "#548CA8", sets.bar.color = "#476072", text.scale=c(1.5,1.5,1.5,1.5,1.5,1.5))
dev.off()

#get interested gene list
genelist_open  <- rownames(resmtx)[which(resmtx[,1]>0 & resmtx[,2]==0 & resmtx[,3]==0)]
genelist_close <- rownames(resmtx)[which(resmtx[,1]==0 & resmtx[,2]==0 & resmtx[,3]>0)]

write.table(genelist_open,"Cdx2KOvsWT_gene_associatedwith_openRegion_uniquely.genelist.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(genelist_close,"Cdx2KOvsWT_gene_associatedwith_closeRegion_uniquely.genelist.txt",row.names=F,col.names=F,quote=F,sep="\t")

#get interested gene list
genelist_open  <- rownames(resmtx)[which(resmtx[,1]>0)]
genelist_close <- rownames(resmtx)[which(resmtx[,3]>0)]
genelist_maintained <- rownames(resmtx)[which(resmtx[,2]>0)]

write.table(genelist_open,"Cdx2KOvsWT_gene_associatedwith_openRegion.genelist.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(genelist_close,"Cdx2KOvsWT_gene_associatedwith_closeRegion.genelist.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(genelist_maintained,"Cdx2KOvsWT_gene_associatedwith_maintainedRegion.genelist.txt",row.names=F,col.names=F,quote=F,sep="\t")


