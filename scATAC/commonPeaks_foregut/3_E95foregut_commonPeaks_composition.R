setwd("/scratch/hz9fq/toronto_project/scATAC/ArchR_workspace/GutProject_foregut")

E95common <- read.table("GutProject_foregut_common_peakset.bed", header=F)

E13dir = "/scratch/hz9fq/toronto_project/E13bulkATACpeak/narrowPeak/"
E13file = c("GSM394318_E13_Stomach_ATAC_r1_m2_peaks.narrowPeak","E135_lung_atac_m2_peaks.narrowPeak")

resVec <- c()

for(i in 1:2){
    E13peak <- read.table(paste0(E13dir,E13file[i]), header=T)

    library(GenomicRanges)
    gr1 <- GRanges(
        seqnames = factor(E95common[,1]),
        ranges = IRanges(start = E95common[,2],end = E95common[,3])
    )
    gr2 <- GRanges(
        seqnames = factor(E13peak[,1]),
        ranges = IRanges(start = E13peak[,2],end = E13peak[,3])
    )

    res <- countOverlaps(gr1, gr2)
    resVec <- c(resVec, res)

}

resmtx <- matrix(resVec, ncol=2)

colnames(resmtx) <- c("E13stomach", "E13lung")

write.table(resmtx,"E95foregut_common_peak_association_matrix.txt",row.names=T,col.names=T,quote=F,sep="\t")

resmtx <- read.table("E95foregut_common_peak_association_matrix.txt",header=T)

E95_lost <- E95common[which(resmtx[,1]==0 & resmtx[,2]==0),]
write.table(E95_lost,"E95_lostpeaks.bed",col.names=F,row.names=F,sep="\t",quote=F)

E95_lost <- E95common[which(resmtx[,1]> 0& resmtx[,2]==0),]
write.table(E95_lost,"E95_lung_uniquely_close.bed",col.names=F,row.names=F,sep="\t",quote=F)

E95_lost <- E95common[which(resmtx[,1]== 0 & resmtx[,2]> 0),]
write.table(E95_lost,"E95_stomach_uniquely_close.bed",col.names=F,row.names=F,sep="\t",quote=F)

E95_lost <- E95common[which(resmtx[,1]> 0 & resmtx[,2]> 0),]
write.table(E95_lost,"E95_maintainedpeaks.bed",col.names=F,row.names=F,sep="\t",quote=F)

barVec <- c(length(which(resmtx[,1]> 0 & resmtx[,2]> 0)),
            length(which(resmtx[,1]==0 & resmtx[,2]==0)),
            length(which(resmtx[,1]> 0 & resmtx[,2]==0)),
            length(which(resmtx[,1]==0 & resmtx[,2]> 0))
            )

bardf <- data.frame(data=barVec)
data_percentage <- apply(bardf, 2, function(x){x*100/sum(x,na.rm=T)})

pal <- c("#fc8d62","#66c2a5","#8da0cb","#e78ac3")

pdf("E95foregut_common_peak_association_stacked_barplot.pdf")
barplot(data_percentage, col=pal , border=pal, xlab="group",xlim=c(0,6),main="",axes=F)
axis(2,at=seq(0,100,length.out=6),labels=c("0%","20%","40%","60%","80%","100%"),las=1)
legend("right",legend=c("closed in stomach uniquely","closed in lung uniquely","closed in stomach and lung","maintained in stomach and lung"),
fill = rev(pal))
text(0.7,(data_percentage[1]/2),barVec[1])
text(0.7,data_percentage[1]+(data_percentage[2]/2),barVec[2])
text(0.7,data_percentage[1]+data_percentage[2]+(data_percentage[3]/2),barVec[3])
text(0.7,data_percentage[1]+data_percentage[2]+data_percentage[3]+(data_percentage[4]/2),barVec[4])
dev.off()

#get repressed genes associated with closed peaks (not used)
setwd("/scratch/hz9fq/toronto_project/scATAC/ArchR_workspace/GutProject_foregut")
genelist <- read.table("/scratch/hz9fq/toronto_project/E95E135_CrossStages/RNA_difgene/difgene_lung_down.genelist.txt",header=F)

genelist <- genelist$V1

t0 <- read.table("/nv/vol190/zanglab/hz9fq/annotation/Mouse/mm10_gene_annotation_geneID_LenOrder_TSS50kb.bed",header=F)
t1 <- t0[which(t0$V5 %in% genelist),]

t2 <- read.table("E95_stomachpeaks.bed",header=F)

library(GenomicRanges)
gr1 <- GRanges(
    seqnames = factor(t1[,1]),
    ranges = IRanges(start = t1[,2],end = t1[,3])
)
gr2 <- GRanges(
    seqnames = factor(t2[,1]),
    ranges = IRanges(start = t2[,2],end = t2[,3])
)

res<-countOverlaps(gr1, gr2)

rep_gene <- genelist[which(res>0)]

write.table(rep_gene, "lung_repressed_gene_associatedwith_stomachUniquePeaks.txt",row.names=F,col.names=F,sep="\t",quote=F)

#repressed gene composition (not used)
genelist <- read.table("/scratch/hz9fq/toronto_project/E95E135_CrossStages/RNA_difgene/difgene_lung_down.genelist.txt",header=F)
genelist <- genelist$V1
gene1 <- read.table("lung_repressed_gene_associatedwith_stomachUniquePeaks.txt",header=F)
gene1 <- gene1$V1
gene2 <- read.table("lung_repressed_gene_associatedwith_lostPeaks.txt",header=F)
gene2 <- gene2$V1

intersectN <- length(intersect(gene1,gene2))

barVec <- rev(c(length(gene1)-intersectN,
            intersectN,
            length(gene2)-intersectN,
            length(genelist)-length(gene1)-length(gene2)+intersectN))
bardf <- data.frame(data=barVec)
data_percentage <- apply(bardf, 2, function(x){x*100/sum(x,na.rm=T)})

pal <- c("#9e9e9e","#66c2a5","#607d8b","#8da0cb")

pdf("E135_lung_repressed_genes_stacked_barplot.pdf")
barplot(data_percentage, col=pal , border=pal, xlab="group",xlim=c(0,6),main="")
legend("right",legend=c("type1 repressed genes", "intersection between type1 and type2", "type2 repressed genes", "type3 repressed genes"),
fill = rev(pal))
text(0.7,(data_percentage[1]/2),barVec[1])
text(0.7,data_percentage[1]+(data_percentage[2]/2),barVec[2])
text(0.7,data_percentage[1]+data_percentage[2]+(data_percentage[3]/2),barVec[3])
text(0.7,data_percentage[1]+data_percentage[2]+data_percentage[3]+(data_percentage[4]/2),barVec[4])
dev.off()


