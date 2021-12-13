#get down-regulated genes
setwd("/scratch/hz9fq/toronto_project/E95E135_CrossStages/RNA_difgene")

#read in FPKM
E95.mtx <- read.table("E95_FPKM_matrix_4Organs.txt",header=T,row.names=1)
E135.mtx <- read.table("E135_FPKM_matrix_4Organs.txt",header=T,row.names=1)

gene.use <- intersect(rownames(E95.mtx),rownames(E135.mtx))
E95.mtx.use <- E95.mtx[gene.use,]
colnames(E95.mtx.use) <- paste0("E95_",colnames(E95.mtx.use))
E135.mtx.use <- E135.mtx[gene.use,]
colnames(E135.mtx.use) <- paste0("E135_",colnames(E135.mtx.use))

organs=c("stomach","lung")

all.mtx <- cbind(E95.mtx.use,E135.mtx.use)
all.mtx.norm <- preprocessCore::normalize.quantiles(as.matrix(all.mtx))
all.mtx.norm <- as.data.frame(all.mtx.norm)
colnames(all.mtx.norm) <- colnames(all.mtx)
rownames(all.mtx.norm) <- rownames(all.mtx)

#get log2 normalized fpkm
all.mtx.log <- log2(all.mtx.norm+1)

#get dif genes
fore <- all.mtx.log[,paste0("E135_","stomach")]
back <- all.mtx.log[,paste0("E95_","stomach")]
fc_stomach <- (fore+0.5)/(back+0.5)

fore <- all.mtx.log[,paste0("E135_","lung")]
back <- all.mtx.log[,paste0("E95_","lung")]
fc_lung <- (fore+0.5)/(back+0.5)

foldchange.df <- data.frame(symbol=rownames(all.mtx.log),fc_stomach=fc_stomach,fc_lung=fc_lung)

foldchange.df <- foldchange.df[order(foldchange.df$fc_stomach),]
stomach_down<- foldchange.df[which(foldchange.df$fc_stomach<0.5),]
stomach_down$E9.5_stomach_expression <- all.mtx.log[stomach_down$symbol,"E95_stomach"]
stomach_down$E13.5_stomach_expression <- all.mtx.log[stomach_down$symbol,"E135_stomach"]
stomach_down$E9.5_lung_expression <- all.mtx.log[stomach_down$symbol,"E95_lung"]
stomach_down$E13.5_lung_expression <- all.mtx.log[stomach_down$symbol,"E135_lung"]
stomach_down <- stomach_down[,c(1,4,5,2,6,7,3)]
write.table(stomach_down,paste0("difgene_stomach_down_for_foregut_analysis.table.txt"),row.names=F,col.names=T,quote=F,sep="\t")

foldchange.df <- foldchange.df[order(foldchange.df$fc_lung),]
lung_down<- foldchange.df[which(foldchange.df$fc_lung<0.5),]
lung_down$E9.5_stomach_expression <- all.mtx.log[lung_down$symbol,"E95_stomach"]
lung_down$E13.5_stomach_expression <- all.mtx.log[lung_down$symbol,"E135_stomach"]
lung_down$E9.5_lung_expression <- all.mtx.log[lung_down$symbol,"E95_lung"]
lung_down$E13.5_lung_expression <- all.mtx.log[lung_down$symbol,"E135_lung"]
lung_down <- lung_down[,c(1,4,5,2,6,7,3)]
write.table(lung_down,paste0("difgene_lung_down_for_foregut_analysis.table.txt"),row.names=F,col.names=T,quote=F,sep="\t")


setwd("/scratch/hz9fq/toronto_project/scATAC/ArchR_workspace/GutProject_foregut/E95RNA/byFC")

#read in gene tables
E135stom_tab <- read.table("/scratch/hz9fq/toronto_project/E95E135_CrossStages/RNA_difgene/difgene_stomach_down_for_foregut_analysis.table.txt",header=T)
rownames(E135stom_tab) <- E135stom_tab$symbol
E135stom <- E135stom_tab$symbol
E135lung_tab <- read.table("/scratch/hz9fq/toronto_project/E95E135_CrossStages/RNA_difgene/difgene_lung_down_for_foregut_analysis.table.txt",header=T)
rownames(E135lung_tab) <- E135lung_tab$symbol
E135lung <- E135lung_tab$symbol

E95common <- unique(c(E135stom,E135lung))

#composition table
compo_mtx <- matrix(c(as.numeric(E95common%in%E135stom),as.numeric(E95common%in%E135lung)),byrow=F,ncol=2)
colnames(compo_mtx) <- c("E135_stom","E135_lung")
rownames(compo_mtx) <- E95common
compo_df <- as.data.frame(compo_mtx)

#!! 1 indicates repression, 0 indicates non-represssion
lost <- rownames(compo_df)[which(compo_df[,1]==1 & compo_df[,2]==1)]
write.table(lost,"E95_lost_genelist.txt",col.names=F,row.names=F,sep="\t",quote=F)

lungRep <- rownames(compo_df)[which(compo_df[,1]== 0& compo_df[,2]==1)]
write.table(lungRep,"E95_lung_uniquely_repressed_genelist.txt",col.names=F,row.names=F,sep="\t",quote=F)

stomRep <- rownames(compo_df)[which(compo_df[,1]== 1 & compo_df[,2]== 0)]
write.table(stomRep,"E95_stomach_uniquely_repressed_genelist.txt",col.names=F,row.names=F,sep="\t",quote=F)

barVec <- c(length(which(compo_df[,1]==1 & compo_df[,2]==1)),
            length(which(compo_df[,1]==0 & compo_df[,2]==1)),
            length(which(compo_df[,1]==1 & compo_df[,2]==0))
            )

bardf <- data.frame(data=barVec)
data_percentage <- apply(bardf, 2, function(x){x*100/sum(x,na.rm=T)})

pal <- c("#66c2a5","#8da0cb","#e78ac3")

#used plot
pdf("E95foregut_common_gene_association_stacked_barplot.pdf")
barplot(data_percentage, col=pal , border=pal, xlab="group",xlim=c(0,6),main="",axes=F)
axis(2,at=seq(0,100,length.out=6),labels=c("0%","20%","40%","60%","80%","100%"),las=1)
legend("right",legend=c("repressed in stomach uniquely","repressed in lung uniquely","repressed in stomach and lung"),
fill = rev(pal))
text(0.7,(data_percentage[1]/2),barVec[1])
text(0.7,data_percentage[1]+(data_percentage[2]/2),barVec[2])
text(0.7,data_percentage[1]+data_percentage[2]+(data_percentage[3]/2),barVec[3])
dev.off()

#refined by closed peaks
t0 <- read.table("/nv/vol190/zanglab/hz9fq/annotation/Mouse/mm10_gene_annotation_geneID_LenOrder_TSS50kb.bed",header=F)
rownames(t0) <- t0$V5
lost[-which(lost%in%rownames(t0))]

#both repressed
t1 <- t0[lost,1:3]

t2 <- read.table("../../E95_lostpeaks.bed",header=F)

library(GenomicRanges)
gr1 <- GRanges(
    seqnames = factor(t1[,1]),
    ranges = IRanges(start = as.numeric(t1[,2]),end = as.numeric(t1[,3]))
)
gr2 <- GRanges(
    seqnames = factor(t2[,1]),
    ranges = IRanges(start = as.numeric(t2[,2]),end = as.numeric(t2[,3]))
)

res1 <- countOverlaps(gr1, gr2)
hit1 <- findOverlaps(gr1, gr2)

refined_gID<- which(res1>0)
hit1.df<-as.data.frame(hit1)

associated_pString.vec <- c()
for(i in refined_gID){
    associated_pID <- hit1.df$subjectHits[which(hit1.df$queryHits==i)]
    associated_pString <- ""
    for(j in associated_pID){
        peakString<- paste0(t2[j,1],":",t2[j,2],"-",t2[j,3])
        if(associated_pString==""){
            associated_pString <- peakString
        }else{
            associated_pString <- paste0(associated_pString,",",peakString)
        }
    }
    associated_pString.vec <- c(associated_pString.vec,associated_pString)
}

refined_gSymbol <- rownames(t1)[refined_gID]
lost_refined.df<- data.frame(Symbol=refined_gSymbol,
           E9.5_stomach_expression=E135stom_tab[refined_gSymbol,"E9.5_stomach_expression"],
           E13.5_stomach_expression=E135stom_tab[refined_gSymbol,"E13.5_stomach_expression"],
           FC_stomach=E135stom_tab[refined_gSymbol,"fc_stomach"],
           E9.5_lung_expression=E135lung_tab[refined_gSymbol,"E9.5_lung_expression"],
           E13.5_lung_expression=E135lung_tab[refined_gSymbol,"E13.5_lung_expression"],
           FC_lung=E135lung_tab[refined_gSymbol,"fc_lung"],
           peaks_associated=associated_pString.vec)

write.table(lost_refined.df, "E95_lost_refined_gene_table.txt", row.names=F, col.names=T, sep="\t", quote=F)

#lung repressed
t1 <- t0[lungRep,1:3]
t2 <- read.table("../../E95_lung_uniquely_close.bed",header=F)

library(GenomicRanges)
gr1 <- GRanges(
    seqnames = factor(t1[,1]),
    ranges = IRanges(start = as.numeric(t1[,2]),end = as.numeric(t1[,3]))
)
gr2 <- GRanges(
    seqnames = factor(t2[,1]),
    ranges = IRanges(start = as.numeric(t2[,2]),end = as.numeric(t2[,3]))
)

res2 <- countOverlaps(gr1, gr2)
hit2 <- findOverlaps(gr1, gr2)

refined_gID<- which(res2>0)
hit2.df<-as.data.frame(hit2)

associated_pString.vec <- c()
for(i in refined_gID){
    associated_pID <- hit2.df$subjectHits[which(hit2.df$queryHits==i)]
    associated_pString <- ""
    for(j in associated_pID){
        peakString<- paste0(t2[j,1],":",t2[j,2],"-",t2[j,3])
        if(associated_pString==""){
            associated_pString <- peakString
        }else{
            associated_pString <- paste0(associated_pString,",",peakString)
        }
    }
    associated_pString.vec <- c(associated_pString.vec,associated_pString)
}

refined_gSymbol <- rownames(t1)[refined_gID]
lungRep_refined.df<- data.frame(Symbol=refined_gSymbol,
           E9.5_stomach_expression=E135lung_tab[refined_gSymbol,"E9.5_stomach_expression"],
           E13.5_stomach_expression=E135lung_tab[refined_gSymbol,"E13.5_stomach_expression"],
           FC_stomach=E135lung_tab[refined_gSymbol,"fc_stomach"],
           E9.5_lung_expression=E135lung_tab[refined_gSymbol,"E9.5_lung_expression"],
           E13.5_lung_expression=E135lung_tab[refined_gSymbol,"E13.5_lung_expression"],
           FC_lung=E135lung_tab[refined_gSymbol,"fc_lung"],
           peaks_associated=associated_pString.vec)

write.table(lungRep_refined.df, "E95_lung_uniquely_repressed_refined_gene_table.txt", row.names=F, col.names=T, sep="\t", quote=F)

#stomach repressed
t1 <- t0[stomRep,1:3]
t2 <- read.table("../../E95_stomach_uniquely_close.bed",header=F)

library(GenomicRanges)
gr1 <- GRanges(
    seqnames = factor(t1[,1]),
    ranges = IRanges(start = as.numeric(t1[,2]),end = as.numeric(t1[,3]))
)
gr2 <- GRanges(
    seqnames = factor(t2[,1]),
    ranges = IRanges(start = as.numeric(t2[,2]),end = as.numeric(t2[,3]))
)

res3 <- countOverlaps(gr1, gr2)
hit3 <- findOverlaps(gr1, gr2)

refined_gID<- which(res3>0)
hit3.df<-as.data.frame(hit3)

associated_pString.vec <- c()
for(i in refined_gID){
    associated_pID <- hit3.df$subjectHits[which(hit3.df$queryHits==i)]
    associated_pString <- ""
    for(j in associated_pID){
        peakString<- paste0(t2[j,1],":",t2[j,2],"-",t2[j,3])
        if(associated_pString==""){
            associated_pString <- peakString
        }else{
            associated_pString <- paste0(associated_pString,",",peakString)
        }
    }
    associated_pString.vec <- c(associated_pString.vec,associated_pString)
}

refined_gSymbol <- rownames(t1)[refined_gID]
stomRep_refined.df<- data.frame(Symbol=refined_gSymbol,
           E9.5_stomach_expression=E135stom_tab[refined_gSymbol,"E9.5_stomach_expression"],
           E13.5_stomach_expression=E135stom_tab[refined_gSymbol,"E13.5_stomach_expression"],
           FC_stomach=E135stom_tab[refined_gSymbol,"fc_stomach"],
           E9.5_lung_expression=E135stom_tab[refined_gSymbol,"E9.5_lung_expression"],
           E13.5_lung_expression=E135stom_tab[refined_gSymbol,"E13.5_lung_expression"],
           FC_lung=E135stom_tab[refined_gSymbol,"fc_lung"],
           peaks_associated=associated_pString.vec)

write.table(stomRep_refined.df, "E95_stomach_uniquely_repressed_refined_gene_table.txt", row.names=F, col.names=T, sep="\t", quote=F)


barVec <- c(length(which(res1==0)),length(which(res1>0)),
            length(which(res2==0)),length(which(res2>0)),
            length(which(res3==0)),length(which(res3>0))
            )

bardf <- data.frame(data=barVec)
data_percentage <- apply(bardf, 2, function(x){x*100/sum(x,na.rm=T)})

pal <- c("white","#66c2a5","white","#8da0cb","white","#e78ac3")

#used plot
pdf("E95foregut_common_gene_association_refined_stacked_barplot.pdf")
barplot(data_percentage, col=pal , border=pal, xlab="group",xlim=c(0,6),main="",axes=F)
axis(2,at=seq(0,100,length.out=6),labels=c("0%","20%","40%","60%","80%","100%"),las=1)
cum_percent <- cumsum(data_percentage)
text(0.7,cum_percent[1]+(data_percentage[2]/2),barVec[2])
text(0.7,cum_percent[3]+(data_percentage[4]/2),barVec[4])
text(0.7,cum_percent[5]+(data_percentage[6]/2),barVec[6])
dev.off()

#write out refined gene list
lost_refined <- lost[which(res1>0)]
write.table(lost,"E95_lost_unrefined_genelist.txt",col.names=F,row.names=F,sep="\t",quote=F)
write.table(lost_refined,"E95_lost_refined_genelist.txt",col.names=F,row.names=F,sep="\t",quote=F)

lungRep_refined <- lungRep[which(res2>0)]
write.table(lungRep,"E95_lung_uniquely_repressed_unrefined_genelist.txt",col.names=F,row.names=F,sep="\t",quote=F)
write.table(lungRep_refined,"E95_lung_uniquely_repressed_refined_genelist.txt",col.names=F,row.names=F,sep="\t",quote=F)

stomRep_refined <- stomRep[which(res3>0)]
write.table(stomRep,"E95_stomach_uniquely_repressed_unrefined_genelist.txt",col.names=F,row.names=F,sep="\t",quote=F)
write.table(stomRep_refined,"E95_stomach_uniquely_repressed_refined_genelist.txt",col.names=F,row.names=F,sep="\t",quote=F)


