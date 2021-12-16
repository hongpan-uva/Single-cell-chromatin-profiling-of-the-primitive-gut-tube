motifoverlaps=paste0("/scratch/hz9fq/toronto_project/E95E135_CrossStages/signal_mergedPeaks_1129/motif_overlap/",c("intestine","pancreas","stomach","lung"),"_peak_motifoverlap.txt")
cor_dir="/scratch/hz9fq/toronto_project/E95E135_CrossStages/signal_mergedPeaks_1129/cor/"
cor_files=c("E95intestine_E135intestine_cor.txt","E95pancreas_E135pancreas_cor.txt","E95stomach_E135stomach_cor.txt","E95lung_E135lung_cor.txt")
plot_dir="/scratch/hz9fq/toronto_project/E95E135_CrossStages/signal_mergedPeaks_1129/cor/plots/"
peaksets=paste0("/scratch/hz9fq/toronto_project/E95E135_CrossStages/mergedPeaks_1129/",c("intestine","pancreas","stomach","lung"),"_mergedPeak.bed")

#prepare FPKM tables
FPKM_E95 = read.table("/scratch/hz9fq/toronto_project/E95_scRNA/E95_FPKM_matrix_4Organs.txt",header=T)
rownames(FPKM_E95) <- toupper(rownames(FPKM_E95))
Exp_E95 <- log2(FPKM_E95+1)
FPKM_col = read.table("/scratch/hz9fq/toronto_project/ALL_RNA_collection/mod_exp/FPKM_collection.txt",header=T)
rownames(FPKM_col) = FPKM_col$Gene.ID
FPKM_E135 <- FPKM_col[,7:14]
FPKM_E135 <- cbind(apply(FPKM_E135[,1:2],1,mean),
                   apply(FPKM_E135[,3:4],1,mean),
                   apply(FPKM_E135[,5:6],1,mean),
                   apply(FPKM_E135[,7:8],1,mean))
colnames(FPKM_E135) <- c("intestine","pancreas","stomach","lung")
rownames(FPKM_E135) <- toupper(rownames(FPKM_E135))
#Exp_E135 <- log2(FPKM_E135+1)

for(set in 1:4){
    motifoverlap <- motifoverlaps[set]
    peakset <- peaksets[set]
    file <- cor_files[set]

    peak.tab <- read.table(peakset,header=F,stringsAsFactors=F)
    motif.tab <- read.table(motifoverlap,header=T,stringsAsFactors=F)
    peakNumber <- apply(motif.tab,2,sum)

    names(peakNumber) <- sapply(names(peakNumber),function(x){
        tmp <- strsplit(x,"_map")
        return(tmp[[1]][1])
    })

    print(file)
    t2 <- read.table(paste0(cor_dir,file),header=T,stringsAsFactors=F)
    rownames(t2) <- t2[,1]

    #set -300 as lower boundary
    t2$wilcox.new[which(t2$wilcox.new<(-300))] <- -300
    t2$wilcox.new[which(t2$wilcox.new>(300))] <- 300

    #wilcox.new
    datatype = "_wilcox.new"
    t2_use <- t2[,c(1,8)]
    colnames(t2_use)[2] <- "Enrichment.Score"
    
    uplimit = max(abs(t2_use[,2]))
    lowlimit = -uplimit

    overlap1 <- intersect(rownames(t2_use),rownames(FPKM_E95))
    overlap2 <- intersect(overlap1,rownames(FPKM_E135))
    t2_use <- t2_use[overlap2,]

    t2_use$peaks<-as.numeric(peakNumber[rownames(t2_use)])
    
    blocks <- seq(min(t2_use$peaks),max(t2_use$peaks),length.out=101)
    blockctrs <- c()
    for(i in 1:(length(blocks)-1)){
        blockctrs <- c(blockctrs,mean(c(blocks[i],blocks[i+1])))
    }
    names(blockctrs) <- 1:100

    getblockctr <- function(x){
     return(names(blockctrs)[which.min(abs(blockctrs-x))])
    }
    t2_use$V4 <- sapply(t2_use$peaks,getblockctr)

    filename <- strsplit(file,"_cor")[[1]][1]
    
    #write.table(t2_use,file=file,col.names=F,row.names=F,sep='\t',quote=F)
    
    horizon_pal=c("#000075","#2E00FF","#9408F7","#C729D6","#FA4AB5","#FF6A95","#FF8B74","#FFAC53","#FFCD32","#FFFF60")
    colors <- colorRampPalette(horizon_pal)(100)
    color.V <- colors[as.numeric(t2_use$V4)]

    t2_use$E95exp <- -1
    count <- 0
    for(i in 1:nrow(t2_use)){
        if(rownames(t2_use)[i]%in%rownames(FPKM_E95)){
            t2_use$E95exp[i] <- FPKM_E95[rownames(t2_use)[i],set]
            count <- count+1
        }
    }
    print(paste0("Number of TF with expression: ",count))
    #t2_use$E95exp <- log2(t2_use$E95exp+1)
    print(summary(t2_use$E95exp))

    t2_use$E135exp <- -1
    count <- 0
    for(i in 1:nrow(t2_use)){
        if(rownames(t2_use)[i]%in%rownames(FPKM_E135)){
            t2_use$E135exp[i] <- FPKM_E135[rownames(t2_use)[i],set]
            count <- count+1
        }
    }
    print(paste0("Number of TF with expression: ",count))
    #t2_use$E135exp <- log2(t2_use$E135exp+1)
    print(summary(t2_use$E135exp))

    #t2_use$rank <- rank(t2_use$wilcox.new,ties.method="first")
    print(max(t2_use[,5]))
    print(max(t2_use[,6]))

    #markers=unique(c(rownames(t2_use)[order(t2_use[,2],decreasing=T)][1:20],rownames(t2_use)[order(t2_use[,2],decreasing=F)][1:20],rownames(t2_use)[order(t2_use[,5],decreasing=T)][1:10],rownames(t2_use)[order(t2_use[,6],decreasing=T)][1:10]))

    #pdf(paste0(plot_dir,filename,datatype,"_byE9Exp.pdf"))
    #par(mar=c(6,6,6,6))
    #plot(t2_use[,2],t2_use[,5],xlab="wilcox_newValue",ylab="E9.5 log2(FPKM+1)",type="p",pch=19,cex=.7,col=color.V,main=paste0(filename,datatype),xlim=c(lowlimit-10,uplimit+10),ylim=c(0,9))
    #labels <- c(min(t2_use$peaks),round((min(t2_use$peaks)+max(t2_use$peaks))/2),max(t2_use$peaks))
    #plotrix::color.legend(360,100,400,250,legend=labels,rect.col=colors,align="rb",gradient="y",cex=0.7)
    #for(i in markers){
    #    text(t2_use[i,2],(t2_use[i,5]+0.2),i,col="black",cex=0.5)
    #}
    #dev.off()

    #pdf(paste0(plot_dir,filename,datatype,"_byE13Exp.pdf"))
    #par(mar=c(6,6,6,6))
    #plot(t2_use[,2],t2_use[,6],xlab="wilcox_newValue",ylab="E13.5 log2(FPKM+1)",type="p",pch=19,cex=.7,col=color.V,main=paste0(filename,datatype),xlim=c(lowlimit-10,uplimit+10),ylim=c(0,9))
    #labels <- c(min(t2_use$peaks),round((min(t2_use$peaks)+max(t2_use$peaks))/2),max(t2_use$peaks))
    #plotrix::color.legend(360,100,400,250,legend=labels,rect.col=colors,align="rb",gradient="y",cex=0.7)
    #for(i in markers){
    #    text(t2_use[i,2],(t2_use[i,6]+0.2),i,col="black",cex=0.5)
    #}
    #dev.off()

    #barplot
    t2_use.order <- t2_use[order(t2_use$Enrichment.Score,decreasing=T),]

    bardf <- t2_use.order[1:20,c(1,2,5,6)]
    write.table(bardf,paste0(plot_dir,filename,datatype,"_top_table.txt"),row.names=F,col.names=T,sep="\t",quote=F)
    
    pdf(paste0(plot_dir,filename,datatype,"_exp_top_barplot.pdf"))
    par(mfrow = c(1:2))
    scoreVec <- bardf[,2]
    names(scoreVec) <- rownames(bardf)
    scoreVec <- rev(-scoreVec)
    par(mar=c(5.1,4.1,4.1,1.1))
    barplot(scoreVec,names.arg=rep("",20),horiz=T,las=1,col="#512D6D",xlim=c(-300,0),axes=FALSE)
    axis(1, at=seq(-0,-300,by=-50), labels=seq(0,300,by=50),cex.axis=0.6)
    legend(-200,26,legend=c("relative motif enrichment score"),fill=c("#512D6D"),bty="n",cex=0.6,xpd=T)

    expmtx <- t(as.matrix(bardf[,c(3:4)]))
    expmtx <- expmtx[2:1,20:1]
    par(mar=c(5.1,2.1,4.1,2.1))
    barplot(expmtx,beside=T,horiz=T,las=1,col=c("#F8485E","#00C1D4"),xlim=c(0,9),cex.names=0.6,axes=FALSE)
    axis(1, at=seq(0,9,by=1), labels=seq(0,9,by=1),cex.axis=0.6)
    legend(3,65,legend=c("E9.5 Expression","E13.5 Expression"),fill=c("#F8485E","#00C1D4"),bty="n",cex=0.6,xpd=T) 
    dev.off()

    bardf <- t2_use.order[(nrow(t2_use.order)-19):nrow(t2_use.order),c(1,2,5,6)]
    write.table(bardf,paste0(plot_dir,filename,datatype,"_bottom_table.txt"),row.names=F,col.names=T,sep="\t",quote=F)
    
    pdf(paste0(plot_dir,filename,datatype,"_exp_bottom_barplot.pdf"))
    par(mfrow = c(1:2))
    scoreVec <- bardf[,2]
    names(scoreVec) <- rownames(bardf)
    scoreVec <- rev(scoreVec)
    par(mar=c(5.1,4.1,4.1,1.1))
    barplot(scoreVec,names.arg=rep("",20),horiz=T,las=1,col="#512D6D",xlim=c(-300,0),axes=FALSE)
    axis(1, at=seq(-0,-300,by=-50), labels=seq(0,-300,by=-50),cex.axis=0.6)
    legend(-200,26,legend=c("relative motif enrichment score"),fill=c("#512D6D"),bty="n",cex=0.6,xpd=T)

    expmtx <- t(as.matrix(bardf[,c(3:4)]))
    expmtx <- expmtx[2:1,20:1]
    par(mar=c(5.1,2.1,4.1,2.1))
    barplot(expmtx,beside=T,horiz=T,las=1,col=c("#F8485E","#00C1D4"),xlim=c(0,9),cex.names=0.6,axes=FALSE)
    axis(1, at=seq(0,9,by=1), labels=seq(0,9,by=1),cex.axis=0.6)
    legend(3,65,legend=c("E9.5 Expression","E13.5 Expression"),fill=c("#F8485E","#00C1D4"),bty="n",cex=0.6,xpd=T)

    dev.off()

}

