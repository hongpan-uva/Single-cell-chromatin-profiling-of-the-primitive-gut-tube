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

    t2_use$value <- -1
    count <- 0
    for(i in 1:nrow(t2_use)){
        if(rownames(t2_use)[i]%in%rownames(FPKM_E135)){
            t2_use$value[i]<- log2((FPKM_E135[rownames(t2_use)[i],set]+0.01)/(FPKM_E95[rownames(t2_use)[i],set]+0.01))
            count <- count+1
        }
    }
    print(paste0("Number of TF with expression: ",count))
    print(summary(t2_use$value))

    markers=unique(c(rownames(t2_use)[order(t2_use[,2],decreasing=T)][1:20],rownames(t2_use)[order(t2_use[,2],decreasing=F)][1:20],rownames(t2_use)[order(t2_use[,5],decreasing=T)][1:10],rownames(t2_use)[order(t2_use[,5],decreasing=F)][1:10]))

    #t2_use$rank <- rank(t2_use$wilcox.new,ties.method="first")
    
    pdf(paste0(plot_dir,filename,datatype,"_byExpLogfc.pdf"))
    par(mar=c(6,6,6,6))
    plot(t2_use[,2],t2_use[,5],xlab="wilcox_newValue",ylab="Expression log(FC)",type="p",pch=19,cex=.7,col=color.V,main=paste0(filename,datatype),xlim=c(lowlimit-10,uplimit+10))
    labels <- c(min(t2_use$peaks),round((min(t2_use$peaks)+max(t2_use$peaks))/2),max(t2_use$peaks))
    plotrix::color.legend(360,100,400,250,legend=labels,rect.col=colors,align="rb",gradient="y",cex=0.7)
    for(i in markers){
        text(t2_use[i,2],(t2_use[i,5]+0.2),i,col="black",cex=0.5)
    }
    dev.off()

}

