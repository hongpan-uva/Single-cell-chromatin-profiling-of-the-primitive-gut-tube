motifoverlap="/scratch/hz9fq/toronto_project/intestine_CrossStages/motif_overlap/peak_motifoverlap3.txt"
cor_dir="/scratch/hz9fq/toronto_project/intestine_CrossStages/signal_forMotifAnalysis3/cor/"
plot_dir="/scratch/hz9fq/toronto_project/intestine_CrossStages/signal_forMotifAnalysis3/cor/wilcox_plots/"
peakset="/scratch/hz9fq/toronto_project/intestine_CrossStages/peaks/mergedPeak_400bp.bed"
cor_files=c("E13intestine_E16colon_cor.txt","E13intestine_E16smallIntestine_cor.txt")

motifoverlap="/scratch/hz9fq/toronto_project/stomach_CrossStages/motif_overlap/peak_motifoverlap3.txt"
cor_dir="/scratch/hz9fq/toronto_project/stomach_CrossStages/signal_forMotifAnalysis3/cor/"
plot_dir="/scratch/hz9fq/toronto_project/stomach_CrossStages/signal_forMotifAnalysis3/cor/wilcox_plots/"
peakset="/scratch/hz9fq/toronto_project/stomach_CrossStages/peaks/mergedPeak_400bp.bed"
cor_files=c("E13stomach_E16forestomach_cor.txt","E13stomach_E16hindstomach_cor.txt")

#prepare FPKM tables
FPKM_col = read.table("/scratch/hz9fq/toronto_project/ALL_RNA_collection/mod_exp/FPKM_collection.txt",header=T)
FPKM_use <- FPKM_col[,7:22]
FPKM_use <- cbind(apply(FPKM_use[,1:2],1,mean),
                   apply(FPKM_use[,5:6],1,mean),
                   apply(FPKM_use[,13:14],1,mean),
                   apply(FPKM_use[,15:16],1,mean),
                   apply(FPKM_use[,9:10],1,mean),
                   apply(FPKM_use[,11:12],1,mean))
colnames(FPKM_use) <- c("E13intestine","E13stomach","E16smallIntestine","E16colon","E16forestomach","E16hindstomach")
rownames(FPKM_use) <- toupper(FPKM_col$Gene.ID)

motif.tab <- read.table(motifoverlap,header=T,stringsAsFactors=F)
peakNumber <- apply(motif.tab,2,sum)

for(set in 1:2){
    file <- cor_files[set]

    print(file)
    t2 <- read.table(paste0(cor_dir,file),header=T,stringsAsFactors=F)
    rownames(t2) <- t2[,1]

    #wilcox.new
    datatype = "_wilcox.new"
    t2_use <- t2[,c(1,8)]

    uplimit = max(abs(t2_use[,2]))
    lowlimit = -uplimit

    overlap1 <- intersect(rownames(t2_use),rownames(FPKM_use))
    t2_use <- t2_use[overlap1,]

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
    background <- strsplit(filename,"_")[[1]][1]
    foreground <- strsplit(filename,"_")[[1]][2]

    horizon_pal=c("#000075","#2E00FF","#9408F7","#C729D6","#FA4AB5","#FF6A95","#FF8B74","#FFAC53","#FFCD32","#FFFF60")
    colors <- colorRampPalette(horizon_pal)(100)
    color.V <- colors[as.numeric(t2_use$V4)]

    t2_use$value <- -1
    count <- 0
    for(i in 1:nrow(t2_use)){
        if(rownames(t2_use)[i]%in%rownames(FPKM_use)){
            t2_use$value[i] <- log2((FPKM_use[rownames(t2_use)[i],foreground]+0.01)/(FPKM_use[rownames(t2_use)[i],background]+0.01))
            count <- count+1
        }
    }
    print(paste0("Number of TF with expression: ",count))
    print(summary(t2_use$value))

    markers=unique(c(rownames(t2_use)[order(t2_use[,2],decreasing=T)][1:20],rownames(t2_use)[order(t2_use[,2],decreasing=F)][1:20],rownames(t2_use)[order(t2_use[,5],decreasing=T)][1:10],rownames(t2_use)[order(t2_use[,5],decreasing=F)][1:10]))

    pdf(paste0(plot_dir,filename,datatype,"_byExpLogfc.pdf"))
    par(mar=c(6,6,6,6))
    plot(t2_use[,2],t2_use[,5],xlab="wilcox_newValue",ylab="Expression log(FC)",type="p",pch=19,cex=.7,col=color.V,main=paste0(filename,datatype),xlim=c(lowlimit-10,uplimit+10))
    labels <- c(min(t2_use$peaks),round((min(t2_use$peaks)+max(t2_use$peaks))/2),max(t2_use$peaks))
    coords <- par("usr")
    width <- coords[2]-coords[1]
    height <- coords[4]-coords[3]
    plotrix::color.legend(coords[2]+width/20,coords[4]-height/4,coords[2]+width/10,coords[4],legend=labels,rect.col=colors,align="rb",gradient="y",cex=0.7)
    for(i in markers){
        text(t2_use[i,2],(t2_use[i,5]+0.2),i,col="black",cex=0.5)
    }
    dev.off()

}
