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

    colnames(t2_use)[2] <- "Enrichment.Score"

    t2_use$E135exp <- FPKM_use[rownames(t2_use),background]
    t2_use$E165exp <- FPKM_use[rownames(t2_use),foreground]

    t2_use.order <- t2_use[order(t2_use$Enrichment.Score,decreasing=T),]

    bardf <- t2_use.order[1:20,c(1,2,5,6)]
    write.table(bardf,paste0(plot_dir,filename,datatype,"_top_table.txt"),row.names=F,col.names=T,sep="\t",quote=F)

    bardf <- t2_use.order[(nrow(t2_use.order)-19):nrow(t2_use.order),c(1,2,5,6)]
    write.table(bardf,paste0(plot_dir,filename,datatype,"_bottom_table.txt"),row.names=F,col.names=T,sep="\t",quote=F)
    
}
