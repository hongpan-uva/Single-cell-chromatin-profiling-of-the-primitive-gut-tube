setwd("/scratch/hz9fq/toronto_project/scATAC/ArchR_workspace/GutProject_foregut")
library(GenomicRanges)

peakset <- read.table("GutProject_foregut_peakset.bed",header=F)
overlapVec <- c()

for(org in c("pharynx","esophagus","lung","stomach")){
    dir = "/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/scATAC/ArchR_workspace/GutProject/PeakCalls/"
    organpeak <- readRDS(paste0(dir,org,"-reproduciblePeaks.gr.rds"))

    gr1 <- GRanges(
        seqnames = factor(peakset[,1]),
        ranges = IRanges(start = peakset[,2],end = peakset[,3])
    )

    res <- countOverlaps(gr1, organpeak)
    overlapVec <- c(overlapVec,res)
}

overlapmtx <- matrix(overlapVec,ncol=4)

colnames(overlapmtx) <- c("pharynx","esophagus","lung","stomach")

commonVec <- apply(overlapmtx,1,function(x){
    return(all(x>0))
})

commonpeak <- peakset[commonVec,]

write.table(commonpeak,"GutProject_foregut_common_peakset.bed",row.names=F,col.names=F,quote=F,sep="\t")

