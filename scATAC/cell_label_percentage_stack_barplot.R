library(ArchR)
set.seed(1)

## Setting default number of Parallel threads to 8
addArchRThreads(threads = 8)
addArchRGenome("mm10")

proj <- loadArchRProject(path = "GutProject")

Cdata <- getCellColData(ArchRProj = proj, select = NULL, drop = FALSE)

tags <- c("pharynx","esophagus","lung","stomach","liver","pancreas","intestine","colon","unidentified")

numberVec <- c()

for(tag in tags){
    cellnames <- rownames(Cdata)[which(Cdata$superClusters==tag)]
    cellnumber <- c(length(which(substr(cellnames,1,7)=="E9.5_r1")),length(which(substr(cellnames,1,7)=="E9.5_r2")))
    numberVec <- c(numberVec,cellnumber)
}

numbermtx <- matrix(numberVec,nrow=2)
rownames(numbermtx) <- c("r1","r2")
colnames(numbermtx) <- tags

propmtx <- apply(numbermtx,2,function(x){
    prop1 <- 100*x[1]/(x[1]+x[2])
    prop2 <- 100*x[2]/(x[1]+x[2])
    return(c(prop1,prop2))
})

pdf("proportional_stack_barplot.pdf")
par(mar=c(8,4,4,8))
barplot(propmtx, col=c("navy","firebrick3") , border="white", xlab="",cex.names=1,cex.axis=1,las=2)
for(i in c(1,2)){
    for(j in 1:ncol(propmtx)){
        if(i==1){
            ycoord <- propmtx[i,j]/2
        }else{
            ycoord <- propmtx[i,j]/2 + propmtx[1,j]
        }
        xcoord <- 1.2*j-0.5
        outstr <- paste0(round(propmtx[i,j],1),"%")
        text(xcoord,ycoord,labels=outstr,col="white",cex=0.7)
    }
}
legend(x=11,y=50,c("replicate_1","replicate_2"),fill=c("navy","firebrick3"),cex=1,bty="n",xpd=T)
dev.off()


