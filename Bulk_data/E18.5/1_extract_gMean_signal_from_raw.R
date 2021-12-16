setwd("/scratch/hz9fq/toronto_project/E185_Cdx2KO/raw_data")

filenames <- c("251486815054_2.txt","251486815053_3.txt","251486815052_3.txt","251486815052_1.txt","251486815053_1.txt","251486815054_1.txt","251486815054_3.txt","251486815052_4.txt","251486815053_4.txt")

for(f in filenames){
    rawdf<-read.delim(f,sep="\t",header=T,quote="",fill=F)
    if(f=="251486815054_2.txt"){
        gMSdf <- cbind(rawdf$Feature.Extraction.Software.ProbeName,rawdf$Feature.Extraction.Software.gMeanSignal)
    }else{
        gMSdf <- cbind(gMSdf,rawdf$Feature.Extraction.Software.gMeanSignal)
    }
}

colnames(gMSdf) <- c("Scan.REF","ileum_Cdx2KO_rep1","ileum_Cdx2KO_rep2","ileum_Cdx2KO_rep3","ileum_WT_rep1","ileum_WT_rep2","ileum_WT_rep3","esophagus_WT_rep1","esophagus_WT_rep2","esophagus_WT_rep3")
gMSdf <- as.data.frame(gMSdf)
for(i in 2:10){
    gMSdf[,i] <- as.numeric(gMSdf[,i])
}
write.table(gMSdf,"gMeanSignal_matrix.txt",col.names=T,row.names=F,quote=F,sep="\t")

#filter probe
gMSdf <- gMSdf[substr(gMSdf$Scan.REF,1,2)=="A_",]

#read in annotation
genetable <- read.table("../A-MEXP-724.adf.edited.txt",header=T,sep="\t")

genenames <- sapply(gMSdf$Scan.REF,function(x){
    genetable$Composite.Element.Name[which(genetable$Reporter.Name==x)][1]
})

#check if there's biomarker mapped to more than one gene symbols
single_map <- c()
multi_map <- c()
for(i in 1:length(genenames)){
    if(length(genenames[[i]]) == 1){
        single_map <- c(single_map,i)
    }else{
        multi_map <- c(multi_map,i)
    }
}

genenames <- unlist(genenames)
#get duplicate probes
dupprobe <- names(table(names(genenames)))[which(table(names(genenames))>1)]

#for genes profiled by duplicate probes (same probes), take the median as gMSdfression
toadd <- gMSdf[1,]
toadd <- toadd[-1,]
todel_1 = c()

for(i in dupprobe){
    tmp <- gMSdf[which(gMSdf$Scan.REF==i),2:10]
    todel_1 <- c(todel_1,rownames(tmp))
    medianExp <- apply(tmp, 2, median)
    medianExp <- as.integer(medianExp)
    medianExp <- c(i, medianExp)
    toadd <- rbind(toadd,as.data.frame(matrix(medianExp,nrow=1)))
}

colnames(toadd) <- colnames(gMSdf)
toaddgene <- genenames[toadd[,1]]

#do deletion for the first time
gMSdf <- gMSdf[-as.numeric(todel_1),]
rownames(gMSdf) <- 1:nrow(gMSdf)
genenames <- genenames[-as.numeric(todel_1)]

#add median gMSdfression back
gMSdf <- rbind(gMSdf, toadd)
genenames <- c(genenames, toaddgene)

#for genes with multi transcriptomes, retain the first one and record the related probe id
multitrans <- names(table(genenames))[which(table(genenames)>1)]
dictionary <- data.frame(Symbol=0,Probe=0)
todel_2 <- c()
todel_probe <- c()
count <- 0
test_probe <- c()
test_gene <- c()

for(i in unique(genenames)){
    if(!i %in% multitrans){
        dictionary <- rbind(dictionary,c(i,names(genenames[which(genenames==i)])))
        count <- count+1
    }else{
        pbs <- names(genenames[which(genenames==i)])
        count <- count+length(pbs)
        dictionary <- rbind(dictionary,c(i,pbs[1]))
        todel_2 <- c(todel_2,rownames(gMSdf)[which(gMSdf$Scan.REF%in%pbs[2:length(pbs)])])
        todel_probe <- c(todel_probe, pbs[2:length(pbs)])
    }
}

dictionary <- dictionary[-1,]
rownames(dictionary) <- dictionary$Probe

#do deletion for the second time
gMSdf <- gMSdf[-as.numeric(todel_2),]
rownames(gMSdf) <- 1:nrow(gMSdf)
genenames <- genenames[-as.numeric(todel_2)]

rownames(gMSdf) <- dictionary[gMSdf$Scan.REF,1]
gMSdf <- gMSdf[,-1]

write.table(gMSdf,"Symbol_gMeanSignal_matrix.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(dictionary,"Symbol_probe_dictionary.txt",row.names=T,col.names=T,sep="\t",quote=F)

gMSdf.use <- gMSdf[,1:6]

pdf("gMS_check_markers.pdf")
checkgenes <- c("Cdx2","Cdx1","Tff3","Sox2","Sox6","Osr1")
colors <- c(rep("#66c2a5",3),rep("#fc8d62",3),rep("#8da0cb",3))
uplimit <- c(15000,2000,21000,2000,800,1000)
for(i in 1:length(checkgenes)){
    plotdf <- gMSdf.use[checkgenes[i],]
    plotvec <- as.numeric(plotdf[1,])
    names(plotvec) <- colnames(gMSdf.use)
    barplot(plotvec,col=colors,main=checkgenes[i],ylim=c(0,uplimit[i]),cex.names=0.5)
}
dev.off()

library(preprocessCore)

gMSmtx <- as.matrix(gMSdf.use)
mode(gMSmtx)<-"numeric"
gMSmtx.norm <- normalize.quantiles(gMSmtx)

rownames(gMSmtx.norm) <- rownames(gMSdf.use)
colnames(gMSmtx.norm) <- colnames(gMSdf.use)

pdf("gMSnorm_check_markers.pdf")
checkgenes <- c("Cdx2","Cdx1","Tff3","Sox2","Sox6","Osr1")
colors <- c(rep("#66c2a5",3),rep("#fc8d62",3),rep("#8da0cb",3))
uplimit <- c(15000,2000,21000,2000,800,1000)
for(i in 1:length(checkgenes)){
    plotdf <- gMSmtx.norm[checkgenes[i],]
    names(plotdf) <- colnames(gMSdf.use)
    barplot(plotdf,col=colors,main=checkgenes[i],ylim=c(0,uplimit[i]),cex.names=0.5)
}
dev.off()

summary(gMSmtx.norm)

gMSdf.use.norm <- as.data.frame(gMSmtx.norm)
write.table(gMSdf.use.norm,"Symbol_gMeanSignal_norm_matrix.txt",row.names=T,col.names=T,sep="\t",quote=F)

