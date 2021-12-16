setwd("/scratch/hz9fq/toronto_project/E185_Cdx2KO/limma_DEG")

library("jsonlite")

genetable <- read.table("Cdx2KOvsWT_DEG_table.txt",header=T)
genelist <- rownames(genetable)
anno <- read.table("/nv/vol190/zanglab/hz9fq/annotation/Mouse/mm10_gene_annotation_geneID_LenOrder_TSS50kb.bed",header=F)

toconfirm <- genelist[which(!genelist%in%anno$V5)]
symbolvec <- c()
toignore <- c()

for(i in toconfirm){
    command1 <- paste0("curl -X 'GET' 'https://mygene.info/v3/query?q=",i,"&fields=symbol&species=mouse'")
    result <- system(command1,intern = TRUE)
    class(result) <- "json"
    resObj <- fromJSON(result)

    #if there's no returned symbol, ignore that gene
    if(class(resObj$hits)!="data.frame"){
        toignore <- c(toignore,i)
        symbolvec <- c(symbolvec,0)
        next
    }
    #if there's more than 1 returned symbols, ignore that gene
    if(nrow(resObj$hits)>1){
        toignore <- c(toignore,i)
        symbolvec <- c(symbolvec,0)
        next
    }

    symbolvec <- c(symbolvec,resObj$hits[1,3])
}

maptable <- data.frame(alias=toconfirm,symbol=symbolvec)
todel <- maptable$alias[which(!maptable$symbol%in%anno$V5)]
toretain <- maptable[which(maptable$symbol%in%anno$V5),]

genelist_retain <- genelist[which(genelist%in%anno$V5)]
mapVec <- c(genelist_retain,toretain$symbol)
names(mapVec) <- c(genelist_retain,toretain$alias)

#delete symbols that are duplicate
dupsymbol <- names(table(mapVec)[which(table(mapVec)>1)])
dupalias <- names(mapVec[which(mapVec%in%dupsymbol)])
todel <- c(todel,dupalias)

#get retained gene table
genetable_retain <- genetable[which(!rownames(genetable)%in%todel),]
rownames(genetable_retain) <- mapVec[rownames(genetable_retain)]

#write out
write.table(genetable_retain,"Cdx2KOvsWT_DEG_table_confirmedHGNCsymbol.txt",col.names=T,row.names=T,sep="\t",quote=F)
write.table(rownames(genetable_retain),"Cdx2KOvsWT_DEG_genelist_confirmedHGNCsymbol.txt",col.names=F,row.names=F,sep="\t",quote=F)



