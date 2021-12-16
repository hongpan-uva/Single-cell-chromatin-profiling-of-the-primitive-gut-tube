#reproduce from here
set.seed(1)

Sox2OE <- read.table("/scratch/hz9fq/toronto_project/E185_Cdx2KO/limma_DEG/Cdx2KOvsWT_DEG_genelist.txt",header=F)
Sox2OE <- Sox2OE$V1

#signature genes boxpot
library("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

library("rjson")
setwd("/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/TCGA/FPKM")
FPKM_dat_use <- read.table("PAAD_FPKM.txt",header=T,row.names=1,sep='\t',stringsAsFactors=F)
Zscore_dat_use <- read.table("PAAD_Zscore.txt",header=T,row.names=1,sep='\t',stringsAsFactors=F)

k="SOX2"
#cutoff
signal <- as.numeric(FPKM_dat_use[k,])
names(signal) <- colnames(FPKM_dat_use)
g1 <- which(signal<1)
g2 <- which(signal>=2)

print(c(length(g1),length(g2)))

glist <- colnames(Zscore_dat_use)
names(glist) <- colnames(Zscore_dat_use)
glist[g1] <- paste0(k,"_low")
glist[g2] <- paste0(k,"_high")    
glist <- factor(glist,levels=c(paste0(k,"_low"),paste0(k,"_high")))

glist <- glist[!is.na(glist)]

#GSEA running enrichment
library("clusterProfiler")

#read in count matrix
setwd("/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/TCGA/count/")
Count_dat <- read.table("PAAD_Count.txt",header=T,stringsAsFactors=F)
#rearrange count matrix by the order of case id in the Zscore matrix
Count_dat <- Count_dat[,colnames(Zscore_dat_use)]

#limma get differential genes
library(DESeq2)
library(edgeR)
set.seed(1)

k="SOX2"

coldata <- data.frame(compare1=factor(glist, levels=c(paste0(k,"_low"),paste0(k,"_high"))))
cnt_dds <- Count_dat[,rownames(coldata)]
dds <- DESeqDataSetFromMatrix(countData = cnt_dds,
                      colData = coldata,
                      design= ~ compare1)

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients

res <- lfcShrink(dds, coef=paste0("compare1_",k,"_high","_vs_",k,"_low"), type="apeglm")
nrow(res[which(res$padj<=0.05 & res$log2FoldChange>=log2(1.5)),])
#filter dif peak
res.f <- res[which(res$padj<=0.05 & res$log2FoldChange>=log2(1.5)),]
res.f <- res.f[order(res.f$log2FoldChange,decreasing=T),]

siglist <- rownames(res.f)


#get the ranked gene list
#16.5
dir="/scratch/hz9fq/toronto_project/ALL_RNA_collection/mod_exp/"
t2 <- read.table(paste0(dir,"FPKM_collection.txt"),header=T,stringsAsFactors=F)
#stomach expression
left_exp <- apply(t2[,15:18],1,mean)
#intestine expression
right_exp <- apply(t2[,19:22],1,mean)
exp_dat <- data.frame(left_exp,right_exp)
rownames(exp_dat) <- t2$Gene.ID

#13.5
dir="/scratch/hz9fq/toronto_project/ALL_RNA_collection/mod_exp/"
t2 <- read.table(paste0(dir,"FPKM_collection.txt"),header=T,stringsAsFactors=F)
#stomach expression
left_exp <- apply(t2[,11:12],1,mean)
#pancreas expression
right_exp <- apply(t2[,9:10],1,mean)
#intestine expression
#right_exp <- apply(t2[,7:8],1,mean)
exp_dat <- data.frame(left_exp,right_exp)
rownames(exp_dat) <- t2$Gene.ID

#9.5
dir="/scratch/hz9fq/toronto_project/E95_scRNA/"
t2 <- read.table(paste0(dir,"E95_Expression_matrix_byOrgan.txt"),header=T,stringsAsFactors=F)
#stomach expression
left_exp <- t2$stomach
#lung expression
#left_exp <- t2$lung
#pancreas expression
right_exp <- t2$pancreas
#intestine expression
#right_exp <- t2$intestine
exp_dat <- data.frame(left_exp,right_exp)
rownames(exp_dat) <- rownames(t2)

exp_dat$fc <- ((exp_dat$left_exp+0.1)/(exp_dat$right_exp+0.1))
exp_dat <- exp_dat[order(exp_dat$fc,decreasing=T),]
geneList <- exp_dat$fc
names(geneList) <- rownames(exp_dat)
#use log fold change so that there are negative values	
geneList<-log2(geneList)


#caculate GSEA from here
#transfer human gene symbol to mouse gene symbol
genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = siglist , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
GeneSet <- unique(genesV2$MGI.symbol)
print(length(GeneSet))

#refine the genelist by Sox2OE
GeneSet2 <- Sox2OE[which(Sox2OE%in%GeneSet)]
length(GeneSet2)
write.table(GeneSet2,"cutoff_1_2_refined/PAAD_SOX2high_signature_mm10.txt",row.names=F,col.names=F,quote=F,sep="\t")

#record human homologues
genesV3 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = GeneSet2 , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
GeneSet2_hg38 <- unique(genesV3$HGNC.symbol)
length(GeneSet2_hg38)
write.table(GeneSet2_hg38,"cutoff_1_2_refined/PAAD_SOX2high_signature_hg38.txt",row.names=F,col.names=F,quote=F,sep="\t")

GeneSet3<-GeneSet2[which(GeneSet2%in%names(geneList))]
TTG <- data.frame(term=rep("term1",length(GeneSet3)),gene=GeneSet3)

#cluster profiler universe GSEA method
uni <- GSEA(
    geneList,
    exponent = 1,
    minGSSize = 10,
    maxGSSize = 5000,
    eps = 1e-100,
    pvalueCutoff = 10,
    pAdjustMethod = "BH",
    TERM2GENE = TTG,
)

print(uni@result)

pdf(paste0("../result_plots/cutoff_1_2_refined/PAAD_GSEA_",k,"_RNAE135.pdf"))
p1 <- gseaplot(uni, geneSetID = 1, by = "runningScore", title=paste0(k," high signature genes"))
print(p1)
dev.off()

