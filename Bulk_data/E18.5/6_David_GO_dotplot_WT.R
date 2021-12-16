setwd("//SOM22F1L5D3/projects/toronto_project/revision/result/E185_Cdx2KO/GOterm")

#t1 <- read.table("Cdx2KOvsWT_gene_associatedwith_openRegion.genelist.GO_BP.txt",header=T,sep="\t")
t1 <- read.table("WTvsCdx2KO_gene_associatedwith_closeRegion.genelist.GO_BP.txt",header=T,sep="\t")

#Benjamini: Benjamini in DAVID requests adjusted p-values by using the linear step-up method of Benjamini and Hochberg (1995). 
term <- c()
term <- c(term,t1$Term[which(t1$Benjamini<0.05)])
term <- unique(term)

df1_vec <- c()
for(t in term){
	if(t%in%t1$Term){
		df1_vec <- c(df1_vec,t1$Count[which(t1$Term==t)],t1$Benjamini[which(t1$Term==t)])
	}
	else{
		df1_vec <- c(df1_vec,0,1)
	}
}

df1 <- as.data.frame(matrix(df1_vec,ncol=2,byrow=T))
colnames(df1) <- c("Count","Benjamini")
df1$Term <- term
df1$GeneType <- "WT-sepecific genes associated with Lost peak"
df1$Term <- factor(df1$Term, levels = rev(df1$Term))

library(ggplot2)
#solarExtra_pal<-c("#3361A5","#248AF3","#14B3FF","#88CEEF","#C1D5DC","#EAD397","#FDB31A","#E42A2A","#A31D1D")
customzied_pal<-c("#3361A5","#14B3FF","#C1D5DC","#D0CDAC","#DFC47B","#EEBC4B","#FDB31A","#E78E1B","#D0681C","#CB350F","#C50202")

p1 <- ggplot(df1) + aes(x=GeneType,y=Term,colour=Benjamini,size=Count) + scale_colour_gradientn(colours=rev(customzied_pal),oob=scales::squish,trans="log10",limits=c(1e-5,1),name="Benjamini") + scale_size(limits=c(0,202),range = c(0, 15),breaks=c(200,150,100,50)) + geom_point() + guides(color = guide_colourbar(order=1), size = guide_legend(order=2)) + theme(axis.text.y = element_text(size = 12))
ggsave("WT_genes_withLostPeaks_GOterm_BP.pdf",height=12,width=12,p1)
