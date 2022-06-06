########################################
################################
################################
#######  get expression for geneList################
################################
################################
########################################
setwd("D:/projects/luad/ER_luad/")
library(limma)
dbFile="mRNA.txt"
targetFile="targetGeneList"
outputFile="targetGeneListExp.txt"

rt=read.table(dbFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

gene=read.table(targetFile, header=F, check.names=F, sep="\t")
sameGene=intersect(as.vector(gene[,1]),rownames(data))
geneExp=data[sameGene,]

out=rbind(ID=colnames(geneExp),geneExp)
write.table(out,file=outputFile,sep="\t",quote=F,col.names=F)