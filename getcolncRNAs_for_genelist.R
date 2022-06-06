########################################
################################
################################
#######  get co-expressed lncRNAs for geneList################
################################
################################
########################################
setwd("D:/projects/luad/ER_luad/")
library(limma)
file1="targetGeneListExp.txt"
file2="lncRNA.txt"
corFilter=0.4                                                        
pvalueFilter=0.001 

rt = read.table(file1,header=T,sep="\t",check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mrna=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mrna=avereps(mrna)
mrna=mrna[rowMeans(mrna)>0.5,]

rt = read.table(file2,header=T,sep="\t",check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
lncRNA=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
lncRNA=avereps(lncRNA)
lncRNA=lncRNA[rowMeans(lncRNA)>0.5,]

outTab=data.frame()

for(i in row.names(lncRNA)){
  if(sd(lncRNA[i,])>0.5){
    for(j in row.names(mrna)){
      x=as.numeric(lncrna[i,])
      y=as.numeric(mrna[j,])
      corT=cor.test(x,y)
      cor=corT$estimate
      pvalue=corT$p.value
      if((cor>corFilter) & (pvalue<pvalueFilter)){
        outTab=rbind(outTab,cbind(mrna=j,lncrna=i,cor,pvalue,Regulation="postive"))
      }
      if((cor< -corFilter) & (pvalue<pvalueFilter)){
        outTab=rbind(outTab,cbind(mrna=j,lncrna=i,cor,pvalue,Regulation="negative"))}
    }}}
write.table(file="cor_correlationship_filtered.txt",outTab,sep="\t",quote=F,row.names=F)
colncrna=unique(as.vector(outTab[,"lncrna"]))
lncrnaExp=lncRNA[colncrna,]
lncrnaExp=rbind(ID=colnames(lncrnaExp),lncrnaExp)
write.table(lncrnaExp,file="colncrnaExp.txt",sep="\t",quote=F,col.names=F)