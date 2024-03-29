########################################
################################
################################
#######  get DElncRNAs ################
################################
################################
########################################
setwd("D:/projects/luad/ER_luad")
library(limma)

inputFile="colncrnaExp.txt"   
fdrFilter=0.05                
logFCfilter=1             

rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2", "1", group)
conNum=length(group[group==1])       
treatNum=length(group[group==0])     
grade=c(rep(1,conNum), rep(2,treatNum))

outTab=data.frame()
for(i in row.names(data)){
  geneName=unlist(strsplit(i,"\\|",))[1]
  geneName=gsub("\\/", "_", geneName)
  rt=rbind(expression=data[i,],grade=grade)
  rt=as.matrix(t(rt))
  wilcoxTest=wilcox.test(expression ~ grade, data=rt)
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  pvalue=wilcoxTest$p.value
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
  }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)

write.table(outTab,file="all_DElncRNAs.xls",sep="\t",row.names=F,quote=F)

outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="DElncRNAs.xls",sep="\t",row.names=F,quote=F)
write.table(outDiff,file="DElncRNAs.txt",sep="\t",row.names=F,quote=F)

exp=rbind(ID=colnames(data[as.vector(outDiff[,1]),]),data[as.vector(outDiff[,1]),])
write.table(exp,file="DE_coexpressed_lncRNAsExp.txt",sep="\t",col.names=F,quote=F)