########################################
################################
################################
#######  uniCox  ################
################################
################################
########################################
library(survival)
pFilter=0.05                                                    
rt=read.table("DE_coexpressed_lncRNAsTime.txt",header=T,sep="\t",check.names=F,row.names=1)     
rt$futime=rt$futime/365                                                   
#rt[,3:ncol(rt)]=log2(rt[,3:ncol(rt)]+1)

outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
  if(sd(rt[,i])<0.01){
    next}
  a=rt[,i]<=median(rt[,i])
  diff=survdiff(Surv(futime, fustat) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  fit=survfit(Surv(futime, fustat) ~ a, data = rt)
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if((pValue<pFilter) & (coxP<pFilter)){
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxP))}
}
write.table(outTab,file="DE_coexpressed_lncRNAs.uniCox.txt",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="DE_coexpressed_lncRNAs.uniSigExp.txt",sep="\t",row.names=F,quote=F)