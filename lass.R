########################################
################################
################################
#######  lasso  ################
################################
################################
########################################

library("glmnet")
library("survival")

coxSigFile="DE_coexpressed_lncRNAs.uniSigExp.txt"       
rt=read.table(coxSigFile,header=T,sep="\t",row.names=1,check.names = F)             

x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))
fit=glmnet(x, y, family = "cox", maxit = 1000)
pdf("lambda.pdf")
plot(fit,xvar="lambda",label=TRUE)
dev.off()
cvfit=cv.glmnet(x, y, family="cox", maxit = 1000)
pdf("minlambda.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(geneCoef,file="geneCoef.txt",sep="\t",quote=F,row.names=F)