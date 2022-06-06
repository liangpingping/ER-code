########################################
################################
################################
#######  relationship between risk and gender stages  ################
################################
################################
########################################
setwd("D:/projects/luad/ER_luad/")   

risk=read.table("trainRisk.txt",header=T,sep="\t",check.names=F)   
cli=read.table("clinical_stage.txt",sep="\t",check.names=F,header=T,quote="")    
merged <- merge(risk,cli,by.x = "id", by.y = "id")

library(ggpubr)
library(ggplot2)

boxplot <- ggboxplot(merged, x = "gender", y = "riskScore",
                     ylab = "riskScore",
                     fill = "gender",palette=c("#4DBBD5FF","#E64B35FF"),notch=TRUE) +
  stat_compare_means(comparisons = list(c("male", "female"))) 
pdf(file="riskscore_gender_box.pdf", width=5, height=4.5)
print(boxplot)
dev.off()