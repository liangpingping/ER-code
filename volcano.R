########################################
################################
################################
####### volcano ################
################################
################################
########################################
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(airway)

file1 <- as.data.frame(read.table('all_DElncRNAs.xls',header=T))
res1 <- file1
vol1 <- EnhancedVolcano(res1,lab = res1$gene,
                        x = 'logFC',y = 'fdr',
                        xlim = c(-12,12),
                        pCutoff = 0.05,FCcutoff = 1,
                        pointSize = 2.0, labSize = 4.0,
                        col=c('black', 'blue', 'green', 'red'),colAlpha = 1,
                        legendLabels=c('NS','Log2FC','P-value','P-value & Log2FC'),
                        legendPosition = 'top')

pdf("volcano_DE_coexpressedlncRNAs.pdf")
plot(vol1)
dev.off()