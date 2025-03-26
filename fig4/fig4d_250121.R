
library(ggplot2)
library(ggbeeswarm)
library(ggrepel)
library(RColorBrewer)
library(data.table)
library(preprocessCore)
library(qusage)
sgbmMgbmLabelPwd = './mgbmSgbmLabel.csv'
sgbmMgbmLabelDf = data.frame(fread(sgbmMgbmLabelPwd))
sgbmMgbmLabelDf = na.omit(sgbmMgbmLabelDf) # filter out the patients without MRI

dim(sgbmMgbmLabelDf)
pwd4rna = './gbm_tcga_data_mrna_affymetrix_microarray.txt'

expressionMat = data.frame(fread(pwd4rna))
colnames(expressionMat) = substr(gsub('[.]', '-', colnames(expressionMat)), 1, 12)

geneNameList = expressionMat$Hugo_Symbol
expressionMat$Hugo_Symbol = NULL
expressionMat$Entrez_Gene_ = NULL
expressionMat4ciberSort = expressionMat
expressionMat4ciberSort$Gene = geneNameList
expressionMatScaled = data.frame(normalize.quantiles(as.matrix(expressionMat)))
colnames(expressionMatScaled) = colnames(expressionMat)

mgbmList = sgbmMgbmLabelDf$ID[sgbmMgbmLabelDf$bigSmall != 'SGBM']

mgbmExpressionDf = expressionMatScaled[, colnames(expressionMatScaled) %in% mgbmList]
dim(mgbmExpressionDf) # there are totally 34 samples have rna
boxplot(mgbmExpressionDf[1:10])
resInferredMigrPwd = './resMeanMigr_tol0.1_fig4F_240912.csv'
resInferredMigrDf = data.frame(fread(resInferredMigrPwd))
meanMigrArr4expr = resInferredMigrDf$migrMean[match(colnames(mgbmExpressionDf),
                                                    resInferredMigrDf$patientId)]
resGeneNameArr = c()
resCorPvalueArr = c()
resCorRhoArr = c()
for (tmpIndex in 1:nrow(mgbmExpressionDf)) {
  tmpCorTest = cor.test(unlist(mgbmExpressionDf[tmpIndex, ]),
                        meanMigrArr4expr, 
                        method = 'spearman')
  
  resGeneNameArr = c(resGeneNameArr, geneNameList[tmpIndex])
  resCorPvalueArr = c(resCorPvalueArr, tmpCorTest$p.value)
  resCorRhoArr = c(resCorRhoArr, tmpCorTest$estimate)
  
}

resExpL2corDf = data.frame(geneName = resGeneNameArr,
                           corPvalue = resCorPvalueArr,
                           corRho = resCorRhoArr)
resDegPwd = './geneList4fig5D_240912.csv'
resDegDf = data.frame(fread(resDegPwd))
resDegDfPos = resDegDf[resDegDf$corRho > 0, ] # cor: expr vs l2-volumn
resDegDfPos = resDegDfPos[order(resDegDfPos$corPvalue), ]
colnames(resDegDf)
resExpL2corDf$degPvalue = resDegDf$degPvalue[match(resExpL2corDf$geneName,
                                                   resDegDf$geneName)]
resExpL2corDf$deglog2FC = resDegDf$deglog2FC[match(resExpL2corDf$geneName,
                                                   resDegDf$geneName)]
resExpL2corDf$label = resExpL2corDf$geneName
resExpL2corDf$label[resExpL2corDf$corPvalue >= 1e-2] = NA

resExpL2corDf$col = 'grey90'
resExpL2corDf$col[resExpL2corDf$corPvalue >= 1e-2] = 'grey90'

resExpL2corDf$col[resExpL2corDf$corRho > 0 & !is.na(resExpL2corDf$label)] = 'red'
resExpL2corDf$col[resExpL2corDf$corRho < 0 & !is.na(resExpL2corDf$label)] = 'blue'

colnames(resExpL2corDf)
ggplot(resExpL2corDf, mapping = aes(x = corRho,
                                    y = -log10(corPvalue),
                                    label = label)) +
  geom_point(col = resExpL2corDf$col, size = 0.5, alpha = 0.8) +
  geom_text_repel(max.overlaps = 15, col = 'grey70', size = 3) +
  geom_hline(yintercept = -log10(0.05), col = 'grey', linetype = 'dashed') +
  theme_classic() +
  theme(text = element_text(size = 15)) +
  labs(x = 'Correlation Rho between expr and avg inferMr',
       y = '-log10(Correlation P value)')

