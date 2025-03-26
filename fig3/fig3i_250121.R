

library(data.table)

nesPwd = './nesScore231003.csv'

nesDf = data.frame(fread(nesPwd))
library(stringr)
nesDf$patientId = str_sub(nesDf$Sample, 1, 12)
nesDf$patientId = gsub('[.]', '-', nesDf$patientId)
mgbmLabelPwd = './mgbmSgbmLabelMerge230907.csv'

mgbmLabelDf = data.frame(fread(mgbmLabelPwd))
mgbmLabelDf = mgbmLabelDf[mgbmLabelDf$bigSmall != 'SGBM',]

colnames(nesDf)
mgbmLabelDf$NES = nesDf$NES[match(mgbmLabelDf$ID, nesDf$patientId)]
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
colnames(mgbmLabelDf)

table(mgbmLabelDf$bigSmall[!is.na(mgbmLabelDf$NES)])
microarrayPwd = './gbm_tcga_data_mrna_affymetrix_microarray.txt'
microarrayDf = data.frame(fread(microarrayPwd))

microarrayDf = microarrayDf[!duplicated(microarrayDf$Hugo_Symbol),]

rownames(microarrayDf) = microarrayDf$Hugo_Symbol

library(stringr)
gsub('[.]', '-', str_sub(colnames(microarrayDf), 1, 12))
mgbmMicroArrDf = microarrayDf[, gsub('[.]', '-', str_sub(colnames(microarrayDf), 1, 12)) %in% mgbmLabelDf$ID]
dim(mgbmMicroArrDf)
nesGeneList = c('S100A10','FOSL2','SPP1','CAV1',
                'ANXA1','VIM','CD44','SERPINH1',
                'LGALS3','CEBPB','ATF5','LGALS1')

length(nesGeneList)

sum(rownames(mgbmMicroArrDf) %in% nesGeneList)
colnames(mgbmMicroArrDf) = gsub('[.]', '-', str_sub(colnames(mgbmMicroArrDf), 1, 12))
resNesScoreRankSum = c()
resNesScoreSum = c()
resPatientId = c()
for (i in 1:ncol(mgbmMicroArrDf)) {
  
  
  tmpMicroArrDf = mgbmMicroArrDf[, i]
  tmpRankSum = sum(rank(tmpMicroArrDf)[rownames(mgbmMicroArrDf) %in% nesGeneList])
  tmpSum = sum(tmpMicroArrDf[rownames(mgbmMicroArrDf) %in% nesGeneList])
  tmpPatientId = colnames(mgbmMicroArrDf)[i]
  
  resNesScoreRankSum = c(resNesScoreRankSum, tmpRankSum)
  resNesScoreSum = c(resNesScoreSum, tmpSum)
  resPatientId = c(resPatientId, tmpPatientId)
  
}

resNesScore4mgbmDf = data.frame(patientId = resPatientId,
                                resNesScoreRankSum = resNesScoreRankSum,
                                resNesScoreSum = resNesScoreSum)

resNesScore4mgbmDf$label = mgbmLabelDf$bigSmall[match(resNesScore4mgbmDf$patientId,
                                                      mgbmLabelDf$ID)]
table(resNesScore4mgbmDf$label)
ggplot(resNesScore4mgbmDf, mapping = aes(x = label,
                                         y = resNesScoreSum/12,
                                         col = label)) +
  geom_boxplot(width = 0.35, outlier.shape = NA) +
  geom_beeswarm(cex = 3, size=0.5, alpha = 0.8) +
  theme_classic() +
  guides(col = 'none') +
  theme(text = element_text(size = 15)) +
  stat_compare_means(comparisons = list(c(1, 2)),
                     method = 'wilcox.test') +
  labs(x = 'MGBM group',
       y = 'NES score (Average Expression)') +
  scale_colour_manual(values=c('#e41a1c', '#377eb8'))

