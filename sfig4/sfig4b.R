
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(fgsea)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(fgsea)
library(clusterProfiler)
library(enrichplot)
library(stringr)
library(reshape2)
library(data.table)
library(ggthemes)
library(ggsci)
setFilePwd = './'
idh = read.delim(file.path(setFilePwd, 'TCGA_GBM_MRI259_IDH.txt'))
mdp = read.delim(file.path(setFilePwd, 'gbm_tcga_data_clinical_patient.txt'), comment.char = '#', na.strings = c('','[Not Available]'))
mds = read.delim(file.path(setFilePwd, 'gbm_tcga_data_clinical_sample.txt'), comment.char = '#')
fg = read.delim(file.path(setFilePwd, 'gbm_tcga_data_clinical_patient.txt'),comment.char = '#')

cl = read.csv(file.path(setFilePwd, 'mgbmSgbmLabel.csv'))
cl = cl[which(cl$modality=='MR'),]
cl$IDH = idh$IDH.status[match(cl$ID, idh$Patient)]

l2b = cl$ID[cl$bigSmall=='Big']
l2s = cl$ID[cl$bigSmall=='Small']
l1 = cl$ID[cl$bigSmall=='SGBM'&cl$IDH=='WT']
realTumorMelt = data.frame(fread('/Users/jiabao/Documents/1_github/mfgCodes/code4figs24010802/fig1/240105_volumeByLesion.csv'))
colnames(realTumorMelt) = c('patient', 'lesion', 'realVolume')
totalVolume = aggregate(realVolume ~ patient, data = realTumorMelt, sum)
orderedPatients = totalVolume$patient[order(totalVolume$realVolume)]
realTumorMelt$patient = factor(realTumorMelt$patient, levels = orderedPatients)
myMax = function(x){
  return(max(x, na.rm = T))
}
newClinicalGbmAll = fread('/Users/jiabao/Documents/1_github/mfgCodes/code4figs24010802/fig1/clinical.tsv')

newClinicalMgbm = newClinicalGbmAll

groupGDC_clinicalNew = newClinicalMgbm[,c('case_submitter_id','days_to_death','days_to_last_follow_up','vital_status')]
mgbm_t3New = apply(groupGDC_clinicalNew[,c('days_to_death','days_to_last_follow_up')],1,myMax)
mgbmDfNew = data.frame('patient' = groupGDC_clinicalNew$case_submitter_id,
                       "time" = mgbm_t3New, 
                       "status" = groupGDC_clinicalNew$vital_status)
survivalDf = unique(mgbmDfNew)
survivalDf$time = as.numeric(survivalDf$time)
sgbmMgbmLabelPwd = '/Users/jiabao/Documents/1_github/mfgCodes/therapySurvival/mgbmSgbmLabel.csv'
sgbmMgbmLabelDf = data.frame(fread(sgbmMgbmLabelPwd))
sgbmMgbmLabelDf = na.omit(sgbmMgbmLabelDf) # filter out the patients without MRI

dim(sgbmMgbmLabelDf)
gema0 = read.delim(file.path(setFilePwd, 'gbm_tcga_data_mrna_affymetrix_microarray.txt'),check.names = F)
gema = gema0[, names(gema0) %in% c('Hugo_Symbol',paste0(c(l2b,l2s),'-01'))]
gema_dif = data.frame(gene = gema$Hugo_Symbol, L2l = 0, L2s = 0, pval_w =1, pval_t =1)
for (i in 1:nrow(gema)){
  
  
  
  tmp = as.data.frame(t(gema[i,-1]))
  names(tmp)[1]='gene'
  tmp$group = ifelse(rownames(tmp) %in% paste0(l2b,'-01'),'L2l','L2s')
  gema_dif$L2l[i] = median(tmp$gene[tmp$group=='L2l'])
  gema_dif$L2s[i] = median(tmp$gene[tmp$group=='L2s'])
  gema_dif$pval_w[i] = wilcox.test(gene ~ group, data= tmp)$p.value
  gema_dif$pval_t[i] = t.test(gene ~ group, data= tmp)$p.value
}
gema_dif$log2FC = gema_dif$L2l - gema_dif$L2s
gema_dif$type = ifelse(gema_dif$pval_t<0.05 & gema_dif$log2FC>0.5,'up',
                       ifelse(gema_dif$pval_t<0.05 & gema_dif$log2FC< -0.5,'down','nosig'))
gema_dif$lab = ifelse(gema_dif$pval_t<0.05 & abs(gema_dif$log2FC)>1,gema_dif$gene,NA)
colnames(gema0) = str_sub(colnames(gema0), 1, 12)

expressionMat = gema0
geneNameList = expressionMat$Hugo_Symbol
expressionMat$Hugo_Symbol = NULL
expressionMat$Entrez_Gene_ = NULL
expressionMatScaled = data.frame(t(apply(expressionMat, 1, scale)))
colnames(expressionMatScaled) = colnames(expressionMat)

sgbmMgbmExpressionDf = gema0[, str_sub(colnames(gema0), 1, 12) %in% c('Hugo_Symbol', sgbmMgbmLabelDf$ID)]
sgbmMgbmExpressionDf[, colnames(sgbmMgbmExpressionDf)[-1]] = expressionMatScaled[, colnames(sgbmMgbmExpressionDf)[-1]]
mgbmClinicDf = data.frame(patientId = colnames(sgbmMgbmExpressionDf[, -1]))
mgbmClinicDf$mriGroup = sgbmMgbmLabelDf$bigSmall[match(mgbmClinicDf$patientId, sgbmMgbmLabelDf$ID)]
mgbmClinicDf$mriGroup = factor(mgbmClinicDf$mriGroup,
                               levels = c('Small', 'SGBM', 'Big'))
mgbmClinicDf$time = survivalDf$time[match(mgbmClinicDf$patientId, survivalDf$patient)]
mgbmClinicDf$status = survivalDf$status[match(mgbmClinicDf$patientId, survivalDf$patient)]
realTumorMeltTotalVolumeDf = data.frame(aggregate(realVolume ~ patient, data = realTumorMelt, FUN = sum))
realTumorMeltL2volumeDf = realTumorMelt[realTumorMelt$lesion =='L2', ]

realTumorNumDf = data.frame(aggregate(realVolume ~ patient, data = realTumorMelt, FUN = length))
colnames(realTumorNumDf) = c('patient', 'lesionNum')
mgbmClinicDf$L2volume = realTumorMeltL2volumeDf$realVolume[match(mgbmClinicDf$patientId, realTumorMeltL2volumeDf$patient)]
mgbmClinicDf$totalVolume = realTumorMeltTotalVolumeDf$realVolume[match(mgbmClinicDf$patientId, realTumorMeltTotalVolumeDf$patient)]
mgbmClinicDf$lesionNum = realTumorNumDf$lesionNum[match(mgbmClinicDf$patientId, realTumorNumDf$patient)]

mgbmClinicDf$relativeL2 = mgbmClinicDf$L2volume/mgbmClinicDf$totalVolume
mgbmExpressionDf = sgbmMgbmExpressionDf[, c('Hugo_Symbol', mgbmClinicDf$patientId)]
colnames(mgbmExpressionDf)[-1] == mgbmClinicDf$patientId
resGeneNameArr = c()
resCorPvalueArr = c()
resCorRhoArr = c()
for (tmpIndex in 1:nrow(mgbmExpressionDf)) {
  tmpCorTest = cor.test(unlist(mgbmExpressionDf[tmpIndex, ][-1]),
                        mgbmClinicDf$L2volume, 
                        method = 'spearman')
  
  resGeneNameArr = c(resGeneNameArr, mgbmExpressionDf$Hugo_Symbol[tmpIndex])
  resCorPvalueArr = c(resCorPvalueArr, tmpCorTest$p.value)
  resCorRhoArr = c(resCorRhoArr, tmpCorTest$estimate)
  
  
  
}
setTmpGeneName = 'UPP1'
setTmpGeneName = 'IGFBP6'
setTmpGeneName = 'IGFBP7'
setTmpGeneName = 'ANXA2'

tmpCorDf = data.frame(expr = unlist(mgbmExpressionDf[mgbmExpressionDf$Hugo_Symbol == setTmpGeneName, ][-1]),
                      volume = mgbmClinicDf$L2volume)

tmpCorDf = na.omit(tmpCorDf)
cor.test(tmpCorDf$expr, tmpCorDf$volume, method = 'spearman')
tmpDf4geneExpr = data.frame(t(expressionMatScaled[geneNameList == setTmpGeneName, ]))
colnames(tmpDf4geneExpr) = 'Expression'
resExpL2corDf = data.frame(geneName = resGeneNameArr,
                           corPvalue = resCorPvalueArr,
                           corRho = resCorRhoArr)

resExpL2corDf = resExpL2corDf[order(resExpL2corDf$corRho, decreasing = T), ]
resExpL2corDf$label = resExpL2corDf$geneName
resExpL2corDf$label[resExpL2corDf$corPvalue > 0.01] = NA
resExpL2corDf$col = 'black'
resExpL2corDf$col[resExpL2corDf$corPvalue > 0.01] = 'grey90'

library(ggrepel)

gema_dif

table(gema_dif$type)
resExpL2corDf$degPvalue = gema_dif$pval_t[match(resExpL2corDf$geneName, gema_dif$gene)]
resExpL2corDf$deglog2FC = gema_dif$log2FC[match(resExpL2corDf$geneName, gema_dif$gene)]
tmpDf = resExpL2corDf[resExpL2corDf$corPvalue < 0.05, ]
tmpDf = tmpDf[tmpDf$deglog2FC > 0.5, ]
tmpDf = tmpDf[tmpDf$degPvalue < 0.05, ]
dim(tmpDf)

cat(tmpDf$geneName)
cat(tmpDf$geneName, sep = '", "')
tmpCorDf = data.frame(L2 = mgbmClinicDf$L2volume,
                      express = unlist(mgbmExpressionDf[mgbmExpressionDf$Hugo_Symbol == 'LOX', ][-1]))

tmpCorDf = na.omit(tmpCorDf)

cor.test(tmpCorDf$express,
         tmpCorDf$L2, method = 'spearman')

cor.test(unlist(mgbmExpressionDf[mgbmExpressionDf$Hugo_Symbol == 'LOX', ][-1]),
         mgbmClinicDf$L2volume, method = 'spearman')
resExpL2corDf$col[resExpL2corDf$geneName %in% tmpDf$geneName] = 'red'
resExpL2corDf$col[resExpL2corDf$geneName %in% gema_dif$gene[gema_dif$type == 'down']] = 'blue'

resExpL2corDf$col[resExpL2corDf$corPvalue >= 0.05] = 'grey90'
ggplot() +
  geom_point(resExpL2corDf[resExpL2corDf$col == 'grey90', ], mapping = aes(x = corRho,
                                                                           y = -log10(corPvalue)),
             col = 'grey90', alpha = 0.5, size = 0.5) +
  
  
  geom_point(resExpL2corDf[resExpL2corDf$col == 'black', ], 
             mapping = aes(x = corRho,
                           y = -log10(corPvalue),
             ),
             col = 'grey', alpha = 0.8, size = 0.5) +
  geom_point(resExpL2corDf[resExpL2corDf$col == 'red', ], 
             mapping = aes(x = corRho,
                           y = -log10(corPvalue),
             ),
             col = resExpL2corDf$col[resExpL2corDf$col == 'red'], alpha = 0.8, size = 0.5) +
  
  geom_point(resExpL2corDf[resExpL2corDf$col == 'blue', ], 
             mapping = aes(x = corRho,
                           y = -log10(corPvalue),
             ),
             col = resExpL2corDf$col[resExpL2corDf$col == 'blue'], alpha = 0.8, size = 0.5) +
  
  geom_hline(yintercept = -log10(0.05), col = 'grey', linetype = 'dashed') +
  geom_text_repel(resExpL2corDf, mapping = aes(x = corRho,
                                               y = -log10(corPvalue),
                                               label = label),
                  max.overlaps = 23, alpha = 0.5, size = 2) +
  theme_classic() +
  theme(text = element_text(size = 15)) +
  labs(x = 'Correlation Rho between expression and L2 volume',
       y = '-log10(Correlation P value)')
resExpL2corDf$geneName[resExpL2corDf$col == 'red']
length(resExpL2corDf$geneName[resExpL2corDf$col == 'red'])

