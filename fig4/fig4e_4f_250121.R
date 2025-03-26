

library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(fgsea)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(clusterProfiler)
library(enrichplot)
library(stringr)
library(reshape2)
library(data.table)
library(ggthemes)
library(ggsci)

library(dplyr)
library(DT)
library(survminer)
library(survival)
library(preprocessCore)
library(matrixStats)
library(circlize)
idhStatusDf = read.delim('./TCGA_GBM_MRI259_IDH.txt')

mgbmSgbmLabelDf = read.csv('./mgbmSgbmLabel.csv')

mgbmSgbmLabelDf = mgbmSgbmLabelDf[which(mgbmSgbmLabelDf$modality=='MR'),]

mgbmSgbmLabelDf$IDH = idhStatusDf$IDH.status[match(mgbmSgbmLabelDf$ID, idhStatusDf$Patient)]

l2largeSampleArr = mgbmSgbmLabelDf$ID[mgbmSgbmLabelDf$bigSmall=='Big']
l2smallSampleArr = mgbmSgbmLabelDf$ID[mgbmSgbmLabelDf$bigSmall=='Small']
sgbmSampleArr = mgbmSgbmLabelDf$ID[mgbmSgbmLabelDf$bigSmall=='SGBM' & mgbmSgbmLabelDf$IDH=='WT']
allIdhwtSampleIdArr = c(l2largeSampleArr, l2smallSampleArr)
length(sgbmSampleArr)
length(allIdhwtSampleIdArr)

myMax = function(x){
  return(max(x, na.rm = T))
}
clinicPwd = './data_clinical_patient.txt'
clinicPwd = './gbm_tcga_data_clinical_patient.txt'
clinicDf = data.frame(fread(clinicPwd, skip = 4))
colnames(clinicDf)
idhwtListPwd = './idhwtGBMList_cell2016_240802.csv'
idhwtListDf = data.frame(fread(idhwtListPwd))

unique(idhwtListDf$Study)
idhwtListDf = idhwtListDf[idhwtListDf$Study == "Glioblastoma multiforme", ]
idhwtListDf = idhwtListDf[idhwtListDf$IDH.status == 'WT', ]
idhwtPatientList = as.character(na.omit(idhwtListDf$Case))
length(idhwtPatientList)
colnames(clinicDf)
survivalDf = clinicDf[, c('PATIENT_ID', 'OS_STATUS', 'OS_MONTHS', 'DFS_STATUS', 'DFS_MONTHS')]

survivalDf[survivalDf == '[Not Available]'] = NA

survivalDf = survivalDf[survivalDf$PATIENT_ID %in% allIdhwtSampleIdArr, ]
dim(survivalDf)
rawRnaPwd = '/Volumes/Expansion/download/1_data/gbm_tcga/data_mrna_affymetrix_microarray.txt' # 12042 probes for 530 samples
rawRnaDf = data.frame(fread(rawRnaPwd))
dim(rawRnaDf)

boxplot(rawRnaDf[, 3:10]) # dist for patient
boxplot(t(rawRnaDf[1:10, -c(1:2)])) # dist for genes
rnaDfOnly = rawRnaDf[, -c(1:2)]
colnames(rnaDfOnly) = str_sub(gsub('[.]', '-', colnames(rnaDfOnly)), 1, 12)
geneNameList = rawRnaDf$Hugo_Symbol
geneNameList[is.na(geneNameList)] = ''

rnaDfOnlyScale = data.frame(normalize.quantiles(as.matrix(rnaDfOnly)))
colnames(rnaDfOnlyScale) = colnames(rnaDfOnly)

rnaDfOnlyScale = rnaDfOnlyScale[!duplicated(geneNameList), ]
rownames(rnaDfOnlyScale) = geneNameList[!duplicated(geneNameList)]

rnaDfOnlyScale = rnaDfOnlyScale[rownames(rnaDfOnlyScale) != '', ]

dim(rnaDfOnlyScale)

boxplot(rnaDfOnlyScale[, 3:10]) # all samples in same distribution
boxplot(t(rnaDfOnlyScale[3:10, ]))
colnames(rnaDfOnlyScale)
sum(c(l2largeSampleArr, l2smallSampleArr) %in% colnames(rnaDfOnlyScale))
sum(sgbmSampleArr %in% colnames(rnaDfOnlyScale))
rnaDfOnlyScaleSubset = rnaDfOnlyScale[, intersect(c(l2largeSampleArr, l2smallSampleArr), colnames(rnaDfOnlyScale))]
dim(rnaDfOnlyScaleSubset)
dim(rnaDfOnlyScale)
b2geneSetPwd = './resCor4geneExpressionWithMrB2_240920.csv'
b2geneSetDf = data.frame(fread(b2geneSetPwd))

b2geneSetDf = b2geneSetDf[order(b2geneSetDf$corPvalue), ]

b2geneSetArr = b2geneSetDf$geneName[b2geneSetDf$col == 'red']
length(b2geneSetArr)
migrTop50geneSetPwd = './resCor4geneExpressionWithMr_top50_qnorm_241030.csv' # qnorm

migrTop50geneSetDf = data.frame(fread(migrTop50geneSetPwd))

migrTop50geneSetArr = migrTop50geneSetDf$geneName[migrTop50geneSetDf$col == 'red']
migrTop50geneSetArrNeg = migrTop50geneSetDf$geneName[migrTop50geneSetDf$col == 'blue']

length(migrTop50geneSetArr)
length(migrTop50geneSetArrNeg)
rnaDfOnlyScaleSubset
rnaDfOnlyScaleSubset4score = rnaDfOnlyScaleSubset[rownames(rnaDfOnlyScaleSubset) %in% migrTop50geneSetArr, ]

resScoreMeanArr = apply(rnaDfOnlyScaleSubset4score, 2, mean)
ssgsea = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
  row_names = rownames(X)
  num_genes = nrow(X)
  gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})
  
  R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')
  
  es = apply(R, 2, function(R_col) {
    gene_ranks = order(R_col, decreasing = TRUE)
    
    es_sample = sapply(gene_sets, function(gene_set_idx) {
      indicator_pos = gene_ranks %in% gene_set_idx
      indicator_neg = !indicator_pos
      
      rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha
      
      step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
      step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)
      
      step_cdf_diff = step_cdf_pos - step_cdf_neg
      
      if (scale) step_cdf_diff = step_cdf_diff / num_genes
      
      if (single) {
        sum(step_cdf_diff)
      } else {
        step_cdf_diff[which.max(abs(step_cdf_diff))]
      }
    })
    unlist(es_sample)
  })
  
  if (length(gene_sets) == 1) es = matrix(es, nrow = 1)
  
  if (norm) es = es / diff(range(es))
  
  rownames(es) = names(gene_sets)
  colnames(es) = colnames(X)
  return(es)
}

geneSets = list(b2geneSet = b2geneSetArr,
                migrPosGeneSet = migrTop50geneSetArr,
                migrNegGeneSet = migrTop50geneSetArrNeg)

sum(rownames(rnaDfOnlyScaleSubset) %in% geneSets$migrPosGeneSet) 
sum(rownames(rnaDfOnlyScaleSubset) %in% geneSets$migrNegGeneSet)
resSsgsea = ssgsea(X = as.matrix(rnaDfOnlyScaleSubset), gene_sets = geneSets)

resSsgseaDf = data.frame(t(resSsgsea))
resSsgseaDf$sample = rownames(resSsgseaDf)

resSsgseaDf = resSsgseaDf[order(resSsgseaDf$migrPosGeneSet, decreasing = F), ]
resSsgseaDf$sample = factor(resSsgseaDf$sample, levels = resSsgseaDf$sample)
rawRnaDf = resSsgseaDf
rawRnaDf$patientId = as.character(rawRnaDf$sample)
rawRnaDf$expression = rawRnaDf$migrPosGeneSet - rawRnaDf$migrNegGeneSet # by gsea
rawRnaDf = data.frame(resScoreMeanArr)
colnames(rawRnaDf) = 'expression'
rawRnaDf$patientId = rownames(rawRnaDf)
resSsgseaDf$meanScoreMigrTop50 = rawRnaDf$expression[match(resSsgseaDf$sample,
                                                           rawRnaDf$patientId)]
rawRnaDf$time = survivalDf$OS_MONTHS[match(rawRnaDf$patientId,
                                           survivalDf$PATIENT_ID)]
rawRnaDf$status = survivalDf$OS_STATUS[match(rawRnaDf$patientId,
                                             survivalDf$PATIENT_ID)]
rawRnaDf$time = as.numeric(rawRnaDf$time)
resExpressionDf4cutoff = rawRnaDf
{
  
  setCutoffRate = 0.5
  
  grp1Patient = resExpressionDf4cutoff$patientId[resExpressionDf4cutoff$expression > quantile(resExpressionDf4cutoff$expression, probs = c(setCutoffRate, 1-setCutoffRate))[1]]
  grp2Patient = resExpressionDf4cutoff$patientId[resExpressionDf4cutoff$expression <= quantile(resExpressionDf4cutoff$expression, probs = c(setCutoffRate, 1-setCutoffRate))[1]]
  
  twoGroupDf = resExpressionDf4cutoff
  twoGroupDf$state = 0
  unique(twoGroupDf$status)
  
  twoGroupDf$state[twoGroupDf$status=='0:LIVING'] = 1
  twoGroupDf$state[twoGroupDf$status=='1:DECEASED'] = 2
  
  twoGroupDf$group4survival[twoGroupDf$patientId %in% grp1Patient] = 'High'
  twoGroupDf$group4survival[twoGroupDf$patientId %in% grp2Patient] = 'Low'
  
  fit = survfit(Surv(time, state) ~ group4survival, data = twoGroupDf)
  
  diff = survdiff(Surv(time, state) ~ group4survival, data = twoGroupDf)
  pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE) # get the t
  
  
  p = ggsurvplot(
    fit, 
    data = twoGroupDf, 
    size = 1,                 # change line size
    palette = "Set1",
    conf.int = F,          # Add confidence interval
    pval = T,              # Add p-value
    risk.table = T,        # Add risk table
    risk.table.col = "strata",# Risk table color by groups
    legend.labs = 
      c("High", "Low"),    # Change legend labels
    risk.table.height = 0.25, # Useful to change when you have multiple groups
    ggtheme = theme_classic() + 
      theme(
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)
      )     # Change ggplot2 theme
  )
  
  print(p)
  
  
  }
ggsurvplot(
  fit, 
  data = twoGroupDf, 
  size = 0.5,                 # change line size
  palette = "Set1",
  conf.int = F,          # Add confidence interval
  pval = F,              # Add p-value
  risk.table = F,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs = 
    c("High", "Low"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_classic() + 
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )     # Change ggplot2 theme
)
rawRnaDf2 = resSsgseaDf
rawRnaDf2$patientId = as.character(rawRnaDf2$sample)
rawRnaDf2$expression = rawRnaDf2$migrPosGeneSet - rawRnaDf2$migrNegGeneSet

rawRnaDf2$time = survivalDf$DFS_MONTHS[match(rawRnaDf2$patientId,
                                             survivalDf$PATIENT_ID)]
rawRnaDf2$status = survivalDf$DFS_STATUS[match(rawRnaDf2$patientId,
                                               survivalDf$PATIENT_ID)]
rawRnaDf2$time = as.numeric(rawRnaDf2$time)

resExpressionDf4cutoff = rawRnaDf2

{
  
  
  setCutoffRate = 0.5
  
  grp1Patient = resExpressionDf4cutoff$patientId[resExpressionDf4cutoff$expression > quantile(resExpressionDf4cutoff$expression, probs = c(setCutoffRate, 1-setCutoffRate))[1]]
  grp2Patient = resExpressionDf4cutoff$patientId[resExpressionDf4cutoff$expression <= quantile(resExpressionDf4cutoff$expression, probs = c(setCutoffRate, 1-setCutoffRate))[1]]
  
  
  twoGroupDf = resExpressionDf4cutoff
  twoGroupDf$state = 0
  unique(twoGroupDf$status)
  
  twoGroupDf$state[twoGroupDf$status =='0:DiseaseFree'] = 1
  twoGroupDf$state[twoGroupDf$status =='1:Recurred/Progressed'] = 2
  
  twoGroupDf$group4survival[twoGroupDf$patientId %in% grp1Patient] = 'High'
  twoGroupDf$group4survival[twoGroupDf$patientId %in% grp2Patient] = 'Low'
  
  fit = survfit(Surv(time, state) ~ group4survival, data = twoGroupDf)
  
  diff = survdiff(Surv(time, state) ~ group4survival, data = twoGroupDf)
  pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE) # get the t
  
  
  p = ggsurvplot(
    fit, 
    data = twoGroupDf, 
    size = 1,                 # change line size
    palette = "Set1",
    conf.int = F,          # Add confidence interval
    pval = T,              # Add p-value
    risk.table = T,        # Add risk table
    risk.table.col = "strata",# Risk table color by groups
    legend.labs = 
      c("High", "Low"),    # Change legend labels
    risk.table.height = 0.25, # Useful to change when you have multiple groups
    ggtheme = theme_classic() + 
      theme(
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)
      )     # Change ggplot2 theme
  )
  
  print(p)
  
  
  
  
}
ggsurvplot(
  fit, 
  data = twoGroupDf, 
  size = 0.5,                 # change line size
  palette = "Set1",
  conf.int = F,          # Add confidence interval
  pval = F,              # Add p-value
  risk.table = F,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs = 
    c("High", "Low"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_classic() + 
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )     # Change ggplot2 theme
)
