
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
gema0 = read.delim(file.path(setFilePwd, 'gbm_tcga_data_mrna_affymetrix_microarray.txt'),
                   check.names = F)
gema = gema0[, names(gema0) %in% c('Hugo_Symbol',paste0(c(l2b,l2s),'-01'))]
colnames(gema0) = str_sub(colnames(gema0), 1, 12)
sgbmMgbmExpressionDf = gema0[, str_sub(colnames(gema0), 1, 12) %in% c('Hugo_Symbol', sgbmMgbmLabelDf$ID)]

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
setFillColor = rev(c('#9B4C4B', '#BC8D78', '#D6D2B0',
                     '#96C0EE', '#C1D8F3', '#485C81'))

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
expressionMat = gema0
geneNameList = expressionMat$Hugo_Symbol
expressionMat$Hugo_Symbol = NULL
expressionMat$Entrez_Gene_ = NULL
expressionMatScaled = data.frame(t(apply(expressionMat, 1, scale)))

colnames(expressionMatScaled) = colnames(expressionMat)
setGeneSet = c("UPP1", "CAPN2", "FOSL1", "DPP4", "MREG", "GADD45A", "TNFRSF12A", "ANG", "FER1L3", "HEBP1", "CNIH3", "SERPINE1", "CA9", "NQO2", "IGFBP6", "SNX10", "MICALL2", "OSBPL3", "CCDC109B", "PI3", "CD44", "ZCCHC6", "BCAT1", "CLEC5A", "MCFD2", "GADD45B", "LOC57228", "ANXA2", "tcag7.1314", "TAGLN2")

expressionMatScaledSubset = expressionMatScaled[geneNameList %in% setGeneSet,
]

expressionMatScaledSubset = expressionMatScaledSubset[, colnames(expressionMatScaledSubset) %in% sgbmMgbmLabelDf$ID]
tmpMeanArr = apply(expressionMatScaledSubset, 2, mean)
tmpVarArr = apply(expressionMatScaledSubset, 2, var)

resGeneSetExpressionDf = data.frame(patientId = names(tmpMeanArr),
                                    mean = tmpMeanArr,
                                    var = tmpVarArr)
resGeneSetExpressionDf$mriGroup = sgbmMgbmLabelDf$bigSmall[match(resGeneSetExpressionDf$patientId, sgbmMgbmLabelDf$ID)]

resGeneSetExpressionDf$mriGroup = factor(resGeneSetExpressionDf$mriGroup,
                                         levels = c('Small', 'SGBM', 'Big'))

table(resGeneSetExpressionDf$mriGroup)
resGeneSetExpressionDf$time = survivalDf$time[match(resGeneSetExpressionDf$patientId, survivalDf$patient)]
resGeneSetExpressionDf$status = survivalDf$status[match(resGeneSetExpressionDf$patientId, survivalDf$patient)]
colnames(resGeneSetExpressionDf)
resExpressionDf4cutoff = resGeneSetExpressionDf[resGeneSetExpressionDf$mriGroup == 'SGBM', ]

resExpressionDf4cutoff = resExpressionDf4cutoff[resExpressionDf4cutoff$time <= 1000, ]

resExpressionDf4cutoff$time = resExpressionDf4cutoff$time/30
{
  
  
  
  grp1Patient = resExpressionDf4cutoff$patientId[resExpressionDf4cutoff$mean > quantile(resExpressionDf4cutoff$mean, probs = c(0.33, 0.67))[2]]
  grp2Patient = resExpressionDf4cutoff$patientId[resExpressionDf4cutoff$mean <= quantile(resExpressionDf4cutoff$mean, probs = c(0.33, 0.67))[1]]
  
  grp1Df = resExpressionDf4cutoff[resExpressionDf4cutoff$patientId %in% grp1Patient,]
  grp2Df = resExpressionDf4cutoff[resExpressionDf4cutoff$patientId %in% grp2Patient,]
  
  twoGroupDf = rbind(grp1Df,
                     grp2Df)
  twoGroupDf$state = 0
  twoGroupDf$state[twoGroupDf$status=='Alive'] = 1
  twoGroupDf$state[twoGroupDf$status=='Dead'] = 2
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
newClinicalGbmAll = readRDS('./processedTcgaDfs240329.rds')
newClinicalMgbm = newClinicalGbmAll
mgbmDfNew = data.frame('patient' = newClinicalMgbm$PATIENT_ID,
                       "time" = as.numeric(as.character(newClinicalMgbm$DFS_MONTHS)), 
                       "status" = as.numeric(as.character(newClinicalMgbm$DFS_STATUS)))
survivalDf = unique(mgbmDfNew)
survivalDf$time = as.numeric(survivalDf$time)

resGeneSetExpressionDf = data.frame(patientId = names(tmpMeanArr),
                                    mean = tmpMeanArr,
                                    var = tmpVarArr)

resGeneSetExpressionDf$mriGroup = sgbmMgbmLabelDf$bigSmall[match(resGeneSetExpressionDf$patientId, sgbmMgbmLabelDf$ID)]

resGeneSetExpressionDf$mriGroup = factor(resGeneSetExpressionDf$mriGroup,
                                         levels = c('Small', 'SGBM', 'Big'))

table(resGeneSetExpressionDf$mriGroup)
resGeneSetExpressionDf$time = survivalDf$time[match(resGeneSetExpressionDf$patientId, survivalDf$patient)]
resGeneSetExpressionDf$status = survivalDf$status[match(resGeneSetExpressionDf$patientId, survivalDf$patient)]
colnames(resGeneSetExpressionDf)
resExpressionDf4cutoff = resGeneSetExpressionDf[resGeneSetExpressionDf$mriGroup == 'SGBM', ]
resExpressionDf4cutoff = na.omit(resExpressionDf4cutoff)

{
  grp1Patient = resExpressionDf4cutoff$patientId[resExpressionDf4cutoff$mean > quantile(resExpressionDf4cutoff$mean, probs = c(0.33, 0.67))[2]]
  grp2Patient = resExpressionDf4cutoff$patientId[resExpressionDf4cutoff$mean <= quantile(resExpressionDf4cutoff$mean, probs = c(0.33, 0.67))[1]]
  
  grp1Df = resExpressionDf4cutoff[resExpressionDf4cutoff$patientId %in% grp1Patient,]
  grp2Df = resExpressionDf4cutoff[resExpressionDf4cutoff$patientId %in% grp2Patient,]
  
  twoGroupDf = rbind(grp1Df,
                     grp2Df)
  twoGroupDf$state = twoGroupDf$status
  
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
    pval = TRUE,              # Add p-value
    risk.table = TRUE,        # Add risk table
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

