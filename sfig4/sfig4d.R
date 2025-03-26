
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
cggaPwd = '/Volumes/Expansion/download/1_data/'

cgga325clinicalPwd = file.path(cggaPwd, 'CGGA.mRNAseq_325_clinical.20200506.txt')
cgga325rnaPwd = file.path(cggaPwd, 'CGGA.mRNAseq_325.RSEM-genes.20200506.txt')

cgga693clinicalPwd = file.path(cggaPwd, 'CGGA.mRNAseq_693_clinical.20200506.txt')
cgga693rnaPwd = file.path(cggaPwd, 'CGGA.mRNAseq_693.RSEM-genes.20200506.txt')
cgga325clinicalDf = data.frame(fread(cgga325clinicalPwd))
cgga325rnaDf = data.frame(fread(cgga325rnaPwd))

dim(cgga325clinicalDf)
dim(cgga325rnaDf)

cgga693clinicalDf = data.frame(fread(cgga693clinicalPwd))
cgga693rnaDf = data.frame(fread(cgga693rnaPwd))

dim(cgga693clinicalDf)
dim(cgga693rnaDf)
cgga325clinicalDf_idhwt = cgga325clinicalDf[cgga325clinicalDf$IDH_mutation_status == 'Wildtype', ]
cgga325clinicalDf_idhwt = cgga325clinicalDf_idhwt[!is.na(cgga325clinicalDf_idhwt$IDH_mutation_status), ]

cgga325clinicalDf_idhwt = cgga325clinicalDf_idhwt[cgga325clinicalDf_idhwt$Histology == 'GBM', ]
cgga325clinicalDf_idhwt = cgga325clinicalDf_idhwt[!is.na(cgga325clinicalDf_idhwt$Histology), ]
table(cgga325clinicalDf_idhwt$IDH_mutation_status)

cgga325rnaDf_idhwt = cgga325rnaDf[, c('Gene_Name', cgga325clinicalDf_idhwt$CGGA_ID)]
dim(cgga325rnaDf_idhwt)
cgga693clinicalDf_idhwt = cgga693clinicalDf[cgga693clinicalDf$IDH_mutation_status == 'Wildtype', ]
cgga693clinicalDf_idhwt = cgga693clinicalDf_idhwt[!is.na(cgga693clinicalDf_idhwt$IDH_mutation_status), ]

cgga693clinicalDf_idhwt = cgga693clinicalDf_idhwt[cgga693clinicalDf_idhwt$Histology == 'GBM', ]
cgga693clinicalDf_idhwt = cgga693clinicalDf_idhwt[!is.na(cgga693clinicalDf_idhwt$Histology), ]
table(cgga693clinicalDf_idhwt$IDH_mutation_status)

cgga693rnaDf_idhwt = cgga693rnaDf[, c('Gene_Name', cgga693clinicalDf_idhwt$CGGA_ID)]
dim(cgga693rnaDf_idhwt)
sum(duplicated(c(cgga325clinicalDf_idhwt$CGGA_ID,
                 cgga693clinicalDf_idhwt$CGGA_ID)))

length(unique(c(cgga325clinicalDf_idhwt$CGGA_ID,
                cgga693clinicalDf_idhwt$CGGA_ID)))
setGeneSet = c("UPP1", "CAPN2", "FOSL1", "DPP4", "MREG", "GADD45A", "TNFRSF12A", "ANG", "FER1L3", "HEBP1", "CNIH3", "SERPINE1", "CA9", "NQO2", "IGFBP6", "SNX10", "MICALL2", "OSBPL3", "CCDC109B", "PI3", "CD44", "ZCCHC6", "BCAT1", "CLEC5A", "MCFD2", "GADD45B", "LOC57228", "ANXA2", "tcag7.1314", "TAGLN2")
cgga325rnaDfsubset = cgga325rnaDf_idhwt[cgga325rnaDf_idhwt$Gene_Name %in% setGeneSet, ]
cgga693rnaDfsubset = cgga693rnaDf_idhwt[cgga693rnaDf_idhwt$Gene_Name %in% setGeneSet, ]
cggaMergeRnaDf = cbind(cgga325rnaDfsubset[, -1], 
                       cgga693rnaDfsubset[, -1])
cggaMergeRnaDf = data.frame(apply(cggaMergeRnaDf, 2, as.numeric))

cggaMergeRnaDf = cggaMergeRnaDf[, colnames(cggaMergeRnaDf) %in% c(cgga325clinicalDf_idhwt$CGGA_ID[cgga325clinicalDf_idhwt$PRS_type == 'Primary'], 
                                                                  cgga693clinicalDf_idhwt$CGGA_ID[cgga693clinicalDf_idhwt$PRS_type == 'Primary'] )]
cggaMergeRnaDfScaled = cggaMergeRnaDf
cggaMergeRnaDfScaled = data.frame(t(apply(cggaMergeRnaDfScaled, 1, scale)))

colnames(cggaMergeRnaDfScaled) = colnames(cggaMergeRnaDf)

resMgbScore4cggaArrMean = apply(cggaMergeRnaDfScaled, 2, mean)
resMgbScore4cggaArrVar = apply(cggaMergeRnaDfScaled, 2, var)
resMarkerExpressionDf = data.frame(patientId = colnames(cggaMergeRnaDfScaled),
                                   mean = resMgbScore4cggaArrMean,
                                   var = resMgbScore4cggaArrVar)

resMarkerExpressionDf$time = coalesce(cgga325clinicalDf_idhwt$OS[match(resMarkerExpressionDf$patientId, cgga325clinicalDf_idhwt$CGGA_ID)],
                                      cgga693clinicalDf_idhwt$OS[match(resMarkerExpressionDf$patientId, cgga693clinicalDf_idhwt$CGGA_ID)])
resMarkerExpressionDf$state = coalesce(cgga325clinicalDf_idhwt$Censor..alive.0..dead.1.[match(resMarkerExpressionDf$patientId, cgga325clinicalDf_idhwt$CGGA_ID)],
                                       cgga693clinicalDf_idhwt$Censor..alive.0..dead.1.[match(resMarkerExpressionDf$patientId, cgga693clinicalDf_idhwt$CGGA_ID)])
resExpressionDf4cutoff = resMarkerExpressionDf

resExpressionDf4cutoff = na.omit(resExpressionDf4cutoff)
resExpressionDf4cutoff = resExpressionDf4cutoff[resExpressionDf4cutoff$time < 2000, ]
resExpressionDf4cutoff$time = resExpressionDf4cutoff$time/30
{
  grp1Patient = resExpressionDf4cutoff$patientId[resExpressionDf4cutoff$mean > quantile(resExpressionDf4cutoff$mean, probs = c(0.2, 0.75))[2]]
  grp2Patient = resExpressionDf4cutoff$patientId[resExpressionDf4cutoff$mean <= quantile(resExpressionDf4cutoff$mean, probs = c(0.2, 0.75))[1]]
  
  
  grp1Df = resExpressionDf4cutoff[resExpressionDf4cutoff$patientId %in% grp1Patient,]
  grp2Df = resExpressionDf4cutoff[resExpressionDf4cutoff$patientId %in% grp2Patient,]
  
  
  
  twoGroupDf = rbind(grp1Df,
                     grp2Df)
  twoGroupDf$group4survival[twoGroupDf$patientId %in% grp1Patient] = 'High'
  twoGroupDf$group4survival[twoGroupDf$patientId %in% grp2Patient] = 'Low'
  
  
  fit = survfit(Surv(time, state) ~ group4survival, data = twoGroupDf)
  
  diff = survdiff(Surv(time, state) ~ group4survival, data = twoGroupDf)
  pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE) # get the t
  
  
  ggsurvplot(
    fit, 
    data = twoGroupDf, 
    size = 1,                 # change line size
    palette = "Set1", 
    conf.int = F,          # Add confidence interval
    pval = T,              # Add p-value
    risk.table = F,        # Add risk table
    risk.table.col = "strata",# Risk table color by groups
    legend.labs = 
      c("High", "Low"),    # Change legend labels
    risk.table.height = 0.25, # Useful to change when you have multiple groups
    ggtheme = theme_classic() + 
      theme(axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12))     # Change ggplot2 theme
  )
  
}

table(twoGroupDf$group4survival)

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
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))     # Change ggplot2 theme
)

