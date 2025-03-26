
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(reshape2)
library(data.table)
library(ggthemes)
library(ggsci)

library(DT)
library(survminer)
library(survival)
myMax = function(x){
  return(max(x, na.rm = T))
}
newClinicalGbmAll = fread('./clinical.project-TCGA-GBM.2022-01-18//clinical.tsv')

newClinicalMgbm = newClinicalGbmAll

groupGDC_clinicalNew = newClinicalMgbm[,c('case_submitter_id','days_to_death','days_to_last_follow_up','vital_status')]
mgbm_t3New = apply(groupGDC_clinicalNew[,c('days_to_death','days_to_last_follow_up')],1,myMax)
mgbmDfNew = data.frame('patient' = groupGDC_clinicalNew$case_submitter_id,
                       "time" = mgbm_t3New, 
                       "status" = groupGDC_clinicalNew$vital_status)
survivalDf = unique(mgbmDfNew)
survivalDf$time = as.numeric(survivalDf$time)

survivalDf$time = survivalDf$time/30 # day to mon
multicentricLabelPwd = './zhangWei2015_multicentricLabel30tcga.csv'

multicentricDf = data.frame(fread(multicentricLabelPwd))
grp1Patient = multicentricDf$patientId[multicentricDf$mri.classification == 'multifocal']
grp2Patient = multicentricDf$patientId[multicentricDf$mri.classification == 'multicentric']
{
  
  
  
  
  
  grp1Df = survivalDf[survivalDf$patient %in% grp1Patient, ]
  grp2Df = survivalDf[survivalDf$patient %in% grp2Patient, ]
  
  
  
  
  
  twoGroupDf = rbind(grp1Df,
                     grp2Df)
  twoGroupDf$state = 0
  twoGroupDf$state[twoGroupDf$status=='Alive'] = 1
  twoGroupDf$state[twoGroupDf$status=='Dead'] = 2
  
  twoGroupDf$group4survival[twoGroupDf$patient %in% grp1Patient] = 'multifocal'
  twoGroupDf$group4survival[twoGroupDf$patient %in% grp2Patient] = 'multicentric'
  
  table(twoGroupDf$group4survival)
  
  
  fit = survfit(Surv(time, state) ~ group4survival, data = twoGroupDf)
  
  diff = survdiff(Surv(time, state) ~ group4survival, data = twoGroupDf)
  pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE) # get the t
  
  
  ggsurvplot(
    fit, 
    data = twoGroupDf, 
    size = 0.5,                 # change line size
    palette = c('#ff7f00', '#fb9a99'), #"Set1",
    conf.int = F,          # Add confidence interval
    pval = F,              # Add p-value
    risk.table = F,        # Add risk table
    risk.table.col = "strata",# Risk table color by groups
    legend.labs = 
      c( "Multicentric", "Multifocal"),    # Change legend labels
    risk.table.height = 0.25, # Useful to change when you have multiple groups
    ggtheme = theme_classic() + 
      theme(
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)
      )     # Change ggplot2 theme
  )
  
}

