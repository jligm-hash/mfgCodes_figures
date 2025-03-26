
library(DT)
library(survminer)
library(survival)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(ggsci)
library(ggrepel)

myMax = function(x){
  return(max(as.numeric(x), na.rm = T))
}
realTumor = data.frame(fread('./240105_volumeByLesion.csv'))
realTumor = data.frame(reshape2::dcast(realTumor, patient ~ lesion))
realTumor$patient = realTumor$patient
realTumor[, paste0('T', 1:6)] = realTumor[, paste0('L', 1:6)]

realTumor[is.na(realTumor)] = 0
newClinicalGbmAll = fread('./clinical.tsv')
newClinicalMgbm = newClinicalGbmAll[case_submitter_id %in% realTumor$patient]
groupGDC_clinicalNew = newClinicalMgbm[,c('case_submitter_id','days_to_death','days_to_last_follow_up','vital_status')]

groupGDC_clinicalNew = data.frame(groupGDC_clinicalNew)
groupGDC_clinicalNew[groupGDC_clinicalNew == "'--"] = 0
mgbm_t3New = apply(groupGDC_clinicalNew[,c('days_to_death','days_to_last_follow_up')],1,myMax)
mgbmDfNew = data.frame('patient' = groupGDC_clinicalNew$case_submitter_id,
                       "time" = mgbm_t3New, 
                       "gbmType" = "group1Mgbm", 
                       "status" = groupGDC_clinicalNew$vital_status)
survivalDf = unique(mgbmDfNew)
survivalDf$time = as.numeric(survivalDf$time)

orderList = c()
cutoffSize = c()
pvalue = c()
for (i in 1:(length(realTumor$T2)-1)){
  
  
  sizeCutoff = realTumor$T2[order(realTumor$T2)][i]
  grp1Patient = realTumor$patient[realTumor$T2 > sizeCutoff]
  grp2Patient = realTumor$patient[realTumor$T2 <= sizeCutoff]
  grp1Df = survivalDf[survivalDf$patient %in% grp1Patient,]
  grp2Df = survivalDf[survivalDf$patient %in% grp2Patient,]
  
  twoGroupDf = rbind(grp1Df,
                     grp2Df)
  twoGroupDf$state = 0
  twoGroupDf$state[twoGroupDf$status=='Alive'] = 0
  twoGroupDf$state[twoGroupDf$status=='Dead'] = 1
  twoGroupDf$gbmType[twoGroupDf$patient %in% grp1Patient] = 'grp1'
  twoGroupDf$gbmType[twoGroupDf$patient %in% grp2Patient] = 'grp2'
  
  diff = survdiff(Surv(time, state) ~ gbmType, data = twoGroupDf, rho = 1)
  tmpPvalue = diff$pvalue # get the t
  
  
  orderList = c(orderList, i)
  cutoffSize = c(cutoffSize, sizeCutoff)
  pvalue = c(pvalue, tmpPvalue)
  
}

pvalueDf = data.frame(orderList = orderList,
                      cutoffSize = cutoffSize,
                      pvalue = pvalue)
orderList2 = c()
cutoffSize2 = c()
pvalue2 = c()
for (i in 1:(length(realTumor$T1)-1)){
  sizeCutoff = realTumor$T1[order(realTumor$T1)][i]
  grp1Patient = realTumor$patient[realTumor$T1 > sizeCutoff]
  grp2Patient = realTumor$patient[realTumor$T1 <= sizeCutoff]
  grp1Df = survivalDf[survivalDf$patient %in% grp1Patient,]
  grp2Df = survivalDf[survivalDf$patient %in% grp2Patient,]
  
  twoGroupDf = rbind(grp1Df,
                     grp2Df)
  twoGroupDf$state = 0
  twoGroupDf$state[twoGroupDf$status=='Alive'] = 0
  twoGroupDf$state[twoGroupDf$status=='Dead'] = 1
  twoGroupDf$gbmType[twoGroupDf$patient %in% grp1Patient] = 'grp1'
  twoGroupDf$gbmType[twoGroupDf$patient %in% grp2Patient] = 'grp2'
  
  diff = survdiff(Surv(time, state) ~ gbmType, data = twoGroupDf)
  tmpPvalue = diff$pvalue # get the t
  
  orderList2 = c(orderList2, i)
  cutoffSize2 = c(cutoffSize2, sizeCutoff)
  pvalue2 = c(pvalue2, tmpPvalue)
  
}

pvalueDf2 = data.frame(orderList = orderList2,
                       cutoffSize = cutoffSize2,
                       pvalue = pvalue2)
orderList3 = c()
cutoffSize3 = c()
pvalue3 = c()
for (i in 1:(length(realTumor$T2)-1)){
  sizeCutoff = realTumor$T3[order(realTumor$T3)][i]
  grp1Patient = realTumor$patient[realTumor$T3 > sizeCutoff]
  grp2Patient = realTumor$patient[realTumor$T3 <= sizeCutoff]
  grp1Df = survivalDf[survivalDf$patient %in% grp1Patient,]
  grp2Df = survivalDf[survivalDf$patient %in% grp2Patient,]
  
  twoGroupDf = rbind(grp1Df,
                     grp2Df)
  twoGroupDf$state = 0
  twoGroupDf$state[twoGroupDf$status=='Alive'] = 0
  twoGroupDf$state[twoGroupDf$status=='Dead'] = 1
  twoGroupDf$gbmType[twoGroupDf$patient %in% grp1Patient] = 'grp1'
  twoGroupDf$gbmType[twoGroupDf$patient %in% grp2Patient] = 'grp2'
  
  diff = survdiff(Surv(time, state) ~ gbmType, data = twoGroupDf)
  tmpPvalue = diff$pvalue # get the t
  
  orderList3 = c(orderList3, i)
  cutoffSize3 = c(cutoffSize3, sizeCutoff)
  pvalue3 = c(pvalue3, tmpPvalue)
  
}

pvalueDf3 = data.frame(orderList = orderList3,
                       cutoffSize = cutoffSize3,
                       pvalue = pvalue3)

pvalueDf$group = 'L2'
pvalueDf2$group = 'L1'
pvalueDf3$group = 'L3'
pvalueDfCombine = rbind(pvalueDf, pvalueDf2, pvalueDf3)
colnames(pvalueDf)

pvalueDfCombine$group = factor(pvalueDfCombine$group,
                               levels = paste0("L", c(1,3,2)))
library(grid)

ggplot(pvalueDfCombine[pvalueDfCombine$orderList< 34, ], mapping = aes(x = orderList/35,
                                                                       y = as.numeric(pvalue),
                                                                       group=group,
                                                                       color=group
)) +
  
  geom_hline(yintercept = 0.05, color = 'black', alpha = 1, linetype = "dotted")+
  geom_line(alpha = 0.7) +
  geom_point(alpha = 0.6, size = 0.8) + 
  labs(x = 'Quantile by tumor volume',
       y = 'P-value of OS differences') +
  scale_color_manual(values = c('#953735', '#D6CA64', '#A76448')) +
  scale_y_continuous(limits = c(1e-3, 1), transform = 'log10') +
  scale_x_continuous(limits = c(0, 1)) +
  labs(color="") +
  theme_classic() +
  theme(text = element_text(size = 15),
  )

