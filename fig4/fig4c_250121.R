
library(data.table)

resBestSimPwd2 = './resRecurrentTime4patient230906.csv'

resBestSimDf2 = data.frame(fread(resBestSimPwd2))
library(data.table)
resDfPwd = './resRecurrentTime4patient230829.csv'
resDf = data.frame(fread(resDfPwd))
table(resDf$resMigrArr)

resDf$recurrentDuration = resDf$resRecurrentTimeArr - resDf$resSurgeryTimeArr
resDf$recurrentDuration4visible = resDf$resRecurrentVisibleTime - resDf$resSurgeryTimeArr

resDf$group[resDf$resMigrArr < 1e-5] = 'Low'
resDf$group[resDf$resMigrArr >= 1e-5] = 'High'

resDfMelt = reshape2::melt(resDf)
colnames(resDfMelt) = c('Group', 'resName', 'Value')
library(data.table)
bestSimPwd = './resAllBestSimNew35.csv'

bestSimDf = data.frame(fread(bestSimPwd))
resDf4bigSmall = resDf
resDf4bigSmall$patientID = NA
for (i in 1:nrow(resDf4bigSmall)) {
  
  
  tmpBestSimDf = resDf4bigSmall[i,]
  
  tmpMigr4tmpBestSim = tmpBestSimDf$resMigrArr
  tmpSeed4tmpBestSim = tmpBestSimDf$resSeedArr
  
  resFilterMigr4bestSimDf = bestSimDf[bestSimDf$migr == tmpMigr4tmpBestSim/1e3,]
  resFilterSeed4bestSimDf = resFilterMigr4bestSimDf[resFilterMigr4bestSimDf$seed == tmpSeed4tmpBestSim, ]
  
  
  
  if(nrow(resFilterSeed4bestSimDf) == 1){
    resDf4bigSmall$patientID[i] = resFilterSeed4bestSimDf$p.ID
  }else{
    print(i)
    print('More than one occurence was observed, will save the first one only')
    resDf4bigSmall$patientID[i] = resFilterSeed4bestSimDf$p.ID[1]
  }
  
}
mgbmPatientLabelPwd = './mgbmSgbmLabel.csv'

mgbmPatientLabelDf = data.frame(fread(mgbmPatientLabelPwd))
mgbmBigSmallDf = mgbmPatientLabelDf[mgbmPatientLabelDf$bigSmall != 'SGBM',]

dim(mgbmBigSmallDf)

mgbmBigSmallDf$migr = bestSimDf$migr[match(mgbmBigSmallDf$ID, bestSimDf$p.ID)]
mgbmBigSmallDf$seed = bestSimDf$seed[match(mgbmBigSmallDf$ID, bestSimDf$p.ID)]

colnames(resDf4bigSmall)
for (i in 1:nrow(mgbmBigSmallDf)) {
  
  tmpSeed = mgbmBigSmallDf$seed[i]
  tmpMigr = mgbmBigSmallDf$migr[i]*1e3
  
  tmpIndex = resDf4bigSmall$resMigrArr == tmpMigr & resDf4bigSmall$resSeedArr == tmpSeed
  
  if(sum(tmpIndex) == 0){
    mgbmBigSmallDf[i, colnames(resDf4bigSmall)] = NA
    
  }else if(sum(tmpIndex) == 1){
    mgbmBigSmallDf[i, colnames(resDf4bigSmall)] = resDf4bigSmall[tmpIndex, ]
  }else{
    mgbmBigSmallDf[i, colnames(resDf4bigSmall)] = resDf4bigSmall[tmpIndex, ][1,]
  }  
  
}
library(RColorBrewer)
library(ggpubr)
library(ggbeeswarm)

table(mgbmBigSmallDf$bigSmall)

my_comparisons = list(c(1, 2))

mgbmBigSmallDf4plt = mgbmBigSmallDf[mgbmBigSmallDf$resTotalLesionArr > 1, ]
mgbmBigSmallDf4plt = na.omit(mgbmBigSmallDf4plt)

colnames(mgbmBigSmallDf4plt)

ggplot(mgbmBigSmallDf4plt, mapping = aes(x = bigSmall,
                                         y = recurrentDuration,
                                         col = bigSmall)) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
  geom_beeswarm(cex = 3, size = 0.5, alpha = 0.8) +
  theme_classic() +
  theme(text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = 'MGBM group',
       y = 'Recurrent time for symptom lesions (Generations)',
       title = 'Remove the largest tumor lesion') +
  scale_color_brewer(palette = 'Set1') +
  guides(col = 'none') +
  stat_compare_means(comparisons = my_comparisons) 

