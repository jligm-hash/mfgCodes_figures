

library(data.table)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)

lesionCenterPwd = './230914_lesionMassCenter221228.csv'
lesionCenterDf = data.frame(fread(lesionCenterPwd))
patientList = unique(lesionCenterDf$Patient)
resL1L2distList = c()
for (tmpPatientId in patientList) {
  
  tmpPatientDf = lesionCenterDf[lesionCenterDf$Patient == tmpPatientId, ]
  
  tmpL1Df = tmpPatientDf[tmpPatientDf$Lesion.ID == 'L1',]
  tmpL2Df = tmpPatientDf[tmpPatientDf$Lesion.ID == 'L2',]
  
  tmpDistance = sqrt((tmpL1Df$diagnostics_Mask.original_CenterOfMass_X - tmpL2Df$diagnostics_Mask.original_CenterOfMass_X)^2 +
                       (tmpL1Df$diagnostics_Mask.original_CenterOfMass_Y - tmpL2Df$diagnostics_Mask.original_CenterOfMass_Y)^2 +
                       (tmpL1Df$diagnostics_Mask.original_CenterOfMass_Z - tmpL2Df$diagnostics_Mask.original_CenterOfMass_Z)^2 )
  
  resL1L2distList = c(resL1L2distList, tmpDistance)
}
resL1L2distDf = data.frame(patientId = patientList,
                           l1l2dist = resL1L2distList)
groupLabelPwd = './mgbmSgbmLabel.csv'
groupLabelDf = data.frame(fread(groupLabelPwd))

resL1L2distDf$group = groupLabelDf$bigSmall[match(resL1L2distDf$patientId, groupLabelDf$ID)]
lesionNumberPwd = './realTCGA211017.csv'
lesionNumDf = data.frame(fread(lesionNumberPwd))
resL1L2distDf$num = lesionNumDf$Num[match(resL1L2distDf$patientId, lesionNumDf$ID)]

table(resL1L2distDf[resL1L2distDf$num <= 3,]$group)
set.seed(100)
ggplot(resL1L2distDf[resL1L2distDf$num <= 3,], mapping = aes(x = group,
                                                             y = l1l2dist,
                                                             col = group)) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.05, size = 0.5, alpha=0.8) +
  theme_classic() +
  guides(col = 'none') +
  scale_color_brewer(palette = 'Set1') +
  stat_compare_means(comparisons = list(c(1,2)), 
                     method = 't.test') +
  theme(text = element_text(size = 15)) +
  labs(x = 'MGBM group',
       y = 'Distance between mass center of L1 and L2')

