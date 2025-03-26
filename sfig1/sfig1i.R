
library(data.table)

clinicalTable = fread('./data_clinical_patient.txt',
                      skip = 4)
colnames(clinicalTable)
library(data.table)
mgbmLabelPwd = './mgbmSgbmLabelMerge230907.csv'
mgbmLabelDf = data.frame(fread(mgbmLabelPwd))

table(mgbmLabelDf$bigSmall)

group1Patient = mgbmLabelDf$ID[mgbmLabelDf$bigSmall == 'Big']
group1Label = 'Big'
group2Patient = mgbmLabelDf$ID[mgbmLabelDf$bigSmall == 'Small']
group2Label = 'Small'
cat(group1Patient, sep = ',')
cat(group2Patient, sep = ',')

length(unique(c(group1Patient, group2Patient)))

targetKPS = clinicalTable[clinicalTable$PATIENT_ID %in% c(group1Patient, group2Patient), ]
targetKPS$group = 'unknown'
targetKPS$group[targetKPS$PATIENT_ID %in% group1Patient] = group1Label
targetKPS$group[targetKPS$PATIENT_ID %in% group2Patient] = group2Label

targetKPS$groupLabel = 'unknown'
targetKPS$groupLabel[targetKPS$group == group1Label] = paste0(group1Label, '=', sum(targetKPS$group == group1Label))
targetKPS$groupLabel[targetKPS$group == group2Label] = paste0(group2Label, '=', sum(targetKPS$group == group2Label))

targetKPS = targetKPS[targetKPS$KARNOFSKY_PERFORMANCE_SCORE != "[Not Available]", ]
targetKPS$KARNOFSKY_PERFORMANCE_SCORE = as.numeric(targetKPS$KARNOFSKY_PERFORMANCE_SCORE )

t.test(data=targetKPS,KARNOFSKY_PERFORMANCE_SCORE ~ group)$p.value
wilcox.test(data=targetKPS,KARNOFSKY_PERFORMANCE_SCORE ~ group)$p.value
pValue = round(wilcox.test(data=targetKPS,KARNOFSKY_PERFORMANCE_SCORE ~ group)$p.value, 
               digits = 4)

library(ggbeeswarm)
table(targetKPS$group)
ggplot(targetKPS, mapping = aes(x = group,
                                y = KARNOFSKY_PERFORMANCE_SCORE,
                                group = group,
                                col = groupLabel)) +
  geom_boxplot(width = 0.3) +
  geom_beeswarm(cex = 3, alpha = 0.8, size = 0.5) +
  scale_colour_manual(values=c("#E41A1C","#377EB8")) +
  guides(col = 'none') +
  theme_classic() +
  theme(text = element_text(size=20),
        plot.subtitle =element_text(hjust=0.5)) +
  labs(fill = '', x = '', y = 'Karnofsky Performance Score',
  )
wilcox.test(targetKPS$KARNOFSKY_PERFORMANCE_SCORE[targetKPS$group == 'Big'],
       targetKPS$KARNOFSKY_PERFORMANCE_SCORE[targetKPS$group == 'Small'])

