
library(ggbeeswarm)
library(ggplot2)
library(data.table)
library(dplyr)
library(progress)

rm(list = ls())
load('./validationCurve4lesion3cutoff0.RData')

resAllRealTumorNumCount4SimCombine = readRDS('./sfig4A_240613.RDS')
colnames(resAllRealTumorNumCount4SimCombine)
table(resAllRealTumorNumCount4SimCombine$label)

ggplot(resAllRealTumorNumCount4SimCombine[resAllRealTumorNumCount4SimCombine$label == as.character(setLesionNum4validation),], 
       mapping = aes(x = log10(cutoff+1),
                     y = Total,
                     group = patientLabel,
                     col = patientLabel
       )) +
  geom_line(alpha=0.8) +
  geom_point(alpha = 0.8, size = 0.5) +
  theme_classic() + 
  theme(text = element_text(size=15),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  labs(x = 'log10(Cutoff cell number for invisible tumor)',
       y = paste0('Prob( only ', setLesionNum4validation, ' lesions | ', setLesionNum4validation, ' visible lesions)'),
       col = 'Virtual patient group',
       title = paste0('Wilcox test p.value = ', signif(wilcox.test(Total ~ patientLabel, data = resAllRealTumorNumCount4SimCombine[resAllRealTumorNumCount4SimCombine$label == as.character(setLesionNum4validation),])$p.value, digits = 4))) +
  guides(col = 'none')

