

library(reshape2)
library(data.table)
library(ggplot2)
library(ggthemes)
library(ggsci)
realTumorMelt = data.frame(fread('./240105_volumeByLesion.csv'))
colnames(realTumorMelt) = c('patient', 'lesion', 'realVolume')

setFillColor = rev(c('#953735', '#A76448', '#D6CA64',
                     '#64A6EE', '#94C0F3', '#214992'))

totalVolume = aggregate(realVolume ~ patient, data = realTumorMelt, sum)
orderedPatients = totalVolume$patient[order(totalVolume$realVolume)]
realTumorMelt$patient = factor(realTumorMelt$patient, levels = orderedPatients)
ggplot(realTumorMelt, aes(fill=factor(lesion, levels = paste0('L',6:1)), 
                          y=realVolume/1000, 
                          x=patient
)) + 
  geom_bar(position="stack", stat="identity", width = 0.85) +
  theme_classic() +
  scale_y_continuous(expand = expansion(add = c(0, 5e3/1e3))) +
  theme(text = element_text(size=15),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.justification = c(1, 1),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()
  ) +
  scale_fill_manual(values = setFillColor) + 
  labs(fill = 'Lesions',
       x = 'Patient ID',
       y = 'Tumor volume') +
  coord_flip()

