
library(reshape2)
library(data.table)
library(ggplot2)
library(ggthemes)
library(ggsci)
library(scales)
realTumorPwd = './realVolume_240412.csv'
realTumorDf = data.frame(fread(realTumorPwd))

ggplot(realTumorDf, 
       aes(fill =  factor(realTumorDf$nTumor, levels = 6:1), 
           x = nTumor)) + 
  geom_bar(position="stack", width = 0.6) +
  scale_fill_manual(values = rev(c('#36BA98', '#E9C46A', '#F4A261', '#E76F51', '#C40C0C'))) +
  scale_y_continuous(expand=c(0, 0), 
                     limits = c(0, 21)) +
  theme_classic() + # coord_flip() +
  theme(text = element_text(size = 15)) +
  labs(fill = '# of total lesions',
       x = '# of total lesions',
       y = '# of patients') +
  guides(fill = 'none')

