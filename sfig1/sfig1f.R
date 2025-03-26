
library(data.table)
library(RColorBrewer)
library(stringr)
multicentricLabelPwd = './zhangWei2015_multicentricLabel30tcga.csv'

mgbmLabelPwd = './mgbmSgbmLabelMerge230907.csv'
multicentricDf = data.frame(fread(multicentricLabelPwd))
mgbmLabelDf = data.frame(fread(mgbmLabelPwd))
mgbmLabelDf = mgbmLabelDf[mgbmLabelDf$bigSmall != 'SGBM', ]
mgbmLabelDf$mriClass = multicentricDf$mri.classification[match(mgbmLabelDf$ID,
                                                               multicentricDf$patientId)]
colnames(mgbmLabelDf)
mgbmLabelDfClean = mgbmLabelDf[, c('ID', 'bigSmall', 'mriClass')]
mgbmLabelDfCleanMelt = reshape2::melt(mgbmLabelDfClean, id = 'ID')
colnames(mgbmLabelDfCleanMelt) = c('ID', 'labelName', 'label')
mgbmLabelDfCleanMelt$label[mgbmLabelDfCleanMelt$label == '<NA>'] = NA
mgbmLabelDfCleanMelt$ID = factor(mgbmLabelDfCleanMelt$ID, 
                                 levels = mgbmLabelDf$ID[order(mgbmLabelDf$bigSmall, mgbmLabelDf$mriClass)])

mgbmLabelDfCleanMelt$label =  str_to_title(mgbmLabelDfCleanMelt$label)

mgbmLabelDfCleanMelt$label = factor(mgbmLabelDfCleanMelt$label,
                                    levels = c('Big', 'Small', 'Multifocal', 'Multicentric'))

mgbmLabelDfCleanMelt$labelName = factor(mgbmLabelDfCleanMelt$labelName,
                                        levels = c('mriClass',
                                                   'bigSmall'
                                        ))

ggplot(mgbmLabelDfCleanMelt, mapping = aes(x = ID,
                                           y = labelName,
                                           fill = label)) +
  geom_tile(col = 'black', size = 0.3) +
  theme_classic2() +
  theme(text = element_text(size = 15),
        axis.ticks = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        
  ) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8", '#fb9a99', '#ff7f00'), 
                    na.value = NA) +
  labs(y = '',
       x = '',
       fill = '') 
table(mgbmLabelDf[, c('bigSmall', 'mriClass')])
fisher.test(table(mgbmLabelDf[, c('bigSmall', 'mriClass')]))

