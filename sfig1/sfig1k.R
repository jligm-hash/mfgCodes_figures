
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
bigSmallLabelPwd = './mgbmSgbmLabelMerge230907.csv'
bigSmallLabelDf = data.frame(fread(bigSmallLabelPwd))
bigSmallLabelDf4mgbm = bigSmallLabelDf[bigSmallLabelDf$bigSmall != 'SGBM', ]

shapeFeatureAllPwd = './mgbmLesionFeature.csv'
shapeFeatureDf = data.frame(fread(shapeFeatureAllPwd))
shapeFeatureDf$bigSmall = bigSmallLabelDf4mgbm$bigSmall[match(shapeFeatureDf$patient, bigSmallLabelDf4mgbm$ID)]
ggplot(shapeFeatureDf, mapping = aes(x = bigSmall,
                                     y = original_shape_SurfaceVolumeRatio_T2,
                                     col = bigSmall)) +
  geom_boxplot(width = 0.3) +
  geom_beeswarm(alpha = 0.8, size = 0.5) +
  theme_classic() +
  theme(text = element_text(size = 15)) +
  scale_color_brewer(palette = 'Set1') +
  guides(col = 'none') +
  stat_compare_means(comparisons = list(c(1, 2)),
                     method = 't.test') +
  labs(x = 'M-GBM group',
       y = 'Surface volume ratio of L2')

