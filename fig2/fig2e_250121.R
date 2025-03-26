
library(ggbeeswarm)
library(ggplot2)
library(rstatix)
load('./resSimAllMigr0.RData')
load('./readSimV134.RData')

colnames(resSimAll)
table(resSimAll$migr)

resSimAll = resSimAll[resSimAll$migr %in% c(1e-6,  5e-05),] # , 1e-4 2.5e-05,

summary(resSimAll$firstMigrNum)

resSimAll = rbind(resSimAll,
                  resSimAll4migr0) 
selectNum4simulation = 20

resSimAll = resSimAll[resSimAll$seed %in% 1:selectNum4simulation, ]
resSimAll$group = as.factor(resSimAll$migr)

resSimAll$invisibleLesionNum = resSimAll$lesionNum - resSimAll$num

table(resSimAll$group)
library(RColorBrewer)
library(ggpubr)

colnames(resSimAll)

resSimAllMeltDf = reshape2::melt(resSimAll[,c('group', 'lesionNum', 'num', 'invisibleLesionNum')])

colnames(resSimAllMeltDf) = c('migrGroup', 'lesionType', 'Number')
library(ggsignif)

colnames(resSimAllMeltDf)

resSimAllMeltDf$grp2 = " "
table(resSimAllMeltDf$lesionType)

colorList = brewer.pal(5, "Set2")
setWidth4boxplot = 0.6
setAlpha4beeswarm = 0.3
colnames(resSimAllMeltDf)

pairwise_wilcox_test(resSimAllMeltDf[resSimAllMeltDf$lesionType == 'lesionNum', ], Number ~ migrGroup, p.adjust.method = "bonferroni")
pairwise_wilcox_test(resSimAllMeltDf[resSimAllMeltDf$lesionType == 'num', ], Number ~ migrGroup, p.adjust.method = "bonferroni")
pairwise_wilcox_test(resSimAllMeltDf[resSimAllMeltDf$lesionType == 'invisibleLesionNum', ], Number ~ migrGroup, p.adjust.method = "bonferroni")
resSimAll$l2l1ratio = resSimAll$T2/resSimAll$T1
resSimAllMeltDf4l1l2 = reshape2::melt(resSimAll[,c('group', 'P1', 'P2',
                                                   'l2l1ratio', 'T1', 'T2')])

colnames(resSimAllMeltDf4l1l2) = c('migrGroup', 'lesion', 'relativeSize')
colorList2 = brewer.pal(5, name = 'Set1')

table(resSimAllMeltDf4l1l2$lesion)
pltList2 = c()
colorList2 = c('#E9C46A', '#E76F51')

pltList2[[1]] =
  ggplot(resSimAllMeltDf4l1l2[resSimAllMeltDf4l1l2$lesion == "P2",], 
         mapping = aes(x = migrGroup, y = relativeSize, col = lesion,
                       fill = lesion
         )) +
  geom_boxplot(col = colorList2[2], position="dodge", fill='white', outlier.shape = NA, width=setWidth4boxplot) +
  geom_beeswarm(col = colorList2[2], alpha = setAlpha4beeswarm, size = 0.5) +
  scale_y_continuous(limits = c(0, 0.5)) +
  
  labs(x = '', y = ''
  ) +
  theme_classic() +
  theme(text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank()
  ) +
  guides(col = 'none', fill = 'none')
pltList2[[3]] =
  ggplot(resSimAllMeltDf4l1l2[resSimAllMeltDf4l1l2$lesion == "P1",],
         mapping = aes(x = migrGroup, y = relativeSize, col = lesion,
                       fill = lesion
         )) +
  geom_boxplot(col = colorList2[1], position="dodge", fill='white', outlier.shape = NA, width=setWidth4boxplot) +
  geom_beeswarm(col = colorList2[1], alpha = setAlpha4beeswarm, size = 0.5) +
  scale_y_continuous(limits = c(0.5, 1)) +
  labs(x = 'Distant migration rate', y = ''
  ) +
  theme_classic() +
  theme(text = element_text(size=15),
        plot.title = element_text(hjust = 0.5)) +
  guides(col = 'none', fill = 'none')
pltList2[[2]] =
  ggplot(resSimAllMeltDf4l1l2[resSimAllMeltDf4l1l2$lesion == "T2",],
         mapping = aes(x = migrGroup, y = relativeSize, col = lesion,
                       fill = lesion
         )) +
  geom_boxplot(col = colorList2[2], position="dodge", fill='white', outlier.shape = NA, width=setWidth4boxplot) +
  geom_beeswarm(col = colorList2[2], alpha = setAlpha4beeswarm, size = 0.5) +
  scale_y_continuous(limits = c(0, 5e5)) +
  labs(x = '', y = ''
  ) +
  theme_classic() +
  theme(text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8)
  ) +
  guides(col = 'none', fill = 'none')
pltList2[[4]] =
  ggplot(resSimAllMeltDf4l1l2[resSimAllMeltDf4l1l2$lesion == "T1",],
         mapping = aes(x = migrGroup, y = relativeSize, col = lesion,
                       fill = lesion
         )) +
  geom_boxplot(col = colorList2[1], position="dodge", fill='white', outlier.shape = NA, width=setWidth4boxplot) +
  geom_beeswarm(col = colorList2[1], alpha = setAlpha4beeswarm, size = 0.5) +
  labs(x = 'Distant migration rate', y = ''
  ) +
  scale_y_continuous(limits = c(5e5, 1e6)) +
  theme_classic() +
  theme(text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 8)
  ) +
  guides(col = 'none', fill = 'none')
ggarrange(plotlist = pltList2, ncol = 2, nrow = 2)
colnames(resSimAllMeltDf4l1l2)
pairwise_wilcox_test(resSimAllMeltDf4l1l2[resSimAllMeltDf4l1l2$lesion == 'P1', ], relativeSize ~ migrGroup, p.adjust.method = "bonferroni")
pairwise_wilcox_test(resSimAllMeltDf4l1l2[resSimAllMeltDf4l1l2$lesion == 'P2', ], relativeSize ~ migrGroup, p.adjust.method = "bonferroni")
pairwise_wilcox_test(resSimAllMeltDf4l1l2[resSimAllMeltDf4l1l2$lesion == 'T1', ], relativeSize ~ migrGroup, p.adjust.method = "bonferroni")
pairwise_wilcox_test(resSimAllMeltDf4l1l2[resSimAllMeltDf4l1l2$lesion == 'T2', ], relativeSize ~ migrGroup, p.adjust.method = "bonferroni")

