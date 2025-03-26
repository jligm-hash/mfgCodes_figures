

library(ggbeeswarm)
library(ggplot2)
library(rstatix)

load('./resSimAllMigr0.RData')
load('./readSimV134.RData')

colnames(resSimAll)
table(resSimAll$migr)

resSimAll = resSimAll[resSimAll$migr %in% c(1e-6, 5e-05),] # , 1e-4 2.5e-05,

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

colorList = c("#C80036", "#FF6969", "#FFBF78")
setWidth4boxplot = 0.5
setAlpha4beeswarm = 0.3
pltList = c()

setWidth4boxplot = 0.5
setAlpha4beeswarm = 0.7

pltList[[1]] = ggplot(resSimAllMeltDf[resSimAllMeltDf$lesionType == 'lesionNum',], mapping = aes(x = migrGroup,
                                                                                                 y = Number,
                                                                                                 col = migrGroup,
                                                                                                 
)) +
  geom_boxplot(width = setWidth4boxplot, position="dodge", fill = 'white', col=colorList[1],
               outlier.shape = NA) + 
  geom_beeswarm(alpha = setAlpha4beeswarm, col=colorList[1], size = 0.5) +
  labs(
    x = '',
    y = 'Number of simulated lesions',
  ) +
  theme_classic() +
  theme(text = element_text(size=20),
        plot.title = element_text(hjust = 0.5),
  ) +
  guides(col='none') 
pltList[[2]] = ggplot(resSimAllMeltDf[resSimAllMeltDf$lesionType == 'num',], mapping = aes(x = migrGroup,
                                                                                           y = Number,
                                                                                           col = migrGroup,
                                                                                           
)) +
  geom_boxplot(width = setWidth4boxplot, position="dodge", fill = 'white', col=colorList[2],
               outlier.shape = NA) + 
  geom_beeswarm(alpha = setAlpha4beeswarm, col=colorList[2], size = 0.5) +
  labs(
    y = '',
    x = 'Distant migration rate',
  ) +
  theme_classic() +
  theme(text = element_text(size=20),
        plot.title = element_text(hjust = 0.5),
  ) +
  guides(col='none') 
pltList[[3]] = ggplot(resSimAllMeltDf[resSimAllMeltDf$lesionType == 'invisibleLesionNum',], mapping = aes(x = migrGroup,
                                                                                                          y = Number,
                                                                                                          col = migrGroup,
                                                                                                          
)) +
  geom_boxplot(width = setWidth4boxplot, position="dodge", fill = 'white', col=colorList[3],
               outlier.shape = NA) + # width = 0.5,
  geom_beeswarm(alpha = setAlpha4beeswarm, col=colorList[3], size = 0.5) +
  
  theme_classic() +
  theme(text = element_text(size=20),
        plot.title = element_text(hjust = 0.5)) +
  guides(col='none') +
  labs(x = "",
       y = "") 

ggarrange(plotlist = pltList, nrow = 1, ncol = 3)
pairwise_wilcox_test(resSimAllMeltDf[resSimAllMeltDf$lesionType == 'lesionNum',], 
                     Number ~ migrGroup, p.adjust.method = "bonferroni")

pairwise_wilcox_test(resSimAllMeltDf[resSimAllMeltDf$lesionType == 'num',], 
                     Number ~ migrGroup, p.adjust.method = "bonferroni")

pairwise_wilcox_test(resSimAllMeltDf[resSimAllMeltDf$lesionType == 'invisibleLesionNum',], 
                     Number ~ migrGroup, p.adjust.method = "bonferroni")

