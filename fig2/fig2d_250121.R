

load(file = './readSimV1345_240702.RData')
dim(resSimAll)

colnames(resSimAll)

resSimAll$migr = resSimAll$migr/1e3 # in the nich level

library(ggplot2)
dim(resSimAll)
colnames(resSimAll)

mid<-mean(resSimAll$migr[resSimAll$firstMigrNum > 0])
lm(firstMigrNum~migr, data = resSimAll[resSimAll$firstMigrNum > 0, ])

cor.test(resSimAll[resSimAll$firstMigrNum > 0, ]$firstMigrNum,
         resSimAll[resSimAll$firstMigrNum > 0, ]$migr)$p.value

cor.test(resSimAll[resSimAll$firstMigrNum > 0, ]$firstMigrGen,
         resSimAll[resSimAll$firstMigrNum > 0, ]$migr)$p.value

ggplot(resSimAll[resSimAll$firstMigrNum > 0, ], 
       mapping = aes(x = migr,
                     y = firstMigrGen
       )) +
  geom_boxplot(mapping = aes(col = migr,
                             group = migr),
               outlier.shape = NA) +
  geom_jitter(width = 0.02, mapping = aes(col = migr,
                                          group = migr),
              alpha = 0.2, size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_gradient(low="#F5E7B2", 
                       high="#A91D3A") +
  
  geom_smooth(method = 'lm', col = '#C68484') +
  theme_classic() +
  theme(text = element_text(size = 15)) +
  guides(color = 'none') +
  labs(x = 'Distant migration rate',
       y = 'Time/Generation')
resMeanValueArr = c()
resMigrValueArr = sort(unique(resSimAll[resSimAll$firstMigrNum > 0, ]$migr))

for (tmpMigr in resMigrValueArr) {
  tmpMean = median(resSimAll$firstMigrGen[resSimAll$firstMigrNum > 0 & resSimAll$migr == tmpMigr])
  resMeanValueArr = c(resMeanValueArr, tmpMean)
}

resMeanMigrGenDf = data.frame(migr = resMigrValueArr,
                              meanMigrGen = resMeanValueArr)

cor.test(resMeanMigrGenDf$migr,
         resMeanMigrGenDf$meanMigrGen)

