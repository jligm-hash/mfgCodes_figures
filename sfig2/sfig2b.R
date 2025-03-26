
growthCurveAllDf = readRDS('./TCGA-27-1830.rds')
load('./TCGA-27-1830.rds.RData')

load('/Users/jiabao/Documents/0_projects/mfg/surgery4mostV134/bestSim4patient/hpc3new230912/TCGA-02-0033.RData')

saveRDS(growthCurveAllDf, file = '/Users/jiabao/Documents/1_github/mfgCodes/code4figs24112001/sfig2/TCGA-02-0033.rds')
load('./TCGA-02-0033.rds.RData')
if (length(unique(growthCurveAllDf$lesion)) > 100 ) {
  growthCurveAllDf2 = growthCurveAllDf[growthCurveAllDf$lesion %in% as.factor(1:100), ]
}else{
  growthCurveAllDf2 = growthCurveAllDf
}
library(ggplot2)
library(scales)
library(ggsci)
library(scico)
nColOfGgsci = 7
p = ggplot(data.frame(x = 1:nColOfGgsci,
                      y = 1:nColOfGgsci),
           mapping = aes(x = x,
                         y = y,
                         col =  as.factor(x))) +
  geom_point() +
  scale_color_startrek()
g = ggplot_build(p)
sci7ColArr = unique(g$data[[1]]["colour"])[[1]]
sci7ColArrOrder = sci7ColArr[c(6,5,4,3,2,1,7)]
cls =c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')
cls[1:6] = sci7ColArrOrder

ggplot(growthCurveAllDf2, mapping = aes(x=as.numeric(time)*1/0.69, 
                                        y=cellNumber*1e3, 
                                        group=lesion, color=lesion)) +
  theme_classic() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                expand = c(0, 0)
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  geom_vline(xintercept = surgeryTime*1/0.69, linetype = "dashed", col = 'red') +
  geom_vline(xintercept = recurrentTime*1/0.69, linetype = "dashed", col = 'blue') +
  geom_vline(xintercept = recurrentVisibleTime*1/0.69, linetype = "dashed", col = 'blue') +
  geom_hline(yintercept = setVisibleCutoff*1e3, linetype = "dashed", col = 'orange') +
  geom_hline(yintercept = setStopCellNum*1e3, linetype = "dashed", col = 'red') +
  labs(x = 'Days',
       y = 'Cell Number',
       color = 'Lesion') +
  scale_colour_manual(values = rainbow(length(unique(growthCurveAllDf2$lesion)))) +
  geom_line(linewidth = 0.5) +
  theme(text = element_text(size=20),
        plot.margin = margin(rep(20,4))) + # set the margin plot
  guides(col = 'none')

