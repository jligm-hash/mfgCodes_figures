
library(ggplot2)
library(scales)
library(ggsci)
library(scico)
library(reshape)
library(igraph)
library(ggrepel)
library(colorspace)
library(rgl)
library(oce)
library(plot3D)
library(RColorBrewer)
cellDF = readRDS('/Users/jiabao/Documents/1_github/mfgCodes/code4figs24112001/sfig2/cellDf.rds')

lesionMother = readRDS('./lesionMotherArr.rds')

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

{
  growthCurveTab = table(cellDF[,c('ancester', 'gen')])
  growthCurveDf = as.data.frame(t(growthCurveTab))
  colnames(growthCurveDf) = c('time', 'lesion', 'cellNumber')
  for (tmpLesion in unique(growthCurveDf$lesion)) {
    growthCurveDf$cellNumber[growthCurveDf$lesion == tmpLesion] =
      cumsum(growthCurveDf$cellNumber[growthCurveDf$lesion == tmpLesion])
  }
  growthCurveDf = growthCurveDf[growthCurveDf$cellNumber>0,]
}

growthCurveDf[nrow(growthCurveDf) + 1,] = c(1,1,1)

growthCurveDf$lesion = paste0('Lesion ',growthCurveDf$lesion)
growthCurveDf$lesion = factor(growthCurveDf$lesion, levels = paste0('Lesion ', 1:10))
plt4growthCurve = ggplot(growthCurveDf, mapping = aes(x=as.numeric(time)*1/0.69, 
                                                      y=cellNumber*1e3, 
                                                      group=lesion, color=lesion)) +
  theme_classic() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                expand = c(0, 0)
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  geom_hline(yintercept = 1e6, linetype = "dashed",col = 'blue') +
  geom_hline(yintercept = 1e9, linetype = "dashed",col = 'red') +
  labs(x = 'Days',
       y = 'Cell Number',
       color = 'Lesion') +
  scale_colour_manual(values = cls) +
  geom_line(linewidth = 0.5) +
  theme(text = element_text(size=20),
        plot.margin = margin(rep(20,4))) + # set the margin plot 
  guides(col = 'none')
print(plt4growthCurve)
phyDf = data.frame(Lesion = 1:length(lesionMother),
                   Mother = lesionMother,
                   lesionId = 1:length(lesionMother))

net = graph.data.frame(phyDf[-1,c('Mother', 'Lesion')],directed = T)
plot(net, 
     edge.arrow.size= 0.5,
     vertex.color=cls,
     vertex.size=20,
     vertex.label.color = 'black')
finalCellNumDf = growthCurveDf[as.numeric(growthCurveDf$time) == max(as.numeric(growthCurveDf$time)),]
cutOff = 1e3
finalCellNumDf = finalCellNumDf[finalCellNumDf$cellNumber> cutOff,]
finalCellNumDf$prop = prop.table(finalCellNumDf$cellNumber)

plt4piePlot = ggplot(finalCellNumDf, aes(x="", y=prop, fill=lesion))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cls[1:dim(finalCellNumDf)[1]]) +
  theme_void() +
  guides(fill = 'none')

print(plt4piePlot)

