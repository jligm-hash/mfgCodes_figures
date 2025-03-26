

library(ggplot2)
library(scales)
library(ggsci)
library(scico)
library(RColorBrewer)

library(reshape)
library(igraph)
library(rgl)
library(oce)
library(plot3D)
library(ggrepel)
dataFile = './bestSim230807.RData'
load(dataFile)

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
  cellDF = data.frame(originX,originY,originZ,ancester,gen)
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
myColor = rev(brewer.pal(9, 'Greys'))

myColor[1] = '#9e0142'
myColor[2] = '#d53e4f'

ggplot(growthCurveDf, mapping = aes(x=as.numeric(time)*1/0.69, 
                                    y=cellNumber*1e3, 
                                    group=lesion, color=lesion)) +
  theme_classic() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                expand = c(0, 0.1),
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  geom_hline(yintercept = 1e6, linetype = "dashed", 
             size = 1.5, col = '#0066B4') +
  geom_hline(yintercept = 1e9, linetype = "dashed", 
             size = 1.5, col = 'red') +
  labs(x = 'Generation',
       y = 'Cell Number',
       color = 'Lesion') +
  scale_colour_manual(values = myColor) +
  geom_line(size = 1, alpha = 0.4) +
  theme(text = element_text(size=20),
        plot.margin = margin(rep(20, 4)), axis.ticks.x = element_blank()) # set the margin plot
{
  
  n_geno = 1
  geno = NULL
  
  geno[1] = n_geno
  
  for (i in 2:length(parent)) {
    parentID = parent[i]
    
    cellMut = ca[i]
    parentMut = ca[parentID]
    
    if(cellMut == parentMut){
      geno[i] = geno[parentID] # the same with parent geno
    } else if ( cellMut > parentMut){
      n_geno = n_geno + 1
      geno[i] = n_geno
    }else{
      print('Error')
      break
    }
    
  }
  
  genoDf = data.frame(geno,ca)
  
  library(plyr)
  counts = ddply(genoDf, .(genoDf$geno, genoDf$ca), nrow)
  counts = counts[order(counts$V1,decreasing = T),]
  colnames(counts) = c('Geno','MutNum','Freq')
  head(counts)
  dim(counts)
  counts$col = 'grey'
  
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
  counts[counts$Freq>= counts$Freq[100],]$col = rainbow(100)
  counts$col[2:7] = sci7ColArrOrder 
  
  counts[counts$MutNum==1,]$col = 'grey90'
  
  colByGenoType = counts
  
  head(colByGenoType)
  }

{
  selectSpace = 1
  index = (spaceArr == selectSpace)
  radius = radiusArr[selectSpace]
  
  cellDF = data.frame(originX,originY,
                      originZ,ancester,gen,
                      spaceArr)
  colnames(cellDF) = c('x', 'y', 'z', 'genIndex', 'time',
                       'space') # add the space
  pltDF = cellDF
  oriDf = cellDF # finalDf
}

{
  lesionNumber = n_space
  
  cls = myColor
  
  colDf = data.frame(colInd = 1:length(cls), colVal = cls)
  colByGenIndex = colDf[match(pltDF$genIndex, colDf$colInd),]$colVal
}

myRepeat = function(shortVector, longVector){
  shortLen = length(shortVector)
  longLen = length(longVector)
  
  res = rep(shortVector, times = round(longLen/shortLen) + 1)
  return(res[1:longLen])
}

{
  myRainbow = myRepeat(rainbow(lesionNumber), originZ)
  maxZ = max(originZ)
  phyDf = data.frame(Lesion = 1:length(lesionMother),
                     Mother = lesionMother,
                     lesionId = 1:length(lesionMother))
  phyDfMelt = melt(phyDf, measure =c('Lesion', 'Mother'))
  phyDfMelt$variable = factor(phyDfMelt$variable,levels = c('Mother','Lesion'))
}

maxRadius = max(abs(pltDF[, c('x', 'y', 'z')]))
selectSpace = 1
{
  index = spaceArr == selectSpace
  colByGenoIndex = colByGenoType[match(geno, colByGenoType$Geno),]$col
  
}
{# visualize top 11 lesions
  xchangArr = 2*maxRadius*c(0,-1,-1,0,1,1,1,0,-1,0,0)
  zchangArr = 2*maxRadius*c(0,0,1,1,1,0,-1,-1,-1,0,0)
  ychangArr = 2*maxRadius*c(0,0,0,0,0,0,0,0,0,1,-1)
  
  maxRadius = max(abs(pltDF[, c('x', 'y', 'z')]))
  pltDFAll = pltDF[, c('x', 'y', 'z')]
  for (i in 1:11) {
    selectSpace = i
    index = spaceArr == selectSpace
    if(i>n_space){
      break
    }
    pltDFAll[index,]$x =  pltDF[index, ]$x + xchangArr[selectSpace]
    pltDFAll[index,]$y =  pltDF[index, ]$y + ychangArr[selectSpace]
    pltDFAll[index,]$z =  pltDF[index, ]$z + zchangArr[selectSpace]
  }
  
  colByGenoIndex = colByGenoType[match(geno, colByGenoType$Geno),]$col
  
  par3d(
    zoom = 0.5, 
    userMatrix = matrix(c(0.9418525,-0.2365354,0.2386722,0,0.104028,-0.4701345,-0.8764427,0,0.3195179,0.8503081,-0.4181913,0,0,0,0,1), ncol = 4), 
    windowRect=c(1920,  225, 2700,  995))
  clear3d()
  
  plot3d(
    pltDFAll,
    type = 'p',
    col = colByGenIndex, # color by lesion
    axes = F,
    xlab = "",
    ylab = "",
    zlab = "",
    xlim = 3*c(-maxRadius, maxRadius),
    ylim = 3*c(-maxRadius, maxRadius),
    zlim = 3*c(-maxRadius, maxRadius),
  )
}
phyDf = data.frame(Lesion = 1:length(lesionMother),
                   Mother = lesionMother,
                   lesionId = 1:length(lesionMother))

net = graph.data.frame(phyDf[-1,c('Mother', 'Lesion')],directed = T)
plot(net, 
     edge.arrow.size= 0.5,
     vertex.color=cls,
     vertex.size=20,
     vertex.label.color = 'black')

