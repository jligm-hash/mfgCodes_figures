
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
rdataPwd = '/Users/jiabao/Documents/0_projects/mfg/bestSimulation240112/rdata'
rdataFileList = list.files('/Users/jiabao/Documents/0_projects/mfg/bestSimulation240112/rdata')

mgbmPatientDf = read.csv('/Users/jiabao/Documents/1_github/mfgCodes/previous/mriLabelAnalysis/mriLabel.csv')
mgbmPatientList = mgbmPatientDf$patientId[mgbmPatientDf$label == 'MGBM']
resPwd = '/Users/jiabao/Documents/0_projects/mfg/bestSimulation240112/newRes/growthCurve'
  {
  
  tmpPatientName = "TCGA-19-5954"
  
  
  print(tmpPatientName)
  tmpDataFileName = rdataFileList[grep(tmpPatientName, rdataFileList)]
  
  dataFile = file.path(rdataPwd, tmpDataFileName)
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

  cls = c("#67001f", "#980043", "#ce1256", "#e7298a", "#df65b0", "#c994c7", "#d4b9da",  "#525252", "#737373", "#969696", "#BDBDBD", "#D9D9D9", "#F0F0F0")
  
  
  cls[sort(lesions$Var1[lesions$Freq < 1e3])] = c("#525252", "#737373", "#969696", "#BDBDBD", "#D9D9D9", "#F0F0F0")[1:sum(lesions$Freq < 1e3)]
  
  {
    cellDF = data.frame(originX,originY,originZ,ancester,gen)
    growthCurveTab = table(cellDF[,c('ancester', 'gen')])
    
    growthCurveDf = as.data.frame(t(growthCurveTab))
    colnames(growthCurveDf) = c('time', 'lesion', 'cellNumber')
    
    
    growthCurveDf = growthCurveDf[as.numeric(as.character(growthCurveDf$lesion)) <= 15, ]
    
    for (tmpLesion in unique(growthCurveDf$lesion)) {
      growthCurveDf$cellNumber[growthCurveDf$lesion == tmpLesion] =
        cumsum(growthCurveDf$cellNumber[growthCurveDf$lesion == tmpLesion])
    }
    growthCurveDf = growthCurveDf[growthCurveDf$cellNumber>0,]
  }
  
  
  if(length(unique(growthCurveDf$lesion)) > 13){
    
    cls[14:length(unique(growthCurveDf$lesion))] = "#F0F0F0"
    
    
  }
  
  
  growthCurveDf[nrow(growthCurveDf) + 1,] = c(1,1,1)
  
  growthCurveDf$lesion = paste0('Lesion ',growthCurveDf$lesion)
  growthCurveDf$lesion = factor(growthCurveDf$lesion, levels = paste0('Lesion ', 1:length(unique(growthCurveDf$lesion))))
  
  xMin4growthCurve = (max(gen) + 0)/0.69 # 228
  xMax4growthCurve = (max(gen) + 0)/0.69 + 12 #240
  ggplot(growthCurveDf, mapping = aes(x=as.numeric(time)*1/0.69, 
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
    labs(x = 'Generation',
         y = 'Cell Number',
         color = 'Lesion') +
    scale_colour_manual(values = cls) +
    geom_line(linewidth = 1, alpha = 0.8) +
    
    annotate("rect", xmin = xMin4growthCurve, xmax = xMax4growthCurve, ymin = 0, ymax = 1e6,
             alpha = .3,
             fill = "grey") +
    annotate("rect", xmin = xMin4growthCurve, xmax = xMax4growthCurve, ymin = 1e6, ymax = 1e9,
             alpha = .3,fill = "red") +
    
    
    theme(text = element_text(size=20),
          plot.margin = margin(rep(20,4))) # set the margin plot
  
  
  
  
  phyDf = data.frame(Lesion = 1:length(lesionMother),
                     Mother = lesionMother,
                     lesionId = 1:length(lesionMother))
  
  net = graph.data.frame(phyDf[-1,c('Mother', 'Lesion')],directed = T)
  
  fullnames=V(net)$name
  fullnames[1:3]
  
  V(net)$familyname=fullnames
  V(net)$fullname=fullnames
  V(net)$name=fullnames # first name
  
  vcol=cls[1:length(fullnames)]
  
  V(net)$color=vcol
  
  
  
  plot(net,
       edge.arrow.size= 0.5,
       vertex.size=20,
       vertex.label.color = 'black')
  
  
  
  finalCellNumDf = growthCurveDf[as.numeric(growthCurveDf$time) == max(as.numeric(growthCurveDf$time)),]
  cutOff = 1e3
  finalCellNumDf = finalCellNumDf[finalCellNumDf$cellNumber> cutOff,]
  finalCellNumDf$prop = prop.table(finalCellNumDf$cellNumber)
  
  ggplot(finalCellNumDf, aes(x="", y=prop, fill=lesion))+
    geom_bar(width = 1, stat = "identity", alpha = 1) +
    coord_polar("y", start=0) +
    scale_fill_manual(values = cls[1:dim(finalCellNumDf)[1]]) +
    theme_void() +
    guides(fill = 'none')
  
  
  
  
  
  
}
