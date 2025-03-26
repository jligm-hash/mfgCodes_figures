
library(data.table)
library(progress)
library(abc)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(ggsci)
library(scico)
pltList = c()
pltFigureCount = 0

{
  
  setTolerance = 0.2
  
  
  bigSmallLabelPwd = '/Users/jiabao/Documents/1_github/mfgCodes/expression/mgbmSgbmLabel.csv'
  bigSmallDf = data.frame(fread(bigSmallLabelPwd))
  bigSmallDf = bigSmallDf[bigSmallDf$bigSmall != 'SGBM',]
  
  
  {
    
    myListToCol = function(x){
      ColDf = data.frame(col = table(x))
      colnames(ColDf) = c('Value', 'Freq')
      ColDf$Col = rainbow(dim(ColDf)[1])
      return(ColDf$Col[match(x, ColDf$Value)])
    }
    
    myListToLog10 = function(x){
      return(ifelse(x==0, 0, log10(x)))
    }
    
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
    
    pwd = '/Users/jiabao/Documents/1_github/mfgCodes/plotAbcHeatmap'
    
    load(file = file.path(pwd, 'readSimV134.RData'))
    dim(resSimAll)
    
    colnames(resSimAll)
    
    resSimAll$migr = resSimAll$migr/1e3 # in the nich level
    
    pwd = '/Users/jiabao/Documents/1_github/mfgCodes/plotAbcHeatmap'
    realSizeDf = data.frame(fread(file.path(pwd, 'realTCGA211017.csv')))
    dim(realSizeDf)
    realSizeDf[is.na(realSizeDf)] = 'Unknown'
    dim(realSizeDf)
    
    resSimAll
    
    obsDat = realSizeDf[,c('Num', paste0('T',1:6))]
    patientName = realSizeDf$ID
    obsDat[1,]
    
    simDat = resSimAll[resSimAll$num>=1,]
    simDat[is.na(simDat)] = 0
    
    
    myMad = function(x){
      return(median(1/qnorm(3/4)*abs(x-median(x))))
    }
    
    myAbcPr = function(setTarget, setSimulation, setLabel, setTol){
      target = setTarget
      sim = setSimulation
      modelLabel1 = as.factor(setLabel)
      tol = setTol
      resMad = c()
      for (i in 1:ncol(sim)) {
        madF = myMad(as.numeric(sim[,i]))
        madF = ifelse(madF==0, 1, madF)
        resMad = c(resMad, madF)
      }
      
      target.scale = target
      sim.scale = sim
      for (i in 1:ncol(sim)) {
        target.scale[,i] = target[,i] / resMad[i]
        sim.scale[,i] = sim[,i] / resMad[i]
      }
      
      resDist = 0
      for (i in 1:ncol(sim.scale)) {
        resDist.tmp = (sim.scale[,i] - target.scale[,i]) ^2
        resDist = resDist + resDist.tmp
      }
      ntol = ceiling(nrow(sim.scale) * tol)
      cutOffN = sort(resDist)[ntol]
      postRindex = resDist <= cutOffN
      resTmpDf = aggregate(rank(-resDist[postRindex]), 
                           by=list(label = modelLabel1[postRindex]), FUN=sum)
      resTmpDf = resTmpDf[order(resTmpDf$label),]
      
      tabTmp = table(modelLabel1[postRindex])
      
      for (i in 1:dim(resTmpDf)[1]) {
        tabTmp[names(tabTmp) == resTmpDf$label[i]] = resTmpDf$x[i]
      }
      return(prop.table(tabTmp))
    }
    
    myAbcPrBestSim = function(setTarget, setSimulation, setLabel, setTol){
      target = setTarget
      sim = setSimulation
      modelLabel1 = as.factor(setLabel)
      tol = setTol
      resMad = c()
      for (i in 1:ncol(sim)) {
        madF = myMad(as.numeric(sim[,i]))
        madF = ifelse(madF==0, 1, madF)
        resMad = c(resMad, madF)
      }
      
      target.scale = target
      sim.scale = sim
      for (i in 1:ncol(sim)) {
        target.scale[,i] = target[,i] / resMad[i]
        sim.scale[,i] = sim[,i] / resMad[i]
      }
      
      resDist = 0
      for (i in 1:ncol(sim.scale)) {
        resDist.tmp = (sim.scale[,i] - target.scale[,i]) ^2
        resDist = resDist + resDist.tmp
      }
      ntol = ceiling(nrow(sim.scale) * tol)
      cutOffN = sort(resDist)[ntol]
      postRindex = resDist <= cutOffN
      return(postRindex)
    }
    
    
    myResMad = function(setTarget, setSimulation){
      target = setTarget
      sim = setSimulation
      resMad = c()
      for (i in 1:ncol(sim)) {
        madF = myMad(as.numeric(sim[,i]))
        madF = ifelse(madF==0, 1, madF)
        resMad = c(resMad, madF)
      }
      return(resMad)
    }
    
    simDatOnly = simDat
    simDatOnly[, paste0('P',1:10)] = 100 * simDatOnly[, paste0('P',1:10)]
    
    unique(simDatOnly$migr)
    unique(simDatOnly$select)
    unique(simDatOnly$distMigrDrvNum)
    
    
    pb = progress_bar$new(total = length(patientName), clear = FALSE)
    rm(resPostAllBest)
    
    for (i in 1:length(patientName)) {
      pb$tick()
      resPost = myAbcPrBestSim(setTarget = obsDat[i,1:7], 
                               setSimulation = simDatOnly[,c('num', paste0('P',1:6))], 
                               setLabel = simDatOnly$migr, 
                               setTol =setTolerance) 
      
      tmpPostRes = simDatOnly[resPost,]
      tmpPostRes$patientId = patientName[i]
      
      
      if(i == 1){
        resPostAllBest = tmpPostRes
      }else{
        resPostAllBest = rbind(resPostAllBest, tmpPostRes)
      }
    }
    
    resPostAllBest
    
    
    resPostAllBest$group = bigSmallDf$bigSmall[match(resPostAllBest$patientId,
                                                     bigSmallDf$ID)]
  }
  
  
  library(ggpubr)
  my_comparison = list(c(1, 2))
  
  colnames(resPostAllBest)
  
  
  resMean = c()
  for (i in 1:length(patientName)) {
    tmpRes = mean(resPostAllBest$firstMigrNum[resPostAllBest$patientId == patientName[i]])
    resMean = c(resMean, tmpRes)
  }
  
  resMeanDf4bigSmall = data.frame(firstMigrNumMean = resMean,
                                  patientId = patientName)
  resMeanDf4bigSmall$group = bigSmallDf$bigSmall[match(resMeanDf4bigSmall$patientId,
                                                       bigSmallDf$ID)]
  
  
  library(ggbeeswarm)
  
  pltFigureCount = pltFigureCount + 1
  
  ggplot(resMeanDf4bigSmall, mapping = aes(x = group,
                                           y = firstMigrNumMean,
                                           col = group)) +
    geom_boxplot(width = 0.3, outlier.shape = NA) +
    geom_beeswarm(cex = 3, alpha = 0.8, size = 0.5) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
    ) +
    stat_compare_means(comparisons = my_comparison) +
    scale_color_brewer(palette = 'Set1') +
    guides(col = 'none') +
    labs(title = paste0('Tolerance = ', setTolerance)) +
    scale_y_continuous(limits = c(0, 5e5))
  
  resMean4migr = c()
  for (i in 1:length(patientName)) {
    tmpRes = mean(resPostAllBest$migr[resPostAllBest$patientId == patientName[i]])
    resMean4migr = c(resMean4migr, tmpRes)
  }
  
  resMeanDf4bigSmall4migr = data.frame(migrMean = resMean4migr,
                                       patientId = patientName)
  
  
  
  
  
  
  resMeanDf4bigSmall4migr$group = bigSmallDf$bigSmall[match(resMeanDf4bigSmall4migr$patientId,
                                                            bigSmallDf$ID)]
  
  pltFigureCount = pltFigureCount + 1
  
  
  ggplot(resMeanDf4bigSmall4migr, mapping = aes(x = group,
                                                y = migrMean,
                                                col = group)) +
    geom_boxplot(width = 0.3, outlier.shape = NA) +
    geom_beeswarm(cex = 3, alpha = 0.8, size = 0.5) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
    )+
    stat_compare_means(comparisons = my_comparison) +
    scale_color_brewer(palette = 'Set1') +
    guides(col = 'none') +
    labs(title = paste0('Tolerance = ', setTolerance)) +
    scale_y_continuous(limits = c(0, 7e-8))
  
}

