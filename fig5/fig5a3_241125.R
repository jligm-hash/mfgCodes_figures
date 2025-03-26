
library(data.table)
library(progress)
library(abc)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(ggsci)
library(scico)

library(scales)
library(reshape)
library(igraph)
library(ggrepel)
library(colorspace)
library(RColorBrewer)
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
load(file = '/Users/jiabao/Documents/1_github/mfgCodes/abc4eachPatient/readSimV134.RData')
dim(resSimAll)

colnames(resSimAll)

resSimAll$migr = resSimAll$migr/1e3 # in the nich level

{
  
  realSizeDf = data.frame(fread('/Users/jiabao/Documents/1_github/mfgCodes/abc4eachPatient/realTCGA211017.csv'))
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
  rm(resPostAllOnly)
  for (i in 1:length(patientName)) {
    pb$tick()
    resPost = myAbcPr(setTarget = obsDat[i,1:7], 
                      setSimulation = simDatOnly[,c('num', paste0('P',1:6))], 
                      setLabel = simDatOnly$migr, 
                      setTol =0.1) 
    sumPost = data.frame(resPost)
    tmpPostRes = data.frame(t(sumPost$Freq),
                            row.names = patientName[i]) 
    colnames(tmpPostRes) = sumPost$Var1
    if(i == 1){
      resPostAllOnly = tmpPostRes
    }else{
      resPostAllOnly = rbind(resPostAllOnly, tmpPostRes)
    }
  }
  
  resDisMigr = resPostAllOnly
  medianList = apply(resDisMigr, 1, max)
  
  maxNameList = c()
  for (i in 1:dim(resDisMigr)[1]) {
    tmpNameValue = names(resDisMigr)[resDisMigr[i,]==medianList[i]][1] # use the value of the 
    maxNameList = c(maxNameList, tmpNameValue)
  }
  
  patientNameListAfterOrder = rownames(resDisMigr)[order(as.numeric(maxNameList),medianList)]
  
  resDisMigrRm = resDisMigr
  migrOrderList2 = c()
  colNameList2 = colnames(resDisMigrRm)
  for (i in 1:length(colNameList2)) {
    tmpMigr = unlist(strsplit(colNameList2[i], '.T'))[1]
    migrOrderList2 = c(migrOrderList2, as.numeric(tmpMigr))
  }
  
  my_col_fun = colorRamp2(c(0, 1), c("white", 'red'))  # sci7ColArrOrder[1]
  
  cn = colnames(resDisMigr)
  rsBarCol = cbind(obsDat[,paste0('T',1:6)])
  ha = rowAnnotation(RS = anno_barplot(rsBarCol,
                                       gp = gpar(fill = sci7ColArrOrder, #c('#8dd3c7','#bebada','#80b1d3','#fdb462','#b3de69','#ffffb3'),
                                                 col = NA))) 
  
  resDisMigrRm = resDisMigr
  Heatmap(resDisMigrRm,
          name = "Probability",
          col = my_col_fun,
          cluster_columns = FALSE,
          column_order = order(migrOrderList2),
          
          cluster_rows = FALSE,
          row_order = rev(as.numeric(factor(levels = rownames(resDisMigr), patientNameListAfterOrder))),
          show_column_dend = F,
          show_row_dend = F,
          rect_gp = gpar(col = "grey90", lwd = 0.5),
          right_annotation = ha,
          column_title = "Migration rate",
          row_names_gp = gpar(fontsize = 7),
          show_column_names = FALSE,
          cell_fun = function(j, i, x, y, width, height, fill) {
            if(resDisMigrRm[i, j] == medianList[[i]])
              grid.points(x,y, pch = 21, size = unit(2, "mm"))
          },
          bottom_annotation = HeatmapAnnotation(text = anno_text(cn, rot = 85,just = "right"))
  )
  
  Heatmap(resDisMigrRm,
          name = "Probability",
          col = my_col_fun,
          cluster_columns = FALSE,
          column_order = order(migrOrderList2),
          cluster_rows = FALSE,
          row_order = rev(as.numeric(factor(levels = rownames(resDisMigr), patientNameListAfterOrder))),
          show_column_dend = F,
          show_row_dend = F,
          rect_gp = gpar(col = "grey90", lwd = 0.5),
          right_annotation = ha,
          column_title = "Migration rate",
          row_names_gp = gpar(fontsize = 7),
          show_column_names = FALSE,
          cell_fun = function(j, i, x, y, width, height, fill) {
            if(resDisMigrRm[i, j] == medianList[[i]])
              grid.points(x,y, pch = 21, size = unit(2, "mm"),
              )
          },
          bottom_annotation = HeatmapAnnotation(text = anno_text(cn, rot = 85,just = "right"))
  )
  
}
rowSplitOrder = data.frame(now = rev(as.numeric(factor(levels = rownames(resDisMigr), patientNameListAfterOrder))),
                           earlyLate = c(rep('High', 22), rep('Low',13)))
Heatmap(resDisMigrRm,
        name = "Probability",
        col = my_col_fun,
        cluster_columns = FALSE,
        column_order = order(migrOrderList2),
        cluster_rows = FALSE,
        row_order = rev(as.numeric(factor(levels = rownames(resDisMigr), patientNameListAfterOrder))),
        row_split = rowSplitOrder$earlyLate[order(rowSplitOrder$now)],
        column_split = factor(c(rep('Low',4), rep('High', 4)),
                              levels = c('Low', 'High')),
        show_column_dend = F,
        show_row_dend = F,
        rect_gp = gpar(col = 'grey',
                       lwd = 0.5),
        right_annotation = ha,
        column_title = "Migration rate",
        row_names_gp = gpar(fontsize = 7),
        show_column_names = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(resDisMigrRm[i, j] == medianList[[i]])
            grid.points(x,y, pch = 21, size = unit(2, "mm"))
        },
        bottom_annotation = HeatmapAnnotation(text = anno_text(cn, rot = 85,just = "right"))
)
  {
  
  i = 28
  tmpInferRes = resDisMigrRm[i,]
  tmpInferDf = data.frame(migrRate = as.numeric(names(tmpInferRes)),
                          prob = unlist(tmpInferRes, use.names = F))
  
  
  
  ggplot(tmpInferDf, mapping = aes(x = migrRate,
                                   y = prob)) +
    
    geom_line() +
    geom_point() +
    geom_ribbon(ymin=-Inf, aes(ymax=prob), fill='red', alpha=0.2) +
    geom_vline(xintercept = tmpInferDf$migrRate[tmpInferDf$prob == max(tmpInferDf$prob)], col = 'red', linetype = 'dashed') +
    
    scale_x_continuous(trans = 'log10') +
    scale_y_continuous(limits = c(0, max(tmpInferDf$prob) * 1.4), expand = c(0, 0)) +
    
    
    theme_classic() +
    labs(x = 'Migration rate (M)',
         y = 'Inferred probability',
         title = rownames(tmpInferRes)) +
    theme(text = element_text(size = 15),
          plot.title = element_text(hjust = 0.5))
  
  
  
}

