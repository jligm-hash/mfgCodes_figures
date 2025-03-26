

library(data.table)
library(progress)
library(abc)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(ggsci)
library(scico)

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
load(file = './readSimV134.RData')
dim(resSimAll)

colnames(resSimAll)

resSimAll$migr = resSimAll$migr/1e3 # in the nich level

{
  
  realSizeDf = data.frame(fread('./realTCGA211017.csv'))
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
  
  my_col_fun = colorRamp2(c(0, 1), c("white", 'red'))  
  
  cn = colnames(resDisMigr)
  rsBarCol = cbind(obsDat[,paste0('T',1:6)])
  ha = rowAnnotation(RS = anno_barplot(rsBarCol,
                                       gp = gpar(fill = sci7ColArrOrder, 
                                                 col = NA))) 
  
  
  
  library(data.table)
  bigSmallPwd = './mgbmSgbmLabel.csv'
  bigSmallDf = data.frame(fread(bigSmallPwd))
  
  patientName
  
  
  bigSmallLabelArr = bigSmallDf$bigSmall[match(patientName, bigSmallDf$ID)]
  bigSmallLabelDf = data.frame(patientId = patientName,
                               bigSmall = bigSmallLabelArr,
                               colCode = NA)
  
  bigSmallLabelDf$colCode[bigSmallLabelDf$bigSmall == 'Big'] =  "#E41A1CFF"
  bigSmallLabelDf$colCode[bigSmallLabelDf$bigSmall == 'Small'] = "#377EB8FF"
  
  
  ha2 = rowAnnotation(Group = anno_simple(bigSmallLabelDf$patientId, pch = 16, 
                                          pt_gp = gpar(col = bigSmallLabelDf$colCode),
                                          gp = gpar(fill = 'white', col = NA)))
  
  resDisMigrRm = resDisMigr
  
  
}
{
  patientNameListAfterOrder = c("TCGA-12-3649", 
                                "TCGA-19-5960", 
                                "TCGA-02-0033",
                                "TCGA-02-0085", 
                                "TCGA-06-1802", 
                                "TCGA-27-1830", 
                                "TCGA-06-0238",
                                "TCGA-02-0003",
                                "TCGA-12-1096",
                                "TCGA-12-1599",
                                "TCGA-06-0188", 
                                "TCGA-19-5951",
                                "TCGA-27-1838", 
                                "TCGA-08-0392", 
                                "TCGA-06-0139", 
                                "TCGA-27-1835", 
                                "TCGA-19-5955", 
                                "TCGA-12-1097", 
                                "TCGA-08-0352",
                                "TCGA-19-1388",
                                "TCGA-19-2620", 
                                "TCGA-76-6193",
                                "TCGA-14-0789", 
                                "TCGA-08-0390",
                                "TCGA-08-0524", 
                                "TCGA-14-1401", 
                                "TCGA-06-0142", 
                                "TCGA-06-0173",
                                "TCGA-08-0349", 
                                "TCGA-76-6664", 
                                "TCGA-19-5954", 
                                "TCGA-06-0133",
                                "TCGA-14-1396",
                                "TCGA-14-1453", 
                                "TCGA-06-0166")
}
rowSplitOrder = data.frame(now = rev(as.numeric(factor(levels = rownames(resDisMigr), patientNameListAfterOrder))),
                           earlyLate = c(rep('High', 22), rep('Low',13)))
ha2_t = columnAnnotation(Group = anno_simple(bigSmallLabelDf$patientId, pch = 16, 
                                             pt_gp = gpar(col = bigSmallLabelDf$colCode),
                                             gp = gpar(fill = 'white', col = NA)))

Heatmap(t(resDisMigrRm)[8:1,],
        name = "Probability",
        col = my_col_fun,
        row_order = order(migrOrderList2),
        row_split = factor(c(rep('Low',4), rep('High', 4)),
                           levels = c('Low', 'High')),
        
        row_title = "Migration rate",
        show_row_dend = F, row_names_side = 'left',
        
        column_order = as.numeric(factor(levels = rownames(resDisMigr), 
                                         patientNameListAfterOrder)),
        column_split = factor(rowSplitOrder$earlyLate[order(rowSplitOrder$now)],
                              levels = c('Low', 'High')),
        
        column_names_gp = gpar(fontsize = 7),
        show_column_dend = F,
        cluster_columns = FALSE,
        rect_gp = gpar(col = 'black',
                       lwd = 0.5),
        top_annotation = ha2_t,
        
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(t(resDisMigrRm)[8:1,][i, j] == medianList[[j]])
            grid.points(x,y, pch = 21, size = unit(2, "mm"))
        },
        
        
)

