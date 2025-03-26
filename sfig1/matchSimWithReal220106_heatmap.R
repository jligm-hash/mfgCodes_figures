library(data.table)
library(progress)
library(abc)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

myListToCol = function(x){
  ColDf = data.frame(col = table(x))
  colnames(ColDf) = c('Value', 'Freq')
  ColDf$Col = rainbow(dim(ColDf)[1])
  return(ColDf$Col[match(x, ColDf$Value)])
}

myListToLog10 = function(x){
  return(ifelse(x==0, 0, log10(x)))
}
pwd = '/Users/jiabao/Documents/workspace/mfg/labMeeting220126/heatmap'

setwd(pwd)

load(file = 'resReadSimv115.RData')
dim(resSimAll)

resSimAll = resSimAll[resSimAll$migr<1e-4,]
dim(resSimAll)

unique(resSimAll$migr)
unique(resSimAll$select)

resSimAll = resSimAll[resSimAll$select == 0.05, ] # according to Nowak's results
dim(resSimAll)

colnames(resSimAll)
realSizeDf = data.frame(fread('realTCGA211017.csv'))
dim(realSizeDf)
realSizeDf[is.na(realSizeDf)] = 'Unknown'
dim(realSizeDf)
resSimAll
{
  

  
  obsDat = realSizeDf[,c('Num', paste0('T',1:6))]
  patientName = realSizeDf$ID
  obsDat[1,]
  
  simDat = resSimAll[resSimAll$num>=1,]
  simDat[is.na(simDat)] = 0
  
  hist(simDat$firstMigrNum)
  
  library(ggplot2)
  library(scales)
  
  ggplot(simDat, mapping = aes(x=firstMigrNum)) +
    geom_vline(aes(xintercept=median(firstMigrNum)),
               color="blue", linetype="dashed", size=1) +
    geom_histogram(aes(y=..density..), colour="black", fill="white")+
    geom_density(alpha=.2, fill="#FF6666") +
    geom_vline(aes(xintercept=median(firstMigrNum)),
               color="blue", linetype="dashed", size=1) +
    theme_minimal() + 
    scale_x_log10(label=trans_format("log10",math_format(10^.x)))+
    labs(x='# of cell niches at the first migration',
         y='Density')
  
  floor(log10(simDat$firstMigrNum+1))+1
  
  
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
    return(prop.table(table(modelLabel1[postRindex])))
  }
  
  
  simDatOnly = simDat
  simDatOnly[, paste0('P',1:10)] = 100 * simDatOnly[, paste0('P',1:10)]
  
  unique(simDatOnly$migr)
  unique(simDatOnly$select)
  unique(simDatOnly$distMigrDrvNum)
  
  simDatOnly$myTime = floor(log10(simDatOnly$firstMigrNum+1))+1
  unique(simDatOnly$myTime)
  
  barpltDf = data.frame(table(floor(log10(simDatOnly$firstMigrNum+1))+1))
  barpltDf$per = signif(barpltDf$Freq/sum(barpltDf$Freq),digits = 2)
  colnames(barpltDf)[1] = 'log10_number'
  
  ggplot(barpltDf, mapping = aes(x=log10_number, y=Freq, fill=log10_number)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25) +
    geom_text(aes(label=per), position=position_dodge(width=0.9), vjust=1.5) +
    labs(x = 'log10(number of cell niches)') +
    theme_minimal()
  
  
  simDatOnly[1, 1:20]
  obsDat[1,]
  
  pb = progress_bar$new(total = length(patientName), clear = FALSE)
  rm(resPostAllOnly)
  for (i in 1:length(patientName)) {
    pb$tick()
    resPost = myAbcPr(setTarget = obsDat[i,1:7], 
                      setSimulation = simDatOnly[,c('num', paste0('P',1:6))], 
                      setLabel = paste0(simDatOnly$migr, '.T', simDatOnly$myTime), 
                      setTol =0.01) 
    sumPost = data.frame(resPost)
    tmpPostRes = data.frame(t(sumPost$Freq),
                            row.names = patientName[i]) #paste0(round(obsDat$T2[i], digits = 2), '.',obsDat$Num[i], '.',patientName[i]))
    colnames(tmpPostRes) = sumPost$Var1
    if(i == 1){
      resPostAllOnly = tmpPostRes
    }else{
      resPostAllOnly = rbind(resPostAllOnly, tmpPostRes)
    }
  }
  
  
  
  rsBarCol = cbind(obsDat[,paste0('T',1:6)])
  ha = rowAnnotation(RS = anno_barplot(rsBarCol,
                                       gp = gpar(fill = c('#8dd3c7','#bebada','#80b1d3','#fdb462','#b3de69','#ffffb3'),
                                                 col = NA))) # 2:7
  
  colBar = data.frame(x=1:6,
                      y=c('#8dd3c7','#bebada','#80b1d3','#fdb462','#b3de69','#ffffb3'))
  ggplot(colBar, mapping = aes(x=x,y=letters[1:6],
                               fill=as.factor(paste0('Lesion ',x)))) +
    geom_bar(stat = 'identity') +
    labs(fill = 'Lesion') +
    scale_fill_manual(values = colBar$y) +
    theme(legend.title = element_text(#size=10, 
      face="bold"))
  
  
  resDisMigr =  resPostAllOnly
  colnames(resDisMigr) # need to check the colname
  
  cn = colnames(resDisMigr)
  medianList = apply(resDisMigr, 1, max)
  
  maxNameList = c()
  for (i in 1:dim(resDisMigr)[1]) {
    tmpNameValue = names(resDisMigr)[resDisMigr[i,]==medianList[i]][1] # use the value of the 
    maxNameList = c(maxNameList, tmpNameValue)
  }
  
  
  
  
  
  
  my_col_fun = colorRamp2(c(0, 1), c("white", "red")) 
  
  resSelMigr = resDisMigr
  hcRes = hclust(dist(resSelMigr)) # use the result of the latest resSelMigr
  
  subClust = cutree(hcRes, k = 2)
  hcRes$labels[subClust==1]
  length(hcRes$labels[subClust==2])
  
  subClust = cutree(hcRes, k = 2)
  
  
  tmpIdName = lapply(names(subClust), strsplit, split='-')
  patientIdList = c()
  for (i4patient in 1:length(tmpIdName)) {
    tmpIdNameSplit = tmpIdName[[i4patient]][[1]]
    tmpId4patient = paste(tmpIdNameSplit[(length(tmpIdNameSplit)-2):length(tmpIdNameSplit)], 
                          collapse = '-')
    patientIdList = c(patientIdList, tmpId4patient)
  }
  
  clustDf = data.frame(clusterNo = subClust,# as.factor(subClust), 
                       orginal = names(subClust),
                       id = patientIdList, # make sure id == original
                       row.names = patientIdList)
  group1PatientNam = clustDf$id[clustDf$clusterNo==1]
  group2PatientNam = clustDf$id[clustDf$clusterNo==2]
  length(group1PatientNam)
  length(group2PatientNam)
  
  clustDf2 = clustDf
  clustDf2 = clustDf2[match(patientName,clustDf2$id),] # reorder the list by RS 1
  clustDf2$RS1 = realSizeDf$T1[match(clustDf2$id, realSizeDf$ID)]
  clustDf2 = clustDf2[order(clustDf2$clusterNo),]
  clustDf2[clustDf2$clusterNo == 1,] = clustDf2[clustDf2$clusterNo == 1,][order(clustDf2[clustDf2$clusterNo == 1,]$RS1),]
  clustDf2[clustDf2$clusterNo == 2,] = clustDf2[clustDf2$clusterNo == 2,][order(clustDf2[clustDf2$clusterNo == 2,]$RS1),]
  rownames(clustDf2) = clustDf2$id
  
  rsBarCol = realSizeDf[match(clustDf2$id, realSizeDf$ID),c(paste0('T',1:6))] #cbind(obsDat[,paste0('T',1:6)])
  
  
  anno4realSizeDf = realSizeDf[match(clustDf2$id, realSizeDf$ID),]
  
  rsBarCol = realSizeDf[match(rownames(resDisMigr), realSizeDf$ID),c(paste0('T',1:6))] #cbind(obsDat[,paste0('T',1:6)])
  anno4realSizeDf = realSizeDf[match(rownames(resDisMigr), realSizeDf$ID),]
  patientClinicalAnno=rowAnnotation(RS = anno_barplot(rsBarCol, 
                                                      gp = gpar(fill = c('#8dd3c7','#bebada','#80b1d3','#fdb462','#b3de69','#ffffb3'),
                                                                col = NA)),
                                    Group=anno4realSizeDf$Group,
                                    PIK3CA=anno4realSizeDf$pik3ca,
                                    Gender=anno4realSizeDf$Gender,
                                    IDH.status=anno4realSizeDf$IDH.status,
                                    MGMT.pro.status=anno4realSizeDf$MGMT.promoter.status,
                                    Original.Subtype=anno4realSizeDf$Original.Subtype,
                                    gp = gpar(col = 'white', lwd=3),
                                    gap = unit(2, "points"),
                                    col=list(Group = structure(c('#33a02c','#fb9a99'), names=c("Big",   "Small")),
                                             PIK3CA=structure(c('#a6cee3','grey90','#1f78b4'), names=c("WT",      "Unknown", "MUT")),
                                             Gender=c('Female'='#beaed4', 'Male'='#fdc086'),
                                             IDH.status=structure(c('#e31a1c', 'grey90'), names=c("WT",      "Unknown")),
                                             MGMT.pro.status=structure(c('grey90', '#ccebc5','#ffed6f'), names=c("Unknown", "Methylated",      "Unmethylated")),
                                             Original.Subtype=structure(c('#80b1d3','#fdb462','#b3de69','#fccde5','grey90'), names=c("Neural",      "Mesenchymal", "Classical",   "Proneural",   "Unknown"))
                                    )
  )
  
  Heatmap(resDisMigr,
          name = "postPr",
          col = my_col_fun,
          show_column_dend = F,
          show_row_dend = F,
          rect_gp = gpar(col = "grey90", lwd = 0.5),
          right_annotation = patientClinicalAnno,
          column_title = "Driver number at the time of migration dissemination",
          row_names_gp = gpar(fontsize = 7),
          show_column_names = FALSE,
          cell_fun = function(j, i, x, y, width, height, fill) {
            if(resDisMigr[i, j] == medianList[[i]])
              grid.points(x,y, pch = 21, size = unit(2, "mm"))
          },
          bottom_annotation = HeatmapAnnotation(text = anno_text(cn, rot = 90,just = "right"))
  )
  
  
  
  resHclust = hclust(dist(resDisMigr))
  resHclust$labels
  cutree(resHclust, k = 4) # 14, 2, 3
  
  resCutree = cutree(resHclust, k = 4) 
  resCutree[resCutree == 4] = 1 # change the label of 06-0166
  
  resCutree[names(resCutree) == 'TCGA-19-1388'] = 1
  
  
  cluster1List = names(resCutree)[resCutree==1]
  cluster2List = names(resCutree)[resCutree==2]
  cluster3List = names(resCutree)[resCutree==3]
  
  
  
  timeListFromMaxNameList = c()
  for (tmpMaxName in maxNameList) {
    timeListFromMaxNameList = c(timeListFromMaxNameList, as.numeric(unlist(strsplit(tmpMaxName, '.T'))[1]))
  }
  
  patientNameListAfterOrder = c(rev(names(medianList)[names(medianList) %in% cluster2List][order(timeListFromMaxNameList[names(medianList) %in% cluster2List], 
                                                                                                 medianList[names(medianList) %in% cluster2List])]),
                                rev(names(medianList)[names(medianList) %in% cluster1List][order(timeListFromMaxNameList[names(medianList) %in% cluster1List], 
                                                                                                 medianList[names(medianList) %in% cluster1List])]),
                                rev(names(medianList)[names(medianList) %in% cluster3List][order(timeListFromMaxNameList[names(medianList) %in% cluster3List], 
                                                                                                 medianList[names(medianList) %in% cluster3List])]))
  
  
  
  
  
  
  migrOrderList = c()
  timeOrderList = c()
  colNameList = colnames(resDisMigr)
  for (i in 1:length(colNameList)) {
    tmpMigr = unlist(strsplit(colNameList[i], '.T'))[1]
    tmpTime = unlist(strsplit(colNameList[i], '.T'))[2]
    migrOrderList = c(migrOrderList, as.numeric(tmpMigr))
    timeOrderList = c(timeOrderList, as.numeric(tmpTime))
  }
  
  Heatmap(resDisMigr,
          name = "postPr",
          col = my_col_fun,
          cluster_columns = FALSE,
          column_order = order(timeOrderList, migrOrderList),
          
          cluster_rows = FALSE,
          row_order = rev(as.numeric(factor(levels = rownames(resDisMigr), patientNameListAfterOrder))),
          show_column_dend = F,
          show_row_dend = F,
          rect_gp = gpar(col = "grey90", lwd = 0.5),
          right_annotation = ha,
          column_title = "Migration probabilities",
          row_names_gp = gpar(fontsize = 7),
          show_column_names = FALSE,
          cell_fun = function(j, i, x, y, width, height, fill) {
            if(resDisMigr[i, j] == medianList[[i]])
              grid.points(x,y, pch = 21, size = unit(2, "mm"))
          },
          bottom_annotation = HeatmapAnnotation(text = anno_text(cn, rot = 80,just = "right"))
  )
  
  
  

  
  resDisMigrRm = resDisMigr[,timeOrderList>2]
  migrOrderList2 = c()
  timeOrderList2 = c()
  colNameList2 = colnames(resDisMigrRm)
  for (i in 1:length(colNameList2)) {
    tmpMigr = unlist(strsplit(colNameList2[i], '.T'))[1]
    tmpTime = unlist(strsplit(colNameList2[i], '.T'))[2]
    migrOrderList2 = c(migrOrderList2, as.numeric(tmpMigr))
    timeOrderList2 = c(timeOrderList2, as.numeric(tmpTime))
  }
  
  Heatmap(resDisMigrRm,
          name = "Probability",
          col = my_col_fun,
          cluster_columns = FALSE,
          column_order = order(timeOrderList2, migrOrderList2),
          
          cluster_rows = FALSE,
          row_order = rev(as.numeric(factor(levels = rownames(resDisMigr), patientNameListAfterOrder))),
          show_column_dend = F,
          show_row_dend = F,
          rect_gp = gpar(col = "grey90", lwd = 0.5),
          right_annotation = ha,
          column_title = "Migration rate with first migration time",
          row_names_gp = gpar(fontsize = 7),
          show_column_names = FALSE,
          cell_fun = function(j, i, x, y, width, height, fill) {
            if(resDisMigrRm[i, j] == medianList[[i]])
              grid.points(x,y, pch = 21, size = unit(2, "mm"))
          },
          bottom_annotation = HeatmapAnnotation(text = anno_text(cn[timeOrderList>2], rot = 85,just = "right"))
  )
  
  
  
}

