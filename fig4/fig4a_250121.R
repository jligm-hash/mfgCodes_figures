
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
    resDist.tmp = (sim.scale[,i] - target.scale[,i])^2
    resDist = resDist + resDist.tmp
  }
  ntol = ceiling(nrow(sim.scale) * tol)
  cutOffN = sort(resDist)[ntol]
  postRindex = resDist <= cutOffN
  resTmpDf = aggregate(rank(-resDist[postRindex]), 
                       by=list(label = modelLabel1[postRindex]), FUN=sum)
  resTmpDf = resTmpDf[order(resTmpDf$label),]
  
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
rm(resPostAllOnly)
for (i in 1:length(patientName)) {
  pb$tick()
  resPost = myAbcPr(setTarget = obsDat[i,1:7], 
                    setSimulation = simDatOnly[,c('num', paste0('P',1:6))], 
                    setLabel = simDatOnly$migr, 
                    setTol =0.1) 
  
  tmpPostRes = simDatOnly[resPost, ]
  tmpPostRes$patientId = patientName[i]
  
  if(i == 1){
    resPostAllOnly = tmpPostRes
  }else{
    resPostAllOnly = rbind(resPostAllOnly, tmpPostRes)
  }
}

resPostAllOnly$invisibleLesionNum = resPostAllOnly$lesionNum - resPostAllOnly$num
resPostAllOnly[, c('invisibleLesionNum', 'patientId')]

resMeanDf = data.frame(aggregate(resPostAllOnly[, c('invisibleLesionNum', 'patientId')], 
                                 by=list(resPostAllOnly$patientId), FUN=mean))

resMeanDf$patientId = resMeanDf$Group.1

sgbmMgbmLabelPwd = './mgbmSgbmLabel.csv'
sgbmMgbmLabelDf = data.frame(fread(sgbmMgbmLabelPwd))
sgbmMgbmLabelDf = na.omit(sgbmMgbmLabelDf) # filter out the patients without MRI

dim(sgbmMgbmLabelDf)
resMeanDf$group = sgbmMgbmLabelDf$bigSmall[match(resMeanDf$patientId,
                                                 sgbmMgbmLabelDf$ID)]
set.seed(100)
ggplot(resMeanDf, mapping = aes(x = group,
                                y = invisibleLesionNum,
                                group = group,
                                col = group)) +
  geom_boxplot(width= 0.3, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.5, alpha = 0.8) +
  stat_compare_means(comparisons = list(c(1, 2)), method = 'wilcox.test') +
  theme_classic() +
  scale_color_brewer(palette = 'Set1') +
  theme(text = element_text(size = 15)) +
  guides(col = 'none') +
  labs(x = '',
       y = '# invisable lesions')

