
library(ggbeeswarm)
library(ggplot2)
library(progress)
library(data.table)

rm(list = ls())
load('./validationCurve4lesion2cutoff0.RData')
sim4patientWith3and3MLesionCombine = sim4patientWith3and3MLesionCombine[sim4patientWith3and3MLesionCombine$lesionNum > 2, ]
sim4migr = sim4patientWith3and3MLesionCombine
cellNumArr = sort(unique(unlist(sim4migr[, paste0('T', 1:6)])))
simCutoffArr = sort(c(0, 1e6, sample(x = 1:length(cellNumArr), size = 1e4)))
pb = progress_bar$new(total = length(simCutoffArr))

rm('res4rocCurveDf')

realVolumePwd = './realVolume_240412.csv'
relativeSizePwd = './realTCGA211017.csv'
realVolumeDf = data.frame(fread(realVolumePwd))
relativeSizeDf = data.frame(fread(relativeSizePwd))

realVolumeDf$patientLabel = ifelse(realVolumeDf$nTumor == 2, '2', '2+')
relativeSizeDf$patientLabel = ifelse(relativeSizeDf$Num == 2, '2', '2+')
trueLabelArr = realVolumeDf$patientLabel
valueArr4roc = realVolumeDf$T2
valueArr4roc = realVolumeDf$T1

trueLabelArr[order(valueArr4roc)]

rm('controlRocDf')
for (tmpCutoff in sort(valueArr4roc)) {
  predictArr = ifelse(valueArr4roc<tmpCutoff, '2+', '2')
  
  contingencyTab = data.frame(table(predictArr, trueLabelArr))
  
  colnames(contingencyTab) = c('simLabel', 'trueLabel', 'Freq')
  contingencyTab
  
  
  if(sum(contingencyTab$simLabel == '2') > 0){
    value4x = 1-contingencyTab$Freq[contingencyTab$simLabel == '2' & contingencyTab$trueLabel == '2']/(sum(contingencyTab$Freq[contingencyTab$trueLabel == '2']))
  }else{
    value4x = 1
  }
  
  
  
  if(sum(contingencyTab$simLabel == '2+') > 0){
    value4y = contingencyTab$Freq[contingencyTab$simLabel == '2+' & contingencyTab$trueLabel == '2+']/(sum(contingencyTab$Freq[contingencyTab$trueLabel == '2+']))
  }else{
    value4y = 0
  }
  
  
  
  tmpDf = data.frame(cutoff = tmpCutoff,
                     value4x = value4x,
                     value4y = value4y)
  
  if(!exists('controlRocDf')){
    controlRocDf = tmpDf
  }else{
    controlRocDf = rbind(controlRocDf, tmpDf)
  }
}
myCalculateControlRocDf = function(trueLabelArr, valueArr4roc, revLabel = T){
  
  trueLabelArr = trueLabelArr
  valueArr4roc = valueArr4roc
  
  rm('controlRocDf_tmp')
  
  for (tmpCutoff in sort(valueArr4roc)) {
    
    print(tmpCutoff)
    if(revLabel){
      predictArr = ifelse(valueArr4roc<tmpCutoff, '2+', '2')
    }else{
      predictArr = ifelse(valueArr4roc<tmpCutoff, '2', '2+')
    }
    
    contingencyTab = data.frame(table(predictArr, trueLabelArr))
    colnames(contingencyTab) = c('simLabel', 'trueLabel', 'Freq')
    
    if(sum(contingencyTab$simLabel == '2') > 0){
      value4x = 1-contingencyTab$Freq[contingencyTab$simLabel == '2' & contingencyTab$trueLabel == '2']/(sum(contingencyTab$Freq[contingencyTab$trueLabel == '2']))
    }else{
      value4x = 1
    }
    
    if(sum(contingencyTab$simLabel == '2+') > 0){
      value4y = contingencyTab$Freq[contingencyTab$simLabel == '2+' & contingencyTab$trueLabel == '2+']/(sum(contingencyTab$Freq[contingencyTab$trueLabel == '2+']))
    }else{
      value4y = 0
    }
    tmpDf = data.frame(cutoff = tmpCutoff,
                       value4x = value4x,
                       value4y = value4y)
    
    print(tmpDf)
    
    if(!exists('controlRocDf_tmp')){
      controlRocDf_tmp = tmpDf
    }else{
      controlRocDf_tmp = rbind(controlRocDf_tmp, tmpDf)
    }
  }
  
  tmpDf = data.frame(cutoff = tmpCutoff,
                     value4x = 1,
                     value4y = 1)
  
  controlRocDf_tmp = rbind(controlRocDf_tmp, tmpDf)
  
  
  return(controlRocDf_tmp)
  
}
controlRocDfT1 = myCalculateControlRocDf(trueLabelArr = realVolumeDf$patientLabel, valueArr4roc = realVolumeDf$T1, revLabel = T)
controlRocDfT2 = myCalculateControlRocDf(realVolumeDf$patientLabel, realVolumeDf$T2, revLabel = T)
controlRocDfP1 = myCalculateControlRocDf(relativeSizeDf$patientLabel, relativeSizeDf$T1, revLabel = T)

controlRocDfP2 = myCalculateControlRocDf(relativeSizeDf$patientLabel, relativeSizeDf$T2, revLabel = T)
load('./rocCurve_lesionNumLargeThan2_240503.RData')

library(ggplotify)
library(plotly)
ggplot() +
geom_abline(slope = 1, intercept = 0, col = 'grey90', linetype = 'dashed') +
  geom_line(controlRocDfP2, mapping = aes(x = value4x,
                                          y = value4y
                                          ),
            linewidth = 0.6, col = '#FFBF78', alpha = 1) +
  geom_line(res4rocCurveDf, mapping = aes(x = value4x,
                                          y = value4y
                                          ),
            linewidth = 1, col = '#B52B65') +
  geom_text() +
  theme_classic() +
  coord_equal() +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = '1 - Specificity',
       y = 'Sensitivity') +
  theme(text = element_text(size = 15))

