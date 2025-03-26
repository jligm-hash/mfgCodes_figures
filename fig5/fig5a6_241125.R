
library(data.table)
library(ggplot2)

bestSimulationPwd = '/Users/jiabao/Documents/1_github/mfgCodes/bestSim4eachPatient/resAllBestSimNew35.csv'
bestSimulationDf = data.frame(fread(bestSimulationPwd))

realPatientPwd = '/Users/jiabao/Documents/1_github/mfgCodes/plotAbcHeatmap/realTCGA211017.csv'
realPatientDf = data.frame(fread(realPatientPwd))
mgbmPatientList = realPatientDf$ID

  {
  
  tmpPatientId = mgbmPatientList[28]
  
  
  tmpPlotDf = data.frame(Sim = unlist(bestSimulationDf[bestSimulationDf$p.ID == tmpPatientId, paste0('P', 1:6)]),
                         Real = unlist(realPatientDf[realPatientDf$ID == tmpPatientId, paste0('T', 1:6)]),
                         lesionId = factor(paste0('T', 1:6), levels = paste0('T', 1:6)))
  
  
  tmpPlotDfMelt = melt(tmpPlotDf, id = 'lesionId')
  colnames(tmpPlotDfMelt) = c('Lesion', 'Type', 'relativeSize')
  
  tmpPlotDfMelt$Type = factor(tmpPlotDfMelt$Type, levels = c('Real', 'Sim'))
  
  ggplot(tmpPlotDfMelt, mapping = aes(x = Lesion,
                                      y = relativeSize,
                                      col = Type,
                                      group = Type,
                                      shape = Type,
                                      fill = Type)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.5) +
    geom_point(mapping = aes(col = Type, group = Type, shape = Type), alpha = 0.5) +
    geom_line(mapping = aes(col = Type, group = Type), alpha = 0.5) +
    theme_classic() +
    labs(col = '',
         x = 'Lesion ID',
         y = 'Relative size (%)',
         shape = '',
         fill = '',
         title = tmpPatientId
         ) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 15)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,100)) +
    scale_fill_manual(values = c('#fc8d59', '#91bfdb')) +
    scale_color_manual(values = c('#fc8d59', '#91bfdb'))
  
  
}


