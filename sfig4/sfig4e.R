
library(survival)
library(survminer)
library(ggplot2)

cl = read.delim('GLASS2019_pairedGBM_with_recur_locations.txt')
cl$PFS_censor =1
my.Surv <- with(cl,Surv(time = Surgical.Interval..mths..1, event = PFS_censor))
expr.surv <- survminer::surv_fit(my.Surv ~ LocalorDistal, data = cl)
log.rank <- survdiff(my.Surv ~ LocalorDistal, rho = 1, data = cl)
ggsurvplot(fit = expr.surv, censor = T, conf.int = F, legend = c(0.7,0.9),surv.median.line='hv',
           surv.scale = "percent", ylab = "Probability", xlab = "Surgical interval (month)",
           font.legend = 8, risk.table = F, palette = "Set1", ggtheme = theme_classic(),data = cl)

setGeneSet1 = c("COL13A1", "DPP4", "SERPINE1", "ANG", "CHST2", "AGPAT2", "GADD45B", "MMP19", "POU2F2", "EMC3",#"TMEM111",
                "HEBP1", "REXO2", "CA9", "NDUFC1", "GKN1", "IGFBP6", "CHCHD7", "MTCP1", "PSORS1C1", "MYOF",#"FER1L3",
                "TNFRSF12A",  "UPP1", "RAB13", "MTMR11", "GUSB",  "SNX10", "CCNO", "CAST")
setGeneSet2 = c("DPP4", "FAM26B", "ZCCHC6", "TMBIM1", "SERPINE1", "NRP2", "CAST", "FER1L3", "OSBPL3",
                "COL13A1", "GADD45B", "BCAT1", "tcag7.1314", "MMP19", "GUSB", "SLC22A18", "CLEC5A", 
                "GALNS", "SQSTM1", "ANG", "AGPAT2", "FKBP2", "POU2F2", "ALK", "MTMR11", "PLK3", "RAB13",
                "TNFRSF12A", "NACAP1", "SLC22A4", "MCFD2", "PTPN14", "MREG", "CNIH3", "ITPR3", "NPEPL1",
                "KLHL21", "LGTN", "UPP1", "HEBP1", "NQO2")
ge = read.delim('/Users/jiabao/Library/CloudStorage/Dropbox/Jiabao/mfgManuscript221028/Fig4/PangliomaPairs/glass20220830_RNAseq_RPKM_all_samples.tsv', check.names = F)
gei = ge[,substr(names(ge),1,15) %in% cl$Sample.barcode]
ger = ge[,substr(names(ge),1,15) %in% cl$Sample.barcode.1]

gei_set = gei[rownames(gei) %in% setGeneSet1,]
gei_set_z = as.data.frame(t(scale(t(log2(gei_set+1)))))
gei_set_md = as.data.frame(apply(gei_set_z,2,median))
names(gei_set_md)[1] = 'Score'
gei_set_md$Case = substr(rownames(gei_set_md),1,12)
gei_set_md$localDistal = cl$LocalorDistal[match(gei_set_md$Case, cl$Case.barcode)]

gei_set_md$surgical_interval = cl$Surgical.Interval..mths..1[match(gei_set_md$Case, cl$Case.barcode)]
gei_set_md$group = ifelse(gei_set_md$Score>median(gei_set_md$Score),'high','low')

gei_set_md$PFS_censor =1
my.Surv <- with(gei_set_md,Surv(time = surgical_interval, event = PFS_censor)) 
expr.surv <- survminer::surv_fit(my.Surv ~ group, data = gei_set_md)
log.rank <- survdiff(my.Surv ~ group, rho = 1, data = gei_set_md)
ggsurvplot(fit = expr.surv, censor = T, conf.int = F, legend = c(0.7,0.9),surv.median.line='hv',
           surv.scale = "percent", ylab = "Probability", xlab = "Surgical interval (month)",
           font.legend = 8, risk.table = F, palette = "Set1", ggtheme = theme_classic(),data = gei_set_md)

table(gei_set_md$group, gei_set_md$localDistal)
gei_set_md$group = factor(gei_set_md$group,levels = c('low','high'))
barplot(100*prop.table(table(gei_set_md$group, gei_set_md$localDistal), margin = 2))
pltDf = data.frame(reshape2::melt(table(gei_set_md$group, gei_set_md$localDistal)))

colnames(pltDf) = c('mGBscoreGroup', 'position', 'value')

myValue2prop = function(x){
  return(x/sum(x))
}
pltDf$proportion[pltDf$position == 'Distal'] = myValue2prop(pltDf$value[pltDf$position == 'Distal'])
pltDf$proportion[pltDf$position == 'Local'] = myValue2prop(pltDf$value[pltDf$position == 'Local'])

pltDf$mGBscoreGroup = factor(pltDf$mGBscoreGroup, 
                             levels = c('high', 'low'))

pltDf$labelPosition = pltDf$proportion
pltDf$labelPosition[c(2, 4)] = pltDf$labelPosition[c(2, 4)] + 2* pltDf$labelPosition[c(1,3)] 
pltDf$labelPosition = pltDf$labelPosition/2
ggplot() +
  geom_bar(pltDf, mapping = aes(x = position,
                                y = proportion,
                                fill = mGBscoreGroup),
           stat = 'identity') +
  
  scale_fill_manual(values = c('#E0E0E0', '#3C3C3C')) +
  theme_classic() +
  
  theme(text = element_text(size = 15),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
  ) +
  
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = '',
       y = 'Proportion of patients',
       fill = 'mGB score') +
  guides(col = 'none')

