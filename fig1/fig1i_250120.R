
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(fgsea)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(fgsea)
library(clusterProfiler)
library(enrichplot)
library(stringr)
gain7loss10labelpwd = './gain7loss10label_240429.csv'

gain7loss10labelDf = data.frame(fread(gain7loss10labelpwd))

idh = read.delim('./TCGA_GBM_MRI259_IDH.txt')
mdp = read.delim('./gbm_tcga_data_clinical_patient.txt',  comment.char = '#',  na.strings = c('', '[Not Available]'))
mds = read.delim('./gbm_tcga_data_clinical_sample.txt',  comment.char = '#')
fg = read.delim('./gbm_tcga_data_clinical_patient.txt', comment.char = '#')

cl = read.csv('./mgbmSgbmLabel.csv')
cl = cl[which(cl$modality == 'MR'), ]
cl$IDH = idh$IDH.status[match(cl$ID,  idh$Patient)]

l2b = cl$ID[cl$bigSmall == 'Big']
l2s = cl$ID[cl$bigSmall == 'Small']
l1 = cl$ID[cl$bigSmall == 'SGBM'&cl$IDH == 'WT']

cna0 = read.delim('./gbm_tcga_data_cna.txt', 
                  check.names = F)
gns = c('EGFR', 'PDGFRA', 'CDK4', 'MDM2', 'PIK3CA', 'CDKN2A', 'PTEN')
cna = cna0[cna0$Hugo_Symbol %in%gns,  names(cna0) %in% c('Hugo_Symbol', paste0(c(l2b, l2s), '-01'))]
mut0 = read.delim('./gbm_tcga_data_mutations.txt')
mut01 = mut0[mut0$Tumor_Sample_Barcode %in% c('Hugo_Symbol', paste0(c(l2b, l2s), '-01')), ]
mut01 = mut01[mut01$Variant_Classification!='Silent',  ]
nmut = as.data.frame(table(mut01$Tumor_Sample_Barcode))
nmut$type = ifelse(nmut$Var1 %in% paste0(l2b, '-01'),  'L2-large',  'L2-small')
gns2 = c('PIK3CB', 'PIK3CG', 'PIK3CA', 'PIK3R1')
mut = mut0[mut0$Hugo_Symbol %in% gns2 & mut0$Tumor_Sample_Barcode %in% c('Hugo_Symbol', paste0(c(l1), '-01')), ]
mut = mut[mut$Variant_Classification!='Silent', ]
length(l1)
length(unique(substr(mut$Tumor_Sample_Barcode, 1, 12)))
num4totalSample = c(13, 7, 208)
num4pi3ksAlteration = c(9, 2, 15) 
label4group = c('l2large', 'l2small', 'sgbm')
df4barPlot = data.frame(group = label4group,
                        num = num4pi3ksAlteration,
                        totalNum = num4totalSample,
                        status = 'altered')

df4barPlot2 = data.frame(group = label4group,
                         num = num4totalSample - num4pi3ksAlteration,
                         totalNum = num4totalSample,
                         status = 'notAltered')
df4barPlotMerge = rbind(df4barPlot, df4barPlot2)
df4barPlotMerge$prop = df4barPlotMerge$num/df4barPlotMerge$totalNum

df4barPlotMerge$status = factor(df4barPlotMerge$status,
                                levels = rev(c('altered', 'notAltered')))

colnames(df4barPlotMerge)

ggplot(df4barPlotMerge, mapping = aes(x = group,
                                      y = prop,
                                      fill = status,
                                      label = num)) +
  geom_col(width = 0.5) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
  theme_classic() +
  theme(text = element_text(size = 15)) +
  scale_fill_manual(values = c('grey90', '#B1AFFF')) + 
  labs(x = '', y = "Proportion",
       fill = '')

prop.test(c(9, 2),
          c(13, 7))
fisher.test(matrix(c(4, 9, 5, 2), ncol = 2))
chisq.test(matrix(c(4, 9, 5, 2), ncol = 2))
mcnemar.test(matrix(c(4, 9, 5, 2), ncol = 2))

