
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

source('./plotfgsea.R')
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

gns2 = c('IDH1', 'PDGFRA', 'RB1', 'NF1', 'EGFR', 'TP53', 'PTEN', 
         'PIK3CB', 'PIK3CG', 'PIK3CA', 'PIK3R1')
mut = mut0[mut0$Hugo_Symbol %in% gns2 & mut0$Tumor_Sample_Barcode %in% c('Hugo_Symbol', paste0(c(l2b, l2s), '-01')), ]
mut = mut[mut$Variant_Classification!='Silent', ]

tmpEgfrDf = mut[mut$Hugo_Symbol == 'EGFR', ]
gain7loss10labelDf4mgbm = gain7loss10labelDf[gain7loss10labelDf$Case %in% c(l2b, l2s), ]
cna4plt = cna
rownames(cna4plt) = cna4plt$Hugo_Symbol
cna4plt=cna4plt[, -1]
cna4plt = cna4plt[gns, paste0(c(l2b, l2s), '-01')]
rant = data.frame(type = c(rep('L2-large',  22),  rep('L2-small',  13)), 
                  row.names =paste0(c(l2b, l2s),  '-01')) 

tmp = as.data.frame(t(cna4plt))

tmp$type = rant$type[match(rownames(tmp),  rownames(tmp))]
tmp = tmp[order(tmp$type,  -tmp$EGFR, -tmp$PDGFRA, -tmp$CDK4, tmp$CDKN2A,  decreasing = F), ]
names(tmp)[names(tmp) == 'PIK3CA'] = 'PIK3CAamp'
names(tmp)[names(tmp) == 'PDGFRA'] = 'PDGFRAamp'
names(tmp)[names(tmp) == 'EGFR'] = 'EGFRamp'
names(tmp)[names(tmp) == 'CDKN2A'] = 'CDKN2Adel'
names(tmp)[names(tmp) == 'CDK4'] = 'CDK4amp'
names(tmp)[names(tmp) == 'MDM2'] = 'MDM2amp'
names(tmp)[names(tmp) == 'PTEN'] = 'PTENloss'

tmp[names(tmp)[grepl('amp', names(tmp))]] = ifelse(tmp[names(tmp)[grepl('amp', names(tmp))]] == 2, 2, 0)
tmp[names(tmp)[grepl('del', names(tmp))]] = ifelse(tmp[names(tmp)[grepl('del', names(tmp))]] == -2, -2, 0)
tmp[names(tmp)[grepl('loss', names(tmp))]] = ifelse(tmp[names(tmp)[grepl('loss', names(tmp))]] == -1, -1, 0)

tmp2 = matrix(data = NA, nrow = length(gns2), ncol = nrow(tmp))
tmp2 = as.data.frame(tmp2)
rownames(tmp2) = gns2; names(tmp2) = rownames(tmp)

for (i in 1:nrow(tmp2)){
  for (j in 1:ncol(tmp2)){
    muti = mut[mut$Hugo_Symbol == rownames(tmp2)[i], ]
    tmp2[i, j] = ifelse(!names(tmp2)[j] %in% mut0$Tumor_Sample_Barcode, NA, 
                        ifelse(!names(tmp2)[j] %in% muti$Tumor_Sample_Barcode, 0, 3))
  }
}

library(stringr)

resGain7loss10Arr = gain7loss10labelDf4mgbm$Chr.7.gain.Chr.10.loss[match(str_sub(colnames(tmp2), 1, 12),
                                                                         gain7loss10labelDf4mgbm$Case)]
resGain7loss10Arr[resGain7loss10Arr == "Gain chr 7 & loss chr 10"] = "4"
resGain7loss10Arr[resGain7loss10Arr == "No combined CNA"] = "0"

resGain7loss10Arr = as.numeric(resGain7loss10Arr)

tmp2['Gain7Loss10', ] = resGain7loss10Arr
pheatmap(rbind(tmp2, t(tmp[, -ncol(tmp)])), border_color = 'black',
         color = c('blue', '#a6cee3', 'white', '#fb9a99', 'red', 'darkgreen', '#D9E5FB'), # '#6495ED'
         annotation_colors = list(
           type = c('L2-large' = "#E41A1C",  
                    'L2-small' = "#377EB8") 
         ), 
         annotation_col = rant, cluster_cols = F, 
         cluster_rows = F, 
         gaps_col = c(22), 
         gaps_row = c(11, 12))
table(tmp$EGFRamp,  tmp$type)

fisher.test(table(tmp$EGFRamp,  tmp$type))
fisher.test(matrix(c(5, 8, 0, 7), nrow=2)) #egfr
fisher.test(matrix(c(5, 2, 4, 9), nrow=2)) #tp53
fisher.test(matrix(c(4, 9, 5, 2), nrow=2)) 

mutDf = t(tmp2)
mutDf = data.frame(mutDf)

rownames(mutDf) = substr(rownames(mutDf), 1, 12)
colnames(mutDf)

mutDf$subgroup = cl$bigSmall[match(rownames(mutDf),
                                   cl$ID)]

table(mutDf[, c('EGFR', 'subgroup')])
table(mutDf[, c('TP53', 'subgroup')])

fisher.test(table(mutDf[, c('EGFR', 'subgroup')]))
fisher.test(table(mutDf[, c('TP53', 'subgroup')]))


