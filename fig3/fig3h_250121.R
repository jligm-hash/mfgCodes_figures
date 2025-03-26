
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
idh = read.delim('./TCGA_GBM_MRI259_IDH.txt')
mdp = read.delim('./gbm_tcga_data_clinical_patient.txt',   comment.char = '#',   na.strings = c('',  '[Not Available]'))
mds = read.delim('./gbm_tcga_data_clinical_sample.txt',   comment.char = '#')
fg = read.delim('./gbm_tcga_data_clinical_patient.txt',  comment.char = '#')

cl = read.csv('./mgbmSgbmLabel.csv')
cl = cl[which(cl$modality == 'MR'),  ]
cl$IDH = idh$IDH.status[match(cl$ID,   idh$Patient)]

l2b = cl$ID[cl$bigSmall == 'Big']
l2s = cl$ID[cl$bigSmall == 'Small']
l1 = cl$ID[cl$bigSmall == 'SGBM'&cl$IDH == 'WT']
cna0 = read.delim('./gbm_tcga_data_cna.txt',  
                  check.names = F)
gns = c('EGFR',  'PDGFRA',  'CDK4',  'MDM2',  'PIK3CA',  'CDKN2A')
cna = cna0[cna0$Hugo_Symbol %in%gns,   names(cna0) %in% c('Hugo_Symbol',  paste0(c(l2b,  l2s),  '-01'))]
mut0 = read.delim('./gbm_tcga_data_mutations.txt')
mut01 = mut0[mut0$Tumor_Sample_Barcode %in% c('Hugo_Symbol',  paste0(c(l2b,  l2s),  '-01')),  ]
mut01 = mut01[mut01$Variant_Classification!='Silent',   ]
nmut = as.data.frame(table(mut01$Tumor_Sample_Barcode))
nmut$type = ifelse(nmut$Var1 %in% paste0(l2b,  '-01'),   'L2-large',   'L2-small')
mds = mds[mds$PATIENT_ID %in% c(l2b, l2s), ]
mds$type = ifelse(mds$PATIENT_ID %in% l2b, 'L2-large', 'L2-small')

ggplot(mds,  aes(x = type ,  y = TMB_NONSYNONYMOUS,
                 col = type))  + 
  geom_boxplot(outlier.shape = NA, width = 0.35) + 
  geom_beeswarm(cex = 3, size = 0.5, alpha = 0.8) +
  stat_compare_means(comparisons = list(c('L2-large', 'L2-small')), label = 'p.format')  + 
  scale_color_manual(values = c('#e31a1c', '#377eb8')) + 
  theme_classic()  + 
  labs(x = '',  y = 'somatic mutation burden') +
  guides(col = 'none') +
  theme(text = element_text(size = 15))

