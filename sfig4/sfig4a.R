
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

setFilePwd = './'
source(file.path(setFilePwd,   'plotfgsea.R'))
idh = read.delim(file.path(setFilePwd,  'TCGA_GBM_MRI259_IDH.txt'))
mdp = read.delim(file.path(setFilePwd,  'gbm_tcga_data_clinical_patient.txt'),  comment.char = '#',  na.strings = c('', '[Not Available]'))
mds = read.delim(file.path(setFilePwd,  'gbm_tcga_data_clinical_sample.txt'),  comment.char = '#')
fg = read.delim(file.path(setFilePwd,  'gbm_tcga_data_clinical_patient.txt'), comment.char = '#')

cl = read.csv(file.path(setFilePwd,  'mgbmSgbmLabel.csv'))
cl = cl[which(cl$modality == 'MR'), ]
cl$IDH = idh$IDH.status[match(cl$ID,  idh$Patient)]

l2b = cl$ID[cl$bigSmall == 'Big']
l2s = cl$ID[cl$bigSmall == 'Small']
l1 = cl$ID[cl$bigSmall == 'SGBM'&cl$IDH == 'WT']
gema0 = read.delim(file.path(setFilePwd,  'gbm_tcga_data_mrna_affymetrix_microarray.txt'), check.names = F)
gema = gema0[,  names(gema0) %in% c('Hugo_Symbol', paste0(c(l2b, l2s), '-01'))]
gema_dif = data.frame(gene = gema$Hugo_Symbol,  L2l = 0,  L2s = 0,  pval_w =1,  pval_t =1)
for (i in 1:nrow(gema)){
  tmp = as.data.frame(t(gema[i, -1]))
  names(tmp)[1]='gene'
  tmp$group = ifelse(rownames(tmp) %in% paste0(l2b, '-01'), 'L2l', 'L2s')
  gema_dif$L2l[i] = median(tmp$gene[tmp$group=='L2l'])
  gema_dif$L2s[i] = median(tmp$gene[tmp$group=='L2s'])
  gema_dif$pval_w[i] = wilcox.test(gene ~ group,  data= tmp)$p.value
  gema_dif$pval_t[i] = t.test(gene ~ group,  data= tmp)$p.value
}
gema_dif$log2FC = gema_dif$L2l-gema_dif$L2s
gema_dif$type = ifelse(gema_dif$pval_t < 0.05 & gema_dif$log2FC > 0.5, 'up', 
                       ifelse(gema_dif$pval_t < 0.05 & gema_dif$log2FC <  -0.5, 'down', 'nosig'))
gema_dif$lab = ifelse((gema_dif$pval_t < 0.05 & abs(gema_dif$log2FC) > 1) | (gema_dif$pval_t < 0.001 ), 
                      gema_dif$gene, NA)
ggplot() + 
  
  geom_vline(xintercept = c(-0.5, 0.5), lty=2, col='#cccccc') + 
  geom_hline(yintercept = -log10(0.05), lty=2, col='#cccccc') + 
  geom_point(gema_dif[gema_dif$type == 'nosig', ],  mapping = aes(x = log2FC,  y = -log10(pval_t), color = type), 
             show.legend = F,
             size = 0.2, alpha = 0.6) + 
  geom_point(gema_dif[gema_dif$type != 'nosig', ],  mapping = aes(x = log2FC,  y = -log10(pval_t), color = type), 
             show.legend = F,
             size = 0.5, alpha = 0.6) + 
  geom_text_repel(gema_dif,  mapping = aes(x = log2FC,  y = -log10(pval_t), label = lab),  size = 2, 
                  min.segment.length = 0, segment.size=0.25) + 
  scale_color_manual(values = c('blue', '#cccccc', 'red')) + 
  scale_size_manual(values =c(1, 0.2, 1)) + 
  scale_y_continuous(expand = c(0, 0), 
                     limits=c(0, 4)) + 
  theme_classic() + 
  labs(x = "L2-large - L2-small",  y = '-log10(P-value)') +
  theme(text = element_text(size = 15))

