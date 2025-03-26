

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
library(clusterProfiler)
library(enrichplot)

source('./plotfgsea.R')
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
gema0 = read.delim('./gbm_tcga_data_mrna_affymetrix_microarray.txt', check.names = F)
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
gema_dif = gema_dif[order(gema_dif$log2FC,  decreasing = T), ]
fcrnk = gema_dif$log2FC
names(fcrnk) = gema_dif$gene  
fcrnk = fcrnk[!is.na(fcrnk)]
fcrnk = fcrnk[!duplicated(names(fcrnk))]
hmpathways <- gmtPathways("h.all.v7.0.symbols.gmt")
fgseaRes <- fgsea(pathways = hmpathways,  stats = fcrnk,  minSize=15, maxSize=500)
c2pathways = gmtPathways("c2.all.v7.5.1.symbols.gmt")
fgseaResc2 <- fgsea(pathways = c2pathways,  stats = fcrnk,  minSize=5, maxSize=500)
kpathways = read.gmt('c2.all.v7.5.1.symbols.gmt')
gsea = GSEA(fcrnk,  TERM2GENE  = kpathways,  eps = 1e-300)

g = gsea@result
g = g[order(g$NES<0,  g$p.adjust), ]
p <- gseaplot2(gsea, geneSetID = g$ID[grepl('VERHAAK_G', g$ID)], 
               color = c('#FF9843', '#FFDD95', '#86A7FC', '#3468C0'), rel_heights = c(1.2, 0.2, 0.6), 
               ES_geom = 'line', pvalue_table = F,  base_size = 12)
p[[1]] = p[[1]] + guides(color = 'none') +
  theme(text = element_text(size = 15))
p

