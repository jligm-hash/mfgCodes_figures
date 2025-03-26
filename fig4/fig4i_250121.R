
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

library(data.table)
library(stringr)
source('./plotfgsea.R')
idh = read.delim('./TCGA_GBM_MRI259_IDH.txt')
mdp = read.delim('./gbm_tcga_data_clinical_patient.txt', comment.char = '#', na.strings = c('','[Not Available]'))
mds = read.delim('./gbm_tcga_data_clinical_sample.txt', comment.char = '#')
fg = read.delim('./gbm_tcga_data_clinical_patient.txt',comment.char = '#')

cl = read.csv('./mgbmSgbmLabel.csv')
cl = cl[which(cl$modality=='MR'),]
cl$IDH = idh$IDH.status[match(cl$ID, idh$Patient)]

l2b = cl$ID[cl$bigSmall=='Big']
l2s = cl$ID[cl$bigSmall=='Small']
l1 = cl$ID[cl$bigSmall=='SGBM'&cl$IDH=='WT']
gema0 = read.delim('gbm_tcga_data_mrna_affymetrix_microarray.txt',check.names = F)
gema = gema0[, names(gema0) %in% c('Hugo_Symbol',paste0(c(l2b,l2s),'-01'))]
gema_dif = data.frame(gene = gema$Hugo_Symbol, L2l = 0, L2s = 0, pval_w =1, pval_t =1)
sum(str_sub(names(gema0), 1, 12) %in% l2b)
sum(str_sub(names(gema0), 1, 12) %in% l2s)
for (i in 1:nrow(gema)){
  tmp = as.data.frame(t(gema[i,-1]))
  names(tmp)[1]='gene'
  tmp$group = ifelse(rownames(tmp) %in% paste0(l2b,'-01'),'L2l','L2s')
  gema_dif$L2l[i] = median(tmp$gene[tmp$group=='L2l'])
  gema_dif$L2s[i] = median(tmp$gene[tmp$group=='L2s'])
  gema_dif$pval_w[i] = wilcox.test(gene ~ group, data= tmp)$p.value
  gema_dif$pval_t[i] = t.test(gene ~ group, data= tmp)$p.value
}
gema_dif$log2FC = gema_dif$L2l-gema_dif$L2s
gema_dif$type = ifelse(gema_dif$pval_t<0.05 & gema_dif$log2FC>0.5,'up',
                       ifelse(gema_dif$pval_t<0.05 & gema_dif$log2FC< -0.5,'down','nosig'))
gema_dif$lab = ifelse(gema_dif$pval_t<0.05 & abs(gema_dif$log2FC)>1,gema_dif$gene,NA)
selectMarkerDf = gema_dif[gema_dif$log2FC > 0.5, ]
selectMarkerDf = selectMarkerDf[selectMarkerDf$pval_t < 0.05, ]

markerGeneList = selectMarkerDf$gene
gema_dif = gema_dif[order(gema_dif$log2FC, decreasing = T),]
fcrnk = gema_dif$log2FC
names(fcrnk) = gema_dif$gene  
fcrnk = fcrnk[!is.na(fcrnk)]
fcrnk = fcrnk[!duplicated(names(fcrnk))]
geneSetListArr = list('geneSet' = markerGeneList)
geneSetTablePwd = './geneSet240327.csv'
geneSetDf = data.frame(fread(geneSetTablePwd))
geneSetDf$Gene_set_name
geneSetDf = geneSetDf[c(2,7), ]
rm('resDf4geneSet')
for (tmpIndex in 1:nrow(geneSetDf)) {
  tmpGeneSetName = geneSetDf$Gene_set_name[tmpIndex]
  tmpGeneSet = strsplit(geneSetDf$core_enrichment[tmpIndex], split = "/")[[1]]
  
  tmpDf = data.frame(term = geneSetDf$Gene_set_name[tmpIndex],
                     gene = tmpGeneSet)
  
  if(!exists('resDf4geneSet')){
    resDf4geneSet = tmpDf
  }else{
    resDf4geneSet = rbind(resDf4geneSet, tmpDf)
  }
  
}

kpathways = resDf4geneSet
gsea = GSEA(fcrnk,  TERM2GENE  = kpathways,  eps = 1e-300)
g = gsea@result
g = g[order(g$NES<0,  g$p.adjust), ]

g
p <- gseaplot2(gsea, geneSetID = g$ID, 
               color = c('#219C90', '#EE9322'), 
               rel_heights = c(1.2, 0.2, 0.6), 
               ES_geom = 'line', pvalue_table = F,  base_size = 12)
p[[1]] = p[[1]] + guides(color = 'none') +
  theme(text = element_text(size = 15))
p

