
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
library(qusage)
library(stringr)

library(reshape2)
library(data.table)
library(ggthemes)
library(ggsci)
allMigrGeneSetpwd = './migration_genesets.v2023.2.Hs.gmt'

allMigrGeneSetArr = read.gmt(allMigrGeneSetpwd)

allMigrGeneSetName = names(allMigrGeneSetArr)
selectedGeneSetPwd = './migrGeneSet240330_selectedMigrOnly.csv'
selectedGeneSetDf = data.frame(fread(selectedGeneSetPwd))

resSelectedMigrGeneSetArr = c()
for (tmpIndex in 1:nrow(selectedGeneSetDf)) {
  tmpName = allMigrGeneSetName[grepl(selectedGeneSetDf$name[tmpIndex], allMigrGeneSetName)]
  
  resSelectedMigrGeneSetArr = c(resSelectedMigrGeneSetArr,
                                tmpName)
}

nrow(selectedGeneSetDf)
length(resSelectedMigrGeneSetArr)
migrGeneSetNameList = allMigrGeneSetArr[names(allMigrGeneSetArr) %in% resSelectedMigrGeneSetArr]
length(migrGeneSetNameList)

setFilterOutName = c('GOBP_NEGATIVE_REGULATION_OF_CELL_MIGRATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS', 
                     'GOBP_CELL_MIGRATION_INVOLVED_IN_GASTRULATION',
                     "GOBP_REGULATION_OF_FIBROBLAST_MIGRATION",
                     "WP_NEURAL_CREST_CELL_MIGRATION_IN_CANCER",
                     "GOBP_AMEBOIDAL_TYPE_CELL_MIGRATION"   
)
migrGeneSetNameList = migrGeneSetNameList[!(names(migrGeneSetNameList) %in% setFilterOutName)]
addGeneSetPwd = './migrationGenes_fromPM30353164.txt'

addGeneSetDf = read.delim(addGeneSetPwd, col.names = 'gene')

migrGeneSetNameList$migrationGenes_fromPM30353164 = addGeneSetDf$gene

length(migrGeneSetNameList)
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
gema0 = read.delim('./gbm_tcga_data_mrna_affymetrix_microarray.txt',check.names = F)

gema = gema0[, names(gema0) %in% c('Hugo_Symbol',paste0(c(l2b,l2s),'-01'))]
gema0 = read.delim('./gbm_tcga_data_mrna_affymetrix_microarray.txt', check.names = F)
gema = gema0[, names(gema0) %in% c('Hugo_Symbol',paste0(c(l2b,l2s),'-01'))]
gema_dif = data.frame(gene = gema$Hugo_Symbol, L2l = 0, L2s = 0, pval_w =1, pval_t =1)
library(stringr)
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

gema_dif = gema_dif[order(gema_dif$log2FC, decreasing = T),]
fcrnk = gema_dif$log2FC
names(fcrnk) = gema_dif$gene  
fcrnk = fcrnk[!is.na(fcrnk)]
fcrnk = fcrnk[!duplicated(names(fcrnk))]
library(qusage)

geneSetC2df = migrGeneSetNameList
library(parallel)
numCores = detectCores()
numCores

myCalculateGeneSet = function(tmpIndex){
  
  geneSetListArr = list(geneSet = geneSetC2df[[tmpIndex]])
  tmpResGsea = fgsea(pathways = geneSetListArr,
                     stats = fcrnk)
  
  tmpGeneSetName = names(geneSetC2df)[tmpIndex]
  
  tmpSize = tmpResGsea$size
  tmpPval = tmpResGsea$pval
  tmpNes = tmpResGsea$NES
  
  if(length(tmpResGsea$size) == 0){tmpSize = -1}
  if(length(tmpResGsea$pval) == 0){tmpPval = -1}
  if(length(tmpResGsea$NES) == 0){tmpNes = -1}
  
  tmpDf = data.frame(geneSetName = tmpGeneSetName,
                     geneNum = tmpSize,
                     pValue = tmpPval,
                     nes = tmpNes)
  
  return(tmpDf)
}
resGsea4c2_para = mclapply(1:length(geneSetC2df),
                           myCalculateGeneSet,
                           mc.cores = round(numCores*2/3))
resGsea4c2Df_para = as.data.frame(do.call(rbind, resGsea4c2_para))
dim(resGsea4c2Df_para)
add_line_breaks_prob <- function(x, max_interval = 40) {
  if(str_length(x) <= max_interval){
    return(x)
  }else{
    length_interval = ceiling(str_length(x)/2)
    x = paste0(substr(x, 1, length_interval), "\n", substr(x, length_interval+1, str_length(x)))
    return(x)
  }
  
}
resGsea4c2Df_para4plt = resGsea4c2Df_para
colnames(resGsea4c2Df_para4plt)
resGsea4c2Df_para4plt$geneSetName = unlist(lapply(resGsea4c2Df_para4plt$geneSetName, add_line_breaks_prob))

cat(tolower(as.character(resGsea4c2Df_para4plt$geneSetName)), sep = '", "')

resGsea4c2Df_para4plt$name4plt = c("GOBP_Cell_Migration", 
                                   "GOBP_Cerebral_Cortex_Radially_
Oriented_Cell_Migration", "GOBP_Positive_Regulation_of_Blood_
Vessel_Endothelial_Cell_Migration", "GOBP_Positive_Regulation_of_Cell_
Migration_Involved_in_Sprouting_Angiogenesis", "GOBP_Positive_Regulation_of_
Mononuclear_Cell_Migration", "HP_Abnormality_of_Neuronal_Migration", "WU_Cell_Migration", "Migration_Set_PM30353164")

resGsea4c2Df_para4plt$name4plt = gsub('_', ' ', resGsea4c2Df_para4plt$name4plt)

resGsea4c2Df_para4plt$name4plt = factor(resGsea4c2Df_para4plt$name4plt,
                                        levels = resGsea4c2Df_para4plt$name4plt[order(resGsea4c2Df_para4plt$nes)])

resGsea4c2Df_para4plt$pValue2 = -log10(resGsea4c2Df_para4plt$pValue)*sign(resGsea4c2Df_para4plt$nes)
ggplot(resGsea4c2Df_para4plt, mapping = aes(x = name4plt,
                                            y = nes,
                                            fill = nes,
                                            
)) +
  geom_bar(stat = 'identity', width = 0.6, just = 0.5, col = 'black') +
  coord_flip() +
  geom_hline(yintercept = 0) +
  theme_classic() +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0,
  ) +
  labs(x = '',
       y = 'Normalized Enrichment Score (NES)',
       fill = 'Enrichment score') +
  theme(text = element_text(size = 15),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
  ) 

