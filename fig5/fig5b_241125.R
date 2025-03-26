
{
  library(data.table)
  newIdhLabelPwd = '/Users/jiabao/Documents/1_github/mfgCodes/expression/mgbmSgbmLabelMerge230907.csv'
  newIdhLabelDf = data.frame(fread(newIdhLabelPwd))
  
  newIdhLabel4sgbmList = newIdhLabelDf$ID[newIdhLabelDf$idhMannual2 == 'WT' & newIdhLabelDf$bigSmall == 'SGBM']
  
  source('/Users/jiabao/OneDrive - HKUST Connect/0_workSpace/thesis/codes/chapter34/labMeeting220126/compareClinical/earlyVSlate/getHeatmap220518.R')
  
  
  library(DT)
  library("survminer")
  library(survival)
  
  myMax = function(x){
    return(max(as.numeric(x), na.rm = T))
  }
  
  
  customize_labels <- function (p, font.title = NULL,
                                font.subtitle = NULL, font.caption = NULL,
                                font.x = NULL, font.y = NULL, font.xtickslab = NULL, font.ytickslab = NULL)
  { # function is from tutorial
    original.p <- p
    if(is.ggplot(original.p)) list.plots <- list(original.p)
    else if(is.list(original.p)) list.plots <- original.p
    else stop("Can't handle an object of class ", class (original.p))
    .set_font <- function(font){
      font <- ggpubr:::.parse_font(font)
      ggtext::element_markdown (size = font$size, face = font$face, colour = font$color)
    }
    for(i in 1:length(list.plots)){
      p <- list.plots[[i]]
      if(is.ggplot(p)){
        if (!is.null(font.title)) p <- p + theme(plot.title = .set_font(font.title))
        if (!is.null(font.subtitle)) p <- p + theme(plot.subtitle = .set_font(font.subtitle))
        if (!is.null(font.caption)) p <- p + theme(plot.caption = .set_font(font.caption))
        if (!is.null(font.x)) p <- p + theme(axis.title.x = .set_font(font.x))
        if (!is.null(font.y)) p <- p + theme(axis.title.y = .set_font(font.y))
        if (!is.null(font.xtickslab)) p <- p + theme(axis.text.x = .set_font(font.xtickslab))
        if (!is.null(font.ytickslab)) p <- p + theme(axis.text.y = .set_font(font.ytickslab))
        list.plots[[i]] <- p
      }
    }
    if(is.ggplot(original.p)) list.plots[[1]]
    else list.plots
  }
  
  
  
  cat(patientNameListAfterOrder, sep = ',')
  
  
  
  
  
  
  
  snvPwd = '/Volumes/Expansion/download/1_data/gbm_tcga/data_mutations.txt'
  
  library(data.table)
  
  snvDf = data.frame(fread(snvPwd))
  
  snvDf$patientID = substr(snvDf$Tumor_Sample_Barcode, 1, 12)
  
  
  gbmLabelPwd = '/Users/jiabao/Documents/1_github/mfgCodes/expression/mgbmSgbmLabel.csv'
  gbmLabelDf = data.frame(fread(gbmLabelPwd))
  
  sgbmPatientList = gbmLabelDf$ID[gbmLabelDf$bigSmall == 'SGBM']
  length(sgbmPatientList)
  
  
  snvDf4sgbm = snvDf[snvDf$patientID %in% sgbmPatientList, ]
  
  sgbmPatientList4snv = unique(snvDf4sgbm$patientID)
  length(sgbmPatientList4snv)
  
  sgbmPatientList4idhMut = unique(snvDf4sgbm$patientID[snvDf4sgbm$Hugo_Symbol == 'IDH1'])
  
  length(sgbmPatientList4idhMut)
  
  
  sgbmPatientList4idhWt = sgbmPatientList4snv[!(sgbmPatientList4snv %in% sgbmPatientList4idhMut)]
  
  setdiff(sgbmPatientList4snv, sgbmPatientList4idhMut)
  
  
  
  library(data.table)
  mgbmLabelPwd = '/Users/jiabao/Documents/1_github/mfgCodes/expression/mgbmSgbmLabel.csv'
  mgbmLabelDf = data.frame(fread(mgbmLabelPwd))
  
  cell2016idhWtPwd = '/Users/jiabao/Documents/1_github/mfgCodes/expression/cell2016IdhWt.csv'
  
  cell2016idhWtDf = data.frame(fread(cell2016idhWtPwd))
  
  
  cell2016idhPwd = '/Users/jiabao/Documents/1_github/mfgCodes/expression/cell2016Idh.csv'
  cell2016idhDf = data.frame(fread(cell2016idhPwd))
  
  cell2016idhPatientList = cell2016idhWtDf$Case
  
  mgbmLabelDf$cell2016idhWt = NA
  mgbmLabelDf$cell2016idhWt[mgbmLabelDf$ID %in% cell2016idhPatientList] = 'WT'
  
  mgbmLabelDf$cell2016idh = cell2016idhDf$IDH.status[match(mgbmLabelDf$ID, cell2016idhDf$Case)]
  
  mgbmLabelDf = mgbmLabelDf[!is.na(mgbmLabelDf$modality), ]
  
  
  colnames(mgbmLabelDf)
  table(mgbmLabelDf[, c("bigSmall", "cell2016idhWt")], useNA = 'ifany')
  
  
  table(mgbmLabelDf$bigSmall)
  
  group1Patient = mgbmLabelDf$ID[mgbmLabelDf$bigSmall == 'Big']
  group1Label = 'Big'
  group2Patient = mgbmLabelDf$ID[mgbmLabelDf$bigSmall == 'Small']
  group2Label = 'Small'
  group3Patient = newIdhLabel4sgbmList
  group3Label = 'sgbmIDHwt'
  
  cat(group1Patient, sep = ',')
  cat(group2Patient, sep = ',')
  cat(group3Patient, sep = ',')
  
  length(unique(c(group1Patient, group2Patient)))
  length(unique(c(group1Patient, group2Patient, group3Patient)))
  
  
  newClinicalGbmAll = fread('/Users/jiabao/Library/CloudStorage/OneDrive-HKUSTConnect/0_workSpace/thesis/codes/chapter34/labMeeting220126/compareClinical/clinical.project-TCGA-GBM.2022-01-18/clinical.tsv')
  
  
  newClinicalTarget = newClinicalGbmAll[case_submitter_id %in% c(group1Patient, group2Patient, group3Patient)]
  groupGDC_newClinicalTarget = newClinicalTarget[,c('case_submitter_id','days_to_death','days_to_last_follow_up','vital_status')]
  
  dim(unique(groupGDC_newClinicalTarget))
  
  survivalTimeList = apply(groupGDC_newClinicalTarget[,c('days_to_death','days_to_last_follow_up')],1,myMax)
  survivalTimeList
  
  groupGDC_newClinicalTarget$time = survivalTimeList
  groupGDC_newClinicalTarget = unique(groupGDC_newClinicalTarget)
  
  colnames(groupGDC_newClinicalTarget)
  
  groupGDC_newClinicalTarget[groupGDC_newClinicalTarget$time < 0, ]
  groupGDC_newClinicalTarget = groupGDC_newClinicalTarget[groupGDC_newClinicalTarget$time >= 0, ]
  dim(groupGDC_newClinicalTarget)
  
  
  
  setdiff(c(group1Patient, group2Patient), groupGDC_newClinicalTarget$case_submitter_id)
  
  setdiff(c(group1Patient, group2Patient, group3Patient), groupGDC_newClinicalTarget$case_submitter_id)
  
  
  
  
  groupGDC_newClinicalTarget$group = 'unknown'
  groupGDC_newClinicalTarget$group[groupGDC_newClinicalTarget$case_submitter_id %in% group1Patient] = group1Label
  groupGDC_newClinicalTarget$group[groupGDC_newClinicalTarget$case_submitter_id %in% group2Patient] = group2Label
  groupGDC_newClinicalTarget$group[groupGDC_newClinicalTarget$case_submitter_id %in% group3Patient] = group3Label
  
  
}
{
  
  twoGroupDf = groupGDC_newClinicalTarget
  twoGroupDf = twoGroupDf[twoGroupDf$group == 'Big', ]
  
  twoGroupDf$state = 0
  twoGroupDf$state[twoGroupDf$vital_status=='Alive'] = 1
  twoGroupDf$state[twoGroupDf$vital_status=='Dead'] = 2
  twoGroupDf$group = twoGroupDf$group
  twoGroupDf$group = factor(twoGroupDf$group, 
                            levels = c(group1Label, group2Label, group3Label))
  
  twoGroupDf$time = twoGroupDf$time/30
  
  
  fit = survfit(Surv(time, state) ~ group, data = twoGroupDf)
  
  
  ggsurv = ggsurvplot(
    fit, 
    data = twoGroupDf, 
    size = 1,                 # change line size
    ylab = 'Probability of Overall Survival',
    xlab = 'Overall Survival (Days)',
    palette = c('#D63821', '#FFFFFF', '#2ca25f'),
    conf.int = F,          # Add confidence interval
    pval.coord = c(50, 0.1),
    risk.table = F,        # Add risk table
    risk.table.col = "strata",# Risk table color by groups
    legend.labs = 'Big', # remove first 7 character
    risk.table.height = 0.25, # Useful to change when you have multiple groups
    break.x.by = 10,          # Set breaks every 10 units
    xlim = c(0, 60),           # Set x-axis limits
    ggtheme = theme_classic() + theme(axis.text.x = element_text(size = 12),
                                      axis.text.y = element_text(size = 12)) 
    
    
  )
  
  ggsurv = customize_labels(
    ggsurv,
    font.title    = c(16, "bold", "black"),
    font.subtitle = c(15, "bold.italic", "black"),
    font.caption  = c(14, "plain", "black"),
    font.x        = c(14, "bold.italic", "black"),
    font.y        = c(14, "bold.italic", "black"),
    font.xtickslab = c(12, "plain", "black")
  ) 
  
  print(ggsurv)
  
  
  
}
{
  
  twoGroupDf = groupGDC_newClinicalTarget
  twoGroupDf = twoGroupDf[twoGroupDf$group == 'Small', ]
  
  twoGroupDf$state = 0
  twoGroupDf$state[twoGroupDf$vital_status=='Alive'] = 1
  twoGroupDf$state[twoGroupDf$vital_status=='Dead'] = 2
  twoGroupDf$group = twoGroupDf$group
  twoGroupDf$group = factor(twoGroupDf$group, 
                            levels = c(group1Label, group2Label, group3Label))
  
  twoGroupDf$time = twoGroupDf$time/30
  
  
  fit = survfit(Surv(time, state) ~ group, data = twoGroupDf)
  
  
  ggsurv = ggsurvplot(
    fit, 
    data = twoGroupDf, 
    size = 1,                 # change line size
    ylab = 'Probability of Overall Survival',
    xlab = 'Overall Survival (Days)',
    palette = c('#224B89', '#2ca25f'),
    conf.int = F,          # Add confidence interval
    pval.coord = c(50, 0.1),
    risk.table = F,        # Add risk table
    risk.table.col = "strata",# Risk table color by groups
    legend.labs = 'Small', # remove first 7 character
    risk.table.height = 0.25, # Useful to change when you have multiple groups
    break.x.by = 10,          # Set breaks every 10 units
    xlim = c(0, 60), 
    ggtheme = theme_classic() + theme(axis.text.x = element_text(size = 12),
                                      axis.text.y = element_text(size = 12)) 
    
  )
  
  ggsurv = customize_labels(
    ggsurv,
    font.title    = c(16, "bold", "black"),
    font.subtitle = c(15, "bold.italic", "black"),
    font.caption  = c(14, "plain", "black"),
    font.x        = c(14, "bold.italic", "black"),
    font.y        = c(14, "bold.italic", "black"),
    font.xtickslab = c(12, "plain", "black")
  )
  
  print(ggsurv)
  
  
}
{
  
  
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
  
  library(reshape2)
  library(data.table)
  library(ggthemes)
  library(ggsci)
  
  library(DT)
  library(survminer)
  library(survival)
  
  
  
  

  source(file.path('/Users/jiabao/Documents/1_github/mfgCodes/code4figs24112001/fig5/plotfgsea_241125.R'))
  
  
  
  
  
  idh = read.delim(file.path(setFilePwd, 'TCGA_GBM_MRI259_IDH.txt'))
  mdp = read.delim(file.path(setFilePwd, 'gbm_tcga_data_clinical_patient.txt'), comment.char = '#', na.strings = c('','[Not Available]'))
  mds = read.delim(file.path(setFilePwd, 'gbm_tcga_data_clinical_sample.txt'), comment.char = '#')
  fg = read.delim(file.path(setFilePwd, 'gbm_tcga_data_clinical_patient.txt'),comment.char = '#')
  
  cl = read.csv(file.path(setFilePwd, 'mgbmSgbmLabel.csv'))
  cl = cl[which(cl$modality=='MR'),]
  cl$IDH = idh$IDH.status[match(cl$ID, idh$Patient)]
  
  l2b = cl$ID[cl$bigSmall=='Big']
  l2s = cl$ID[cl$bigSmall=='Small']
  l1 = cl$ID[cl$bigSmall=='SGBM'&cl$IDH=='WT']
  
  
  
  
  
  
  realTumorMelt = data.frame(fread('/Users/jiabao/Documents/1_github/mfgCodes/code4figs24010802/fig1/240105_volumeByLesion.csv'))
  colnames(realTumorMelt) = c('patient', 'lesion', 'realVolume')
  
  
  setFillColor = rev(c('#9B4C4B', '#BC8D78', '#D6D2B0',
                       '#96C0EE', '#C1D8F3', '#485C81'))
  
  totalVolume = aggregate(realVolume ~ patient, data = realTumorMelt, sum)
  orderedPatients = totalVolume$patient[order(totalVolume$realVolume)]
  realTumorMelt$patient = factor(realTumorMelt$patient, levels = orderedPatients)
  
  
  
  
  
  
  
  
  
  
  
  myMax = function(x){
    return(max(x, na.rm = T))
  }
  
  
  newClinicalGbmAll = fread('/Users/jiabao/Documents/1_github/mfgCodes/code4figs24010802/fig1/clinical.tsv')
  
  newClinicalMgbm = newClinicalGbmAll
  
  groupGDC_clinicalNew = newClinicalMgbm[,c('case_submitter_id','days_to_death','days_to_last_follow_up','vital_status')]
  mgbm_t3New = apply(groupGDC_clinicalNew[,c('days_to_death','days_to_last_follow_up')],1,myMax)
  mgbmDfNew = data.frame('patient' = groupGDC_clinicalNew$case_submitter_id,
                         "time" = mgbm_t3New, 
                         "status" = groupGDC_clinicalNew$vital_status)
  survivalDf = unique(mgbmDfNew)
  survivalDf$time = as.numeric(survivalDf$time)
  
  
  
  sgbmMgbmLabelPwd = '/Users/jiabao/Documents/1_github/mfgCodes/therapySurvival/mgbmSgbmLabel.csv'
  sgbmMgbmLabelDf = data.frame(fread(sgbmMgbmLabelPwd))
  sgbmMgbmLabelDf = na.omit(sgbmMgbmLabelDf) # filter out the patients without MRI
  
  dim(sgbmMgbmLabelDf)
  
  
  
  
  
  
  gema0 = read.delim(file.path(setFilePwd, 'gbm_tcga_data_mrna_affymetrix_microarray.txt'),check.names = F)
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
  
  markerGeneList = selectMarkerDf$gene
  
  
  
  
  gema_dif = gema_dif[order(gema_dif$log2FC, decreasing = T),]
  fcrnk = gema_dif$log2FC
  names(fcrnk) = gema_dif$gene  
  fcrnk = fcrnk[!is.na(fcrnk)]
  fcrnk = fcrnk[!duplicated(names(fcrnk))]
  
  
  
  
  
  
  
  geneSetTablePwd = '/Users/jiabao/Documents/1_github/mfgCodes/predictMgbmByExpression24032502/geneSet240327.csv'
  
  
  geneSetDf = data.frame(fread(geneSetTablePwd))
  
  dim(geneSetDf)
  
  
  
  set.seed(2024)
  tmpIndex = 4
  
  source(file.path('/Users/jiabao/Documents/1_github/mfgCodes/code4figs24112001/fig5/plotfgsea_241125.R'))
  
  {
    tmpGeneSetName = geneSetDf$Gene_set_name[tmpIndex]
    tmpGeneSet = strsplit(geneSetDf$core_enrichment[tmpIndex], split = "/")[[1]]
    
    geneSetListArr = list(geneSet = tmpGeneSet)
    plotfgsea(pathway = geneSetListArr[['geneSet']],
              stats = fcrnk, 
              fgseaRes = fgsea(pathways = geneSetListArr,
                               stats = fcrnk),
              gene.set.name = 'geneSet',
              gene.set.title.name = tmpGeneSetName, posClass = 'L2-large',
              negClass = 'L2-small', setCol4curve = 'red', setCex = 1.2)
  }
  
  
  
}
{
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
  
  
  
  
  
  
  
  
  allMigrGeneSetpwd = '/Users/jiabao/Documents/1_github/mfgCodes/predictMgbmByExpression24032502/migration_genesets.v2023.2.Hs.gmt'
  
  allMigrGeneSetArr = read.gmt(allMigrGeneSetpwd)
  
  allMigrGeneSetName = names(allMigrGeneSetArr)
  
  
  
  selectedGeneSetPwd = '/Users/jiabao/Documents/1_github/mfgCodes/predictMgbmByExpression24032502/migrGeneSet240330_selectedMigrOnly.csv'
  
  
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
  
  
  addGeneSetPwd = '/Users/jiabao/Library/CloudStorage/Dropbox/Jiabao/mfgManuscript221028/Fig4/TCGA/migrationGenes_fromPM30353164.txt'
  
  addGeneSetDf = read.delim(addGeneSetPwd, col.names = 'gene')
  
  migrGeneSetNameList$migrationGenes_fromPM30353164 = addGeneSetDf$gene
  
  length(migrGeneSetNameList)
  
  
  
  
  
  setFilePwd = '/Users/jiabao/Documents/1_github/mfgCodes/predictMgbmByExpression24032502/TCGA'
  
  idh = read.delim(file.path(setFilePwd, 'TCGA_GBM_MRI259_IDH.txt'))
  mdp = read.delim(file.path(setFilePwd, 'gbm_tcga_data_clinical_patient.txt'), comment.char = '#', na.strings = c('','[Not Available]'))
  mds = read.delim(file.path(setFilePwd, 'gbm_tcga_data_clinical_sample.txt'), comment.char = '#')
  fg = read.delim(file.path(setFilePwd, 'gbm_tcga_data_clinical_patient.txt'),comment.char = '#')
  
  cl = read.csv(file.path(setFilePwd, 'mgbmSgbmLabel.csv'))
  cl = cl[which(cl$modality=='MR'),]
  cl$IDH = idh$IDH.status[match(cl$ID, idh$Patient)]
  
  l2b = cl$ID[cl$bigSmall=='Big']
  l2s = cl$ID[cl$bigSmall=='Small']
  l1 = cl$ID[cl$bigSmall=='SGBM'&cl$IDH=='WT']
  
  
  
  
  
  
  
  
  
  
  gema0 = read.delim(file.path(setFilePwd, 'gbm_tcga_data_mrna_affymetrix_microarray.txt'),check.names = F)
  
  gema = gema0[, names(gema0) %in% c('Hugo_Symbol',paste0(c(l2b,l2s),'-01'))]
  
  
  
  
  
  
  
  
  library(matrixStats)
  library(circlize)
  library(ComplexHeatmap)
  
  ssgsea = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
    row_names = rownames(X)
    num_genes = nrow(X)
    gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})
    
    R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')
    
    es = apply(R, 2, function(R_col) {
      gene_ranks = order(R_col, decreasing = TRUE)
      
      es_sample = sapply(gene_sets, function(gene_set_idx) {
        indicator_pos = gene_ranks %in% gene_set_idx
        indicator_neg = !indicator_pos
        
        rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha
        
        step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
        step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)
        
        step_cdf_diff = step_cdf_pos - step_cdf_neg
        
        if (scale) step_cdf_diff = step_cdf_diff / num_genes
        
        if (single) {
          sum(step_cdf_diff)
        } else {
          step_cdf_diff[which.max(abs(step_cdf_diff))]
        }
      })
      unlist(es_sample)
    })
    
    if (length(gene_sets) == 1) es = matrix(es, nrow = 1)
    
    if (norm) es = es / diff(range(es))
    
    rownames(es) = names(gene_sets)
    colnames(es) = colnames(X)
    return(es)
  }
  
  
  data = gema
  data[1:5, 1:5]
  
  data = data[!duplicated(data$Hugo_Symbol), ]
  
  rownames(data) = data$Hugo_Symbol
  data$Hugo_Symbol = NULL
  
  colnames(data) = str_sub(colnames(data), 1, 12)
  
  data[1:5,1:5]
  dim(data)
  
  data = as.matrix(data)
  
  
  gene_sets = migrGeneSetNameList
  
  system.time(assign('res', ssgsea(data, gene_sets, scale = TRUE, norm = FALSE)))
  
  
  
  res1 = t(res)
  head(res1)
  
  
  mat = (res - rowMeans(res))/(rowSds(as.matrix(res)))[row(res)]
  
  mat = mat[, order(colnames(mat))]
  
  
  sgbmMgbmLabelPwd = '/Users/jiabao/Documents/1_github/mfgCodes/therapySurvival/mgbmSgbmLabel.csv'
  sgbmMgbmLabelDf = data.frame(fread(sgbmMgbmLabelPwd))
  sgbmMgbmLabelDf = na.omit(sgbmMgbmLabelDf) # filter out the patients without MRI
  
  
  mgbmLabelDf = sgbmMgbmLabelDf[sgbmMgbmLabelDf$bigSmall != 'SGBM', ]
  rownames(mgbmLabelDf) = mgbmLabelDf$ID
  
  mgbmLabelDf = mgbmLabelDf[rownames(mgbmLabelDf) %in% colnames(data) , ]
  
  
  bestSimPwd = '/Users/jiabao/Documents/1_github/mfgCodes/bestSim4eachPatient/resAllBestSimNew35rerun230822.csv'
  
  bestSimDf = data.frame(fread(bestSimPwd))
  
  bestSimDf$group = mgbmLabelDf$bigSmall[match(bestSimDf$p.ID,
                                               mgbmLabelDf$ID)]
  
  bestSimDf = na.omit(bestSimDf) # 34 patients with rna
  
  
  newLabelByl2l1ratioPwd = '/Users/jiabao/Documents/1_github/mfgCodes/code4figs24010802/fig1/mgbmL2LargeSmallLabelByL2L1ratio.csv'
  newLabelByl2l1ratioDf = data.frame(fread(newLabelByl2l1ratioPwd))
  
  
  multicentricPwd = '/Users/jiabao/Documents/1_github/mfgCodes/predictMgbmByExpression24032502/zhangWei2015_multicentricLabel30tcga.csv'
  multicentricDf = data.frame(fread(multicentricPwd))
  
  
  
  
  
  df = mgbmLabelDf[, c('bigSmall', 'highLow')]
  df$l2l1ratioLoS = newLabelByl2l1ratioDf$group[match(rownames(df),
                                                      newLabelByl2l1ratioDf$patient)]
  
  df$mriClass = multicentricDf$mri.classification[match(rownames(df),
                                                        multicentricDf$patientId)]
  df$mriClass[is.na(df$mriClass)] = "NA"
  
  df = df[, c('bigSmall', 'l2l1ratioLoS', 'mriClass', 'highLow')]
  ha = HeatmapAnnotation(df = df,
                         col = list(bigSmall = c("Big" = '#e41a1c', "Small" = '#377eb8'), 
                                    l2l1ratioLoS = c('Large' = '#dd1c77', 'Small' = '#7fcdbb'),
                                    mriClass = c('multifocal' = '#edf8b1', 'multicentric' = '#67000d', 'NA' = 'white'),
                                    highLow = c("High" = "orange", "Low" = "lightblue")
                                    
                         ), 
                         gap = unit(1, "mm"))
  
  
  
  mat2 = mat
  rownames(mat2)
  
  add_line_breaks <- function(x, length_interval) {
    x = paste0(substr(x, 1, length_interval), "\n", substr(x, length_interval+1, str_length(x)))
    return(x)
  }
  
  
  add_line_breaks_prob <- function(x, max_interval = 40) {
    if(str_length(x) <= max_interval){
      return(x)
    }else{
      length_interval = ceiling(str_length(x)/2)
      x = paste0(substr(x, 1, length_interval), "\n", substr(x, length_interval+1, str_length(x)))
      return(x)
    }
    
  }
  
  
  
  
  rownames(mat2) =unlist(lapply(rownames(mat2), add_line_breaks_prob))
  
  
  rownames(mat2)
  
  ht = Heatmap(mat2, top_annotation = ha, 
               col = colorRamp2(c(-2,0,2),
                                c("blue", "white", "red")),
               clustering_method_columns = "ward.D",
               
               cluster_rows = T, 
               show_row_dend = T, 
               show_column_dend = T,
               rect_gp = grid::gpar(col = "white"), show_heatmap_legend = F,
               
  )
  
  
  
  matDf = data.frame(mat)
  matDf$pathwayName = rownames(matDf)
  
  
  matDfMelt = reshape2::melt(matDf, values = 'pathwayName')
  colnames(matDfMelt) = c('pathwayName', 'patientId', 'ssgseaScore')
  matDfMelt$patientId = gsub('[.]', '-', matDfMelt$patientId)
  
  matDfMelt$group = mgbmLabelDf$bigSmall[match(matDfMelt$patientId, 
                                               mgbmLabelDf$ID)]
  
  
  
  newLabelByl2l1ratioDf 
  
  
  matDfMelt$group = newLabelByl2l1ratioDf$group[match(matDfMelt$patientId, 
                                                      newLabelByl2l1ratioDf$patient)]
  
  patientOrderList = mgbmLabelDf$ID[order(mgbmLabelDf$bigSmall)]
  
  
  tmpRes = hclust(dist(t(matDf[, 1:34])))
  resCutreeList = cutree(tmpRes, k = 2)
  
  cat(gsub('[.]', '-', names(resCutreeList)[resCutreeList == 1]), sep = "','")
  
  cat(gsub('[.]', '-', names(resCutreeList)[resCutreeList == 2]), sep = "','")
  
  
  matDfMelt$pathwayName = factor(matDfMelt$pathwayName,
                                 levels = tmpRes$labels[tmpRes$order])
  
  matDfMelt$patientId = factor(matDfMelt$patientId,
                               levels = patientOrderList)
  
  
  
  gema0 = read.delim(file.path(setFilePwd, 'gbm_tcga_data_mrna_affymetrix_microarray.txt'),check.names = F)
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
  
  gema_dif$log2FC = gema_dif$L2s-gema_dif$L2l
  
  
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
  
  
  
  
  names(geneSetC2df)
  
  tmpIndex = 6 # 2 6
  
  set.seed(2024)
  {
    tmpGeneSetName = as.character(resGsea4c2Df_para4plt$name4plt)[tmpIndex] # names(geneSetC2df)[tmpIndex]
    tmpGeneSet = unlist(geneSetC2df[tmpIndex])
    
    geneSetListArr = list(geneSet = tmpGeneSet)
    plotfgsea(pathway = geneSetListArr[['geneSet']],
              stats = fcrnk, 
              fgseaRes = fgsea(pathways = geneSetListArr,
                               stats = fcrnk),
              gene.set.name = 'geneSet',
              gene.set.title.name = tmpGeneSetName,
              negClass = 'L2-large', posClass = 'L2-small', setCol4curve = '#2D4B83', setCex = 1.2)
  }
  
  
  
  
}
somaticMutDf = data.frame(TP53 = c(4/13, 5/7),
                          EGFR = c(5/13, 0/7),
                          EGFRamp = c(13/22, 4/13),
                          PI3Ks = c(9/13, 2/7), # PI3K pathway
                          group = c('Large', 'Small'))
somaticMutDf = data.frame(t(somaticMutDf))
colnames(somaticMutDf) = c('Large', 'Small')
somaticMutDf$id = rownames(somaticMutDf)

somaticMutDf = somaticMutDf[1:4, ]
somaticMutDfMelt = reshape2::melt(somaticMutDf)
somaticMutDfMelt$Large = as.numeric(somaticMutDfMelt$Large)
somaticMutDfMelt$Small = as.numeric(somaticMutDfMelt$Small)

ggplot(somaticMutDfMelt, mapping = aes(x = id,
                                       y = Large,
                                       label = signif(Large, digits = 2))) +
  geom_bar(stat = 'identity', width = 0.9, fill = 'red', alpha=0.7) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  geom_text(hjust = 0) +
  theme_classic() +
  coord_flip() +
  theme(text = element_text(size = 15), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(face = c("bold.italic"))) +
  labs(x = '')
ggplot(somaticMutDfMelt, mapping = aes(x = id,
                                       y = Small,
                                       label = signif(Small, digits = 2))) +
  geom_bar(stat = 'identity', width = 0.9, fill = '#2D4B83', alpha=0.7) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  geom_text(hjust = 0) +
  theme_classic() +
  coord_flip() +
  theme(text = element_text(size = 15), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(face = c("bold.italic"))) +
  labs(x = '')
library(data.table)
library(progress)
library(abc)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(ggsci)
library(scico)
pltFigureCount = 0

setTolerance = 0.1

{
  
  
  bigSmallLabelPwd = '/Users/jiabao/Documents/1_github/mfgCodes/expression/mgbmSgbmLabel.csv'
  bigSmallDf = data.frame(fread(bigSmallLabelPwd))
  bigSmallDf = bigSmallDf[bigSmallDf$bigSmall != 'SGBM',]
  
  
  {
    
    myListToCol = function(x){
      ColDf = data.frame(col = table(x))
      colnames(ColDf) = c('Value', 'Freq')
      ColDf$Col = rainbow(dim(ColDf)[1])
      return(ColDf$Col[match(x, ColDf$Value)])
    }
    
    myListToLog10 = function(x){
      return(ifelse(x==0, 0, log10(x)))
    }
    
    nColOfGgsci = 7
    p = ggplot(data.frame(x = 1:nColOfGgsci,
                          y = 1:nColOfGgsci),
               mapping = aes(x = x,
                             y = y,
                             col =  as.factor(x))) +
      geom_point() +
      scale_color_startrek()
    
    g = ggplot_build(p)
    sci7ColArr = unique(g$data[[1]]["colour"])[[1]]
    sci7ColArrOrder = sci7ColArr[c(6,5,4,3,2,1,7)]
    
    pwd = '/Users/jiabao/Documents/1_github/mfgCodes/plotAbcHeatmap'
    setwd(pwd)
    
    load(file = 'readSimV134.RData')
    dim(resSimAll)
    
    colnames(resSimAll)
    
    resSimAll$migr = resSimAll$migr/1e3 # in the nich level
    
    
    realSizeDf = data.frame(fread('realTCGA211017.csv'))
    dim(realSizeDf)
    realSizeDf[is.na(realSizeDf)] = 'Unknown'
    dim(realSizeDf)
    
    resSimAll
    
    obsDat = realSizeDf[,c('Num', paste0('T',1:6))]
    patientName = realSizeDf$ID
    obsDat[1,]
    
    simDat = resSimAll[resSimAll$num>=1,]
    simDat[is.na(simDat)] = 0
    
    
    myMad = function(x){
      return(median(1/qnorm(3/4)*abs(x-median(x))))
    }
    
    myAbcPr = function(setTarget, setSimulation, setLabel, setTol){
      target = setTarget
      sim = setSimulation
      modelLabel1 = as.factor(setLabel)
      tol = setTol
      resMad = c()
      for (i in 1:ncol(sim)) {
        madF = myMad(as.numeric(sim[,i]))
        madF = ifelse(madF==0, 1, madF)
        resMad = c(resMad, madF)
      }
      
      target.scale = target
      sim.scale = sim
      for (i in 1:ncol(sim)) {
        target.scale[,i] = target[,i] / resMad[i]
        sim.scale[,i] = sim[,i] / resMad[i]
      }
      
      resDist = 0
      for (i in 1:ncol(sim.scale)) {
        resDist.tmp = (sim.scale[,i] - target.scale[,i]) ^2
        resDist = resDist + resDist.tmp
      }
      ntol = ceiling(nrow(sim.scale) * tol)
      cutOffN = sort(resDist)[ntol]
      postRindex = resDist <= cutOffN
      resTmpDf = aggregate(rank(-resDist[postRindex]), 
                           by=list(label = modelLabel1[postRindex]), FUN=sum)
      resTmpDf = resTmpDf[order(resTmpDf$label),]
      
      tabTmp = table(modelLabel1[postRindex])
      
      for (i in 1:dim(resTmpDf)[1]) {
        tabTmp[names(tabTmp) == resTmpDf$label[i]] = resTmpDf$x[i]
      }
      return(prop.table(tabTmp))
    }
    
    myAbcPrBestSim = function(setTarget, setSimulation, setLabel, setTol){
      target = setTarget
      sim = setSimulation
      modelLabel1 = as.factor(setLabel)
      tol = setTol
      resMad = c()
      for (i in 1:ncol(sim)) {
        madF = myMad(as.numeric(sim[,i]))
        madF = ifelse(madF==0, 1, madF)
        resMad = c(resMad, madF)
      }
      
      target.scale = target
      sim.scale = sim
      for (i in 1:ncol(sim)) {
        target.scale[,i] = target[,i] / resMad[i]
        sim.scale[,i] = sim[,i] / resMad[i]
      }
      
      resDist = 0
      for (i in 1:ncol(sim.scale)) {
        resDist.tmp = (sim.scale[,i] - target.scale[,i]) ^2
        resDist = resDist + resDist.tmp
      }
      ntol = ceiling(nrow(sim.scale) * tol)
      cutOffN = sort(resDist)[ntol]
      postRindex = resDist <= cutOffN
      return(postRindex)
    }
    
    
    myResMad = function(setTarget, setSimulation){
      target = setTarget
      sim = setSimulation
      resMad = c()
      for (i in 1:ncol(sim)) {
        madF = myMad(as.numeric(sim[,i]))
        madF = ifelse(madF==0, 1, madF)
        resMad = c(resMad, madF)
      }
      return(resMad)
    }
    
    simDatOnly = simDat
    simDatOnly[, paste0('P',1:10)] = 100 * simDatOnly[, paste0('P',1:10)]
    
    unique(simDatOnly$migr)
    unique(simDatOnly$select)
    unique(simDatOnly$distMigrDrvNum)
    
    
    pb = progress_bar$new(total = length(patientName), clear = FALSE)
    rm(resPostAllBest)
    
    for (i in 1:length(patientName)) {
      pb$tick()
      resPost = myAbcPrBestSim(setTarget = obsDat[i,1:7], 
                               setSimulation = simDatOnly[,c('num', paste0('P',1:6))], 
                               setLabel = simDatOnly$migr, 
                               setTol =setTolerance) 
      
      tmpPostRes = simDatOnly[resPost,]
      tmpPostRes$patientId = patientName[i]
      
      
      if(i == 1){
        resPostAllBest = tmpPostRes
      }else{
        resPostAllBest = rbind(resPostAllBest, tmpPostRes)
      }
    }
    
    resPostAllBest
    
    
    resPostAllBest$group = bigSmallDf$bigSmall[match(resPostAllBest$patientId,
                                                     bigSmallDf$ID)]
  }
  
  
  library(ggpubr)
  my_comparison = list(c(1, 2))
  
  colnames(resPostAllBest)
  
  
  resMean = c()
  for (i in 1:length(patientName)) {
    tmpRes = mean(resPostAllBest$firstMigrNum[resPostAllBest$patientId == patientName[i]])
    resMean = c(resMean, tmpRes)
  }
  
  resMeanDf4bigSmall = data.frame(firstMigrNumMean = resMean,
                                  patientId = patientName)
  resMeanDf4bigSmall$group = bigSmallDf$bigSmall[match(resMeanDf4bigSmall$patientId,
                                                       bigSmallDf$ID)]
  
  
  
  
  
}
library(ggbeeswarm)
resMean4migr = c()
for (i in 1:length(patientName)) {
  tmpRes = mean(resPostAllBest$migr[resPostAllBest$patientId == patientName[i]])
  resMean4migr = c(resMean4migr, tmpRes)
}

resMeanDf4bigSmall4migr = data.frame(migrMean = resMean4migr,
                                     patientId = patientName)
resMeanDf4bigSmall4migr$group = bigSmallDf$bigSmall[match(resMeanDf4bigSmall4migr$patientId,
                                                          bigSmallDf$ID)]
resMean4invisibleLesionNum = c()
for (i in 1:length(patientName)) {
  tmpRes = mean(resPostAllBest$lesionNum[resPostAllBest$patientId == patientName[i]] - 
                  resPostAllBest$num[resPostAllBest$patientId == patientName[i]]
                
  )
  resMean4invisibleLesionNum = c(resMean4invisibleLesionNum, tmpRes)
}

resMeanDf4bigSmall4invisible = data.frame(invisibleMean = resMean4invisibleLesionNum,
                                          patientId = patientName)
resMeanDf4bigSmall4invisible$group = bigSmallDf$bigSmall[match(resMeanDf4bigSmall4invisible$patientId,
                                                               bigSmallDf$ID)]
resMeanDf4bigSmall$obsName = 'firstMigrMean'
resMeanDf4bigSmall4migr$obsName = 'meanMigr'
resMeanDf4bigSmall4invisible$obsName = 'invisibleNum'
library(ggridges)
ggplot(resMeanDf4bigSmall[resMeanDf4bigSmall$group == 'Big', ], mapping = aes(
  x = firstMigrNumMean,
  col = group)) +
  geom_density(fill = 'red', alpha = 0.5) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 15),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )+
  scale_color_brewer(palette = 'Set1') +
  guides(col = 'none') +
  scale_x_continuous(limits = c(1e4, 3e5)) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
  labs(x = 'Cell number at first migration',
       y = '')
ggplot(resMeanDf4bigSmall[resMeanDf4bigSmall$group == 'Small', ], mapping = aes(
  x = firstMigrNumMean,
  col = group)) +
  geom_density(fill = '#2D4B83', color = '#2D4B83', alpha = 0.5) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 15),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )+
  scale_color_brewer(palette = 'Set1') +
  guides(col = 'none') +
  scale_x_continuous(limits = c(1e4, 3e5)) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
  labs(x = 'Cell number at first migration',
       y = '')
ggplot(resMeanDf4bigSmall4invisible[resMeanDf4bigSmall4invisible$group == 'Big', ], mapping = aes(
  x = invisibleMean,
  col = group)) +
  geom_density(fill = 'red', alpha = 0.5) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 15),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )+
  scale_color_brewer(palette = 'Set1') +
  guides(col = 'none') +
  scale_x_continuous(limits = c(0.2, 35)) +
  
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
  labs(x = 'Mean invisible lesion number',
       y = '')
ggplot(resMeanDf4bigSmall4invisible[resMeanDf4bigSmall4invisible$group == 'Small', ], mapping = aes(
  x = invisibleMean,
  col = group)) +
  geom_density(fill = '#2D4B83', color = '#2D4B83', alpha = 0.5) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 15),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )+
  scale_color_brewer(palette = 'Set1') +
  guides(col = 'none') +
  scale_x_continuous(limits = c(0.2, 35)) +
  
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
  labs(x = 'Mean invisible lesion number',
       y = '')

