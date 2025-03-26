
library(data.table)
library(DT)
library("survminer")
library(survival)
newIdhLabelPwd = './mgbmSgbmLabelMerge230907.csv'
newIdhLabelDf = data.frame(fread(newIdhLabelPwd))

newIdhLabel4sgbmList = newIdhLabelDf$ID[newIdhLabelDf$idhMannual2 == 'WT' & newIdhLabelDf$bigSmall == 'SGBM']
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
library(data.table)
mgbmLabelPwd = './mgbmSgbmLabelWithCell2016idhWt.csv'
mgbmLabelDf = data.frame(fread(mgbmLabelPwd))
mgbmLabelDf = mgbmLabelDf[!is.na(mgbmLabelDf$modality), ]
mgbmLesionNumPwd = './realVolume_240412.csv'
mgbmLesionNumDf = data.frame(fread(mgbmLesionNumPwd))

mgbmLesionNumDf$Patient = mgbmLesionNumDf$patient
mgbmL2patientList = mgbmLesionNumDf$Patient[mgbmLesionNumDf$nTumor == 2]
mgbmL3456patientList = mgbmLesionNumDf$Patient[mgbmLesionNumDf$nTumor > 2]
group1Patient = mgbmLabelDf$ID[mgbmLabelDf$ID %in% mgbmL2patientList]
group1Label = 'L2'
group2Patient = mgbmLabelDf$ID[mgbmLabelDf$ID %in% mgbmL3456patientList]
group2Label = 'L3456'

newClinicalGbmAll = fread('./clinical.project-TCGA-GBM.2022-01-18/clinical.tsv')
newClinicalTarget = newClinicalGbmAll[case_submitter_id %in% c(group1Patient, group2Patient)]
groupGDC_newClinicalTarget = newClinicalTarget[,c('case_submitter_id','days_to_death','days_to_last_follow_up','vital_status')]
survivalTimeList = apply(groupGDC_newClinicalTarget[,c('days_to_death','days_to_last_follow_up')],1,myMax)
survivalTimeList

groupGDC_newClinicalTarget$time = survivalTimeList
groupGDC_newClinicalTarget = unique(groupGDC_newClinicalTarget)
groupGDC_newClinicalTarget[groupGDC_newClinicalTarget$time < 0, ]
groupGDC_newClinicalTarget = groupGDC_newClinicalTarget[groupGDC_newClinicalTarget$time >= 0, ]
setdiff(c(group1Patient, group2Patient), groupGDC_newClinicalTarget$case_submitter_id)

groupGDC_newClinicalTarget$group = 'unknown'
groupGDC_newClinicalTarget$group[groupGDC_newClinicalTarget$case_submitter_id %in% group1Patient] = group1Label
groupGDC_newClinicalTarget$group[groupGDC_newClinicalTarget$case_submitter_id %in% group2Patient] = group2Label

twoGroupDf = groupGDC_newClinicalTarget
twoGroupDf$state = 0
twoGroupDf$state[twoGroupDf$vital_status=='Alive'] = 1
twoGroupDf$state[twoGroupDf$vital_status=='Dead'] = 2
twoGroupDf$group = twoGroupDf$group
twoGroupDf$group = factor(twoGroupDf$group, 
                          levels = c(group1Label, group2Label))

twoGroupDf$time = twoGroupDf$time/30
fit = survfit(Surv(time, state) ~ group, data = twoGroupDf)

diff = survdiff(Surv(time, state) ~ group, data = twoGroupDf)

pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE) # get the t
log.rank <- survdiff(Surv(time, state) ~ group, data = twoGroupDf, rho = 0)
mantle.cox <- survdiff(Surv(time, state) ~ group, data = twoGroupDf, rho = 1)

model <- summary(coxph(Surv(time, state) ~ group, data = twoGroupDf))
HR <- round(model$conf.int[1],2)
HR.lower <- round(model$conf.int[3],2)
HR.upper <- round(model$conf.int[4],2)
log.rank.p <- round(1 - pchisq(log.rank$chi, df = 1), 4)
mantle.cox.p <- round(1 - pchisq(mantle.cox$chi, df = 1), 4)
c(HR, HR.lower, HR.upper, log.rank.p, mantle.cox.p)
ggsurv = ggsurvplot(
  fit, 
  data = twoGroupDf, 
  size = 1,                 # change line size
  ylab = 'Probability of Overall Survival',
  xlab = 'Overall Survival (Months)',
  palette = c('#102C57', '#E68369'),
  conf.int = F,          # Add confidence interval
  pval = F,              # Add p-value
  pval.coord = c(25, 0.8),
  risk.table = F,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs = substring(names(fit$strata),7), # remove first 7 character
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_classic() + theme(axis.text.x = element_text(size = 12),
                                    axis.text.y = element_text(size = 12))     # Change ggplot2 theme
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

ggsurvplot(
  fit, 
  data = twoGroupDf, 
  size = 0.5,                 # change line size
  ylab = 'Probability of Overall Survival',
  xlab = 'Overall Survival (Months)',
  palette = c('#102C57', '#E68369'),
  conf.int = F,          # Add confidence interval
  pval = F,              # Add p-value
  pval.coord = c(25, 0.8),
  risk.table = F,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs = substring(names(fit$strata),7), # remove first 7 character
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_classic() + theme(axis.text.x = element_text(size = 12),
                                    axis.text.y = element_text(size = 12))     # Change ggplot2 theme
)

