
library(data.table)
newIdhLabelPwd = './mgbmSgbmLabelMerge230907.csv'
newIdhLabelDf = data.frame(fread(newIdhLabelPwd))

newIdhLabel4sgbmList = newIdhLabelDf$ID[newIdhLabelDf$idhMannual2 == 'WT' & newIdhLabelDf$bigSmall == 'SGBM']
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
library(data.table)
mgbmLabelPwd = './mgbmSgbmLabelMerge230907.csv'
mgbmLabelDf = data.frame(fread(mgbmLabelPwd))

table(mgbmLabelDf$bigSmall)

group1Patient = mgbmLabelDf$ID[mgbmLabelDf$bigSmall == 'Big']
group1Label = 'Big'
group2Patient = mgbmLabelDf$ID[mgbmLabelDf$bigSmall == 'Small']
group2Label = 'Small'
group3Patient = newIdhLabel4sgbmList
group3Label = 'SGBM'
cat(group1Patient, sep = ',')
cat(group2Patient, sep = ',')
cat(group3Patient, sep = ',')

length(unique(c(group1Patient, group2Patient)))
length(unique(c(group1Patient, group2Patient, group3Patient)))
newClinicalGbmAll = data.frame(fread('./data_clinical_patient.txt', skip = 4))

unique(newClinicalGbmAll$DFS_STATUS)
newClinicalGbmAll$DFS_STATUS[newClinicalGbmAll$DFS_STATUS == "1:Recurred/Progressed"] = 1
newClinicalGbmAll$DFS_STATUS[newClinicalGbmAll$DFS_STATUS == "0:DiseaseFree"] = 0
newClinicalGbmAll$DFS_STATUS[newClinicalGbmAll$DFS_STATUS == "[Not Available]"] = NA
newClinicalTarget = newClinicalGbmAll[newClinicalGbmAll$PATIENT_ID %in% c(group1Patient, group2Patient, group3Patient), ]
groupGDC_newClinicalTarget = newClinicalTarget[,c('PATIENT_ID', 'DFS_MONTHS', 'DFS_STATUS')]

dim(unique(groupGDC_newClinicalTarget))

survivalTimeList = groupGDC_newClinicalTarget$DFS_MONTHS
survivalTimeList[survivalTimeList == "[Not Available]"] = NA
survivalTimeList = as.numeric(survivalTimeList)
survivalTimeList

groupGDC_newClinicalTarget$time = survivalTimeList
groupGDC_newClinicalTarget = unique(groupGDC_newClinicalTarget)
groupGDC_newClinicalTarget = groupGDC_newClinicalTarget[!is.na(groupGDC_newClinicalTarget$time),]

colnames(groupGDC_newClinicalTarget)

dim(groupGDC_newClinicalTarget)

setdiff(c(group1Patient, group2Patient), groupGDC_newClinicalTarget$PATIENT_ID)
setdiff(c(group1Patient, group2Patient, group3Patient), groupGDC_newClinicalTarget$PATIENT_ID)

groupGDC_newClinicalTarget$group = 'unknown'
groupGDC_newClinicalTarget$group[groupGDC_newClinicalTarget$PATIENT_ID %in% group1Patient] = group1Label
groupGDC_newClinicalTarget$group[groupGDC_newClinicalTarget$PATIENT_ID %in% group2Patient] = group2Label
groupGDC_newClinicalTarget$group[groupGDC_newClinicalTarget$PATIENT_ID %in% group3Patient] = group3Label
twoGroupDf = groupGDC_newClinicalTarget
table(twoGroupDf$DFS_STATUS)
twoGroupDf$state = 0
twoGroupDf$state[twoGroupDf$DFS_STATUS == 0] = 1 # live 
twoGroupDf$state[twoGroupDf$DFS_STATUS ==1] = 2 # dead
twoGroupDf$group = factor(twoGroupDf$group,
                          levels = c(group1Label, group2Label, group3Label))

fit = survfit(Surv(time, state) ~ group, data = twoGroupDf)

diff = survdiff(Surv(time, state) ~ group, data = twoGroupDf)
pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE) # get the t

unique(twoGroupDf$group)

diff = survdiff(Surv(time, state) ~ group, data = twoGroupDf[twoGroupDf$group %in% c('Big', 'SGBM'), ])
pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE) 

diff = survdiff(Surv(time, state) ~ group, data = twoGroupDf[twoGroupDf$group %in% c('Small', 'SGBM'), ])
pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE) 

diff = survdiff(Surv(time, state) ~ group, data = twoGroupDf[twoGroupDf$group %in% c('Big', 'Small'), ])
pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE) 
ggsurv = ggsurvplot(
  fit, 
  data = twoGroupDf, 
  size = 1,                 # change line size
  ylab = 'Probability of Disease Free Survival',
  xlab = ' Disease Free Survival (Months)',
  palette = c('#D63821', '#224B89', '#2ca25f'),
  conf.int = F,          # Add confidence interval
  pval.coord = c(50, 0.1),
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
  ylab = 'Probability of Disease Free Survival',
  xlab = ' Disease Free Survival (Months)',
  palette = c('#D63821', '#224B89', '#2ca25f'),
  conf.int = F,          # Add confidence interval
  pval.coord = c(50, 0.1),
  risk.table = F,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs = substring(names(fit$strata),7), # remove first 7 character
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_classic() + theme(axis.text.x = element_text(size = 12),
                                    axis.text.y = element_text(size = 12))     # Change ggplot2 theme
)


