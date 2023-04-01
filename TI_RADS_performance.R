# This script evaluate performance of TI-RADS classifier as reported by radiologists 
# for thyroid nodules in the UW trainin/validation set


library(ggplot2)
library(ROCR)
# library(pROC)
library(readxl)

rm(list=ls())


# change path to point to the Supplementary_Tables.xlsx from the paper
meta <- read_excel('/media/nikitos/A820EB5420EB2850/Users/nikit/OneDrive - The University of Colorado Denver/Thyroid_AI/manuscript/CNN_PRS_draft/Supplementary_Tables.xlsx', 
                   sheet = 2)

meta_tr <- meta[!is.na(meta$Ti_RADS_points) & !is.na(meta$Ti_RADS_score), ]


# Kernel Density Plot for Ti-RADS points

colors <- c("BENIGN" = "#FFFFFF", "MALIGNANT" = "#333333")
ggplot(data = meta_tr, aes(x = Ti_RADS_points, fill = Machine_learning_class)) + geom_density(size = 1, alpha = 0.3) + 
  scale_fill_manual(values = colors, name = 'Class') + 
  theme(panel.background = element_rect(fill = 'white', colour = 'white', size = 0)) + 
#  labs(x = 'Ti-RADS Points', y = 'Density', size = 16) +
  scale_x_continuous(breaks = seq(0,13,2), name = "TI-RADS points") +
  scale_y_continuous(breaks = seq(0,0.25,0.05), name = "Density") +
  theme(axis.text.x  = element_text(angle=0, vjust=8, hjust = 0.5, size=16), axis.ticks.x = element_blank()) +
  theme(axis.title.x = element_text(vjust = 6, size = 16)) +
  theme(axis.text.y  = element_text(angle=0, size=16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16))

  
# Kernel Density Plot for Ti-RADS score

colors <- c("BENIGN" = "#FFFFFF", "MALIGNANT" = "#333333")
ggplot(data = meta_tr, aes(x = Ti_RADS_score, fill = Machine_learning_class)) + 
  geom_density(size = 1, alpha = 0.3, adjust = 2) + 
  scale_fill_manual(values = colors, name = 'Class') + 
  theme(panel.background = element_rect(fill = 'white', colour = 'white', size = 0)) + 
  #  labs(x = 'Ti-RADS Points', y = 'Density', size = 16) +
  scale_x_continuous(breaks = seq(0,5,1), name = "Ti-RADS category") +
  scale_y_continuous(breaks = seq(0,0.4,0.1), limits = c(0,0.5), name = "Density") +
  theme(axis.text.x  = element_text(angle=0, vjust=8, hjust = 0.5, size=16), axis.ticks.x = element_blank()) +
  theme(axis.title.x = element_text(vjust = 6, size = 16)) +
  theme(axis.text.y  = element_text(angle=0, size=16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16))


# constructing ROC for TI-RADS points

pred <-  prediction(meta_tr$Ti_RADS_points, meta_tr$Machine_learning_class)
perf <- performance(pred, "tpr", "fpr")
plot(perf)

auc <- performance(pred, "auc")
text(x = 0.6, y = 0.2, labels = paste0("ROC AUC = ", round(auc@y.values[[1]], digits = 3)))

# ROC for TI_RADS category


pred <-  prediction(meta_tr$Ti_RADS_score, meta_tr$Machine_learning_class)
perf <- performance(pred, "tpr", "fpr")
plot(perf)

auc <- performance(pred, "auc")
text(x = 0.6, y = 0.2, labels = paste0("ROC AUC = ", round(auc@y.values[[1]], digits = 3)))





