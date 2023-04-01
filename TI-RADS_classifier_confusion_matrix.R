# this script is calculating confusion matrix with confidence intervals for TI-RADS classifier
# source data is Table 1 of the Supplementary data


# library(data.table)
library(caret)
library(boot)
library(readxl)

# change path to point to the Supplementary_Tables.xlsx from the paper
dat <- read_excel('/media/nikitos/A820EB5420EB2850/Users/nikit/OneDrive - The University of Colorado Denver/Thyroid_AI/manuscript/CNN_PRS_draft/Supplementary_Tables.xlsx', 
                   sheet = 2)

# subsetting to samples with sufficient data for TI-RADS classification in the ultrasound report
dat <- dat[!is.na(dat$Ti_RADS_score), ]

# Estimating maximal size of the nodule
tmp <- dat[, c("Size_AP", "Size_longitudinal", "Size_transverse")]
dat$maxsize <- apply(tmp, MARGIN = 1, max)

# # calculating recommendation to proceed with biopsy based on the TI-RADS score and nodule size
# dat$biopsy <- NA
# dat$biopsy[dat$Ti_RADS_score == 5 & dat$maxsize >= 1]  <- "Yes"
# dat$biopsy[dat$Ti_RADS_score == 4 & dat$maxsize >= 1.5]  <- "Yes"
# dat$biopsy[dat$Ti_RADS_score == 3 & dat$maxsize >= 2.5]  <- "Yes"
# dat$biopsy[is.na(dat$biopsy)] <- "No"
# 
# dat1 <- fread("Supplementary_Table 1.csv")
# dat2 <- dat
# dat2 <- dat2[, c("Nodule_ID", "biopsy")]
# dat1 <- merge(dat1, dat2, by = "Nodule_ID", all.x = T, sort = F)
# write.csv(dat1, "Supplementary_Table_1_with_biopsy_rec.csv", row.names = F)
# 

# calculating recommendation to proceed with biopsy based on the TI-RADS score and nodule size
dat$biopsy <- NA
dat$biopsy[dat$Ti_RADS_score == 5 & dat$maxsize >= 1]  <- "MALIGNANT"
dat$biopsy[dat$Ti_RADS_score == 4 & dat$maxsize >= 1.5]  <- "MALIGNANT"
dat$biopsy[dat$Ti_RADS_score == 3 & dat$maxsize >= 2.5]  <- "MALIGNANT"
dat$biopsy[is.na(dat$biopsy)] <- "BENIGN"

table(dat$Machine_learning_class)

dat$biopsy <- factor(dat$biopsy, levels=c("MALIGNANT", "BENIGN"))
dat$Machine_learning_class <- factor(dat$Machine_learning_class, levels=c("MALIGNANT", "BENIGN"))

# calculating confusion matrix
cm <- confusionMatrix(data=dat$biopsy, reference = dat$Machine_learning_class)

cm$byClass["Sensitivity"]

sens <- function(df, indices){
  df <- df[indices,]
  return(confusionMatrix(data=df$biopsy, reference = df$Machine_learning_class)$byClass["Sensitivity"])
}

sens(dat)

sensitivity <- boot(data = dat, statistic = sens, R= 10000)
boot.ci(sensitivity, type="bca")


spec <- function(df, indices){
  df <- df[indices,]
  return(confusionMatrix(data=df$biopsy, reference = df$Machine_learning_class)$byClass["Specificity"])
}

spec(dat)

specificity <- boot(data = dat, statistic = spec, R= 10000)
boot.ci(specificity, type="bca")


npv <- function(df, indices){
  df <- df[indices,]
  return(confusionMatrix(data=df$biopsy, reference = df$Machine_learning_class)$byClass["Neg Pred Value"])
}

npv(dat)

negpv <- boot(data = dat, statistic = npv, R= 10000)
boot.ci(negpv, type="bca")




ppv <- function(df, indices){
  df <- df[indices,]
  return(confusionMatrix(data=df$biopsy, reference = df$Machine_learning_class)$byClass["Pos Pred Value"])
}

ppv(dat)

pospv <- boot(data = dat, statistic = ppv, R= 10000)
boot.ci(pospv, type="bca")

