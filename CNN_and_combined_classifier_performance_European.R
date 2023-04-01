# Calculating confusion matrix of deep learning thyroid nodule classifier
# 95% confidence interval is computed using BCa bootstrapping method
# INDIVIDUALS OF EUROPEAN GENETIC ANCESTRY ONLY

library(caret)
library(boot)
library(readxl)
library(pROC)


# calculates threshold that produces desired sensitvity (e.g. >95% to match the sensitivity of FNA biopsy)
# calculates confusion matrix
calcConfusionMatrix <- function (rocobj, desiredSpecificity){
  
  
  #plotting threshold at >95% sensitivity
  sensitivity <- rocobj$sensitivities[sum(rocobj$sensitivities > desiredSpecificity)]
  specificity <- rocobj$specificities[sum(rocobj$sensitivities > desiredSpecificity)]
  threshold <- rocobj$thresholds[sum(rocobj$sensitivities > desiredSpecificity)]
  
  tp <- sum(rocobj$cases > threshold)
  fp <- sum(rocobj$controls > threshold)
  tn <- sum(rocobj$controls <= threshold)
  fn <- sum(rocobj$cases <= threshold)
  
  ppv <- tp/(tp + fp)
  npv <- tn/(tn+fn)
  
  cm <- data.frame (sensitivity = sensitivity,
                    specificity = specificity,
                    threshold = threshold,
                    ppv = ppv,
                    npv = npv)
  
  return(list(threshold=threshold, cm=cm))
}

# calculates 95% CI for confusion matric using BCa bootstrapping
# df should contain two columns: df$prediction and df$reference
bootstrapConfusionMatrix <- function(df) {
  
  sens <- function(df, indices){
    df <- df[indices,]
    return(confusionMatrix(data=df$prediction, reference = df$reference)$byClass["Sensitivity"])
  }
  
  spec <- function(df, indices){
    df <- df[indices,]
    return(confusionMatrix(data=df$prediction, reference = df$reference)$byClass["Specificity"])
  }
  
  npv <- function(df, indices){
    df <- df[indices,]
    return(confusionMatrix(data=df$prediction, reference = df$reference)$byClass["Neg Pred Value"])
  }
  
  ppv <- function(df, indices){
    df <- df[indices,]
    return(confusionMatrix(data=df$prediction, reference = df$reference)$byClass["Pos Pred Value"])
  }
  
  sensitivity <- boot(data = df, statistic = sens, R= 5000)
  specificty <- boot(data = df, statistic = spec, R= 5000)
  NegPredValue <- boot(data = df, statistic = npv, R= 5000)
  PosPredValue <- boot(data = df, statistic = ppv, R= 5000)
  
  return(list(sens = sens(df),
              sensci = boot.ci(sensitivity, type="bca"),
              spec = spec(df),
              specdf = boot.ci(specificty, type="bca"),
              npv = npv(df),
              npvci = boot.ci(NegPredValue, type="bca"),
              ppv = ppv(df),
              ppvci = boot.ci(PosPredValue, type="bca")))
  
}


######  confusion matrix for  crossvalidated predictions on the training/validation data  ######

# change path to point to the Supplementary_Tables.xlsx from the paper
dat <- read_excel('/media/nikitos/A820EB5420EB2850/Users/nikit/OneDrive - The University of Colorado Denver/Thyroid_AI/manuscript/CNN_PRS_draft/Supplementary_Tables.xlsx',
                  sheet = 2)

rocobj1 <- plot.roc(dat$Machine_learning_class, dat$Prob_malign_by_nodule,
                    main="",
                    percent=TRUE,
                    col="#1c61b6")
rocobj1$auc

thresholds <- calcConfusionMatrix(rocobj1, 95)

dat$biopsy <- NA
dat$biopsy[dat$Prob_malign_by_nodule >= thresholds$threshold]  <- "MALIGNANT"
dat$biopsy[dat$Prob_malign_by_nodule < thresholds$threshold]   <- "BENIGN"
bootstrapConfusionMatrix(df = data.frame(prediction = factor(dat$biopsy, levels = c("MALIGNANT", "BENIGN")),
                                         reference = factor(dat$Machine_learning_class, levels = c("MALIGNANT", "BENIGN"))))


################################################################
######  confusion matrix for out of sample test data set  ######
################################################################

# change path to point to the Supplementary_Tables.xlsx from the paper
dat <- read_excel('/media/nikitos/A820EB5420EB2850/Users/nikit/OneDrive - The University of Colorado Denver/Thyroid_AI/manuscript/CNN_PRS_draft/Supplementary_Tables.xlsx',
                  sheet = 3)

dat <- dat[dat$UMAP_RACE == "EUR", ]

# calculating AUROC
rocobj <- plot.roc(dat$class, dat$Prob_malign_by_nodule,
                   main="",
                   percent=TRUE,
                   col="#1c61b6")
rocobj$auc

dat$biopsy <- NA
dat$biopsy[dat$Prob_malign_by_nodule >= thresholds$threshold]  <- "MALIGNANT"
dat$biopsy[dat$Prob_malign_by_nodule < thresholds$threshold]   <- "BENIGN"
bootstrapConfusionMatrix(df = data.frame(prediction = factor(dat$biopsy, levels = c("MALIGNANT", "BENIGN")),
                                         reference = factor(dat$ML_class, levels = c("MALIGNANT", "BENIGN"))))





###############################
######  CNN + crude PRS  ######
###############################

set.seed(123)
nfolds = 5

ctrl <- trainControl(method = "cv", number = nfolds, classProbs = TRUE, summaryFunction = twoClassSummary)
formula <- "ML_class ~ Prob_malign_by_nodule + PRSscaled + mega1"

logit.CV <- train(as.formula(formula), data = dat, method = "glm", trControl = ctrl, family = "binomial", metric = "ROC")
# logit.CV <- train(as.formula(formula), data = dat, method = "cforest", trControl = ctrl, metric = "ROC")
# logit.CV <- train(as.formula(formula), data = dat, method = "AdaBag", trControl = ctrl, metric = "ROC")
# logit.CV <- train(as.formula(formula), data = dat, method = "adaboost", trControl = ctrl, metric = "ROC")

dat$pred_class <- predict(logit.CV, newdata = dat, "prob")$MALIGNANT


# distribution of probabilities of the DL classifier and DL+PRS classifier are different 
# to be able to compare confusion matrix metrics, for each model selecting threshold that achieved the sensitivity of ~ 0.95
# matching that of the FNA. AUROC calculation are not affected by the threshold choice

rocobj <- plot.roc(dat$ML_class, dat$pred_class,
                   main="",
                   percent=TRUE,
                   col="#1c61b6")

rocobj$auc

thresholds <- calcConfusionMatrix(rocobj, 95)

dat$biopsy <- NA
dat$biopsy[dat$pred_class >= thresholds$threshold]  <- "MALIGNANT"
dat$biopsy[dat$pred_class < thresholds$threshold]   <- "BENIGN"
bootstrapConfusionMatrix(df = data.frame(prediction = factor(dat$biopsy, levels = c("MALIGNANT", "BENIGN")),
                                         reference = factor(dat$ML_class, levels = c("MALIGNANT", "BENIGN"))))




# CNN alone or in combination with PRS plus PCs
pdf("/media/nikitos/A820EB5420EB2850/Users/nikit/OneDrive - The University of Colorado Denver/Thyroid_AI/manuscript/CNN_PRS_draft/figures/CNN_plus_PRS_EUR.pdf",
    width = 5, height = 5)
rocobj1 <- plot.roc(dat$class, dat$Prob_malign_by_nodule,
                    main="",
                    percent=TRUE,
                    col="#1c61b6")



rocobj2 <- lines.roc(dat$class, dat$pred_class, 
                     percent=TRUE, 
                     col="#FF0000")

cm1 <- calcConfusionMatrix(rocobj1, 95)
cm2 <- calcConfusionMatrix(rocobj2, 95)

segments(100, cm2$cm$sensitivity, cm2$cm$specificity, cm2$cm$sensitivity, col = "green")
segments(cm2$cm$specificity, cm2$cm$sensitivity, cm2$cm$specificity, 0, col = "green")


testobj <- roc.test(rocobj1, rocobj2, method = "delong", paired = T, alternative = "less")
text(80, 40, labels=paste("DeLong p-value =", format.pval(testobj$p.value, digits = 3)), adj=c(0, .5))
text(80, 50, labels=paste("AUROC CNN =", format(rocobj1$auc/100, digits = 3)), adj=c(0, .5))
text(80, 45, labels=paste("AUROC CNN + PRS =", format(rocobj2$auc/100, digits = 3)), adj=c(0, .5))
legend("bottomright", legend=c("CNN", "CNN + PRS "), col=c("#1c61b6", "#FF0000"), lwd=2)
dev.off()




###############################
######  CNN + PRS + PC   ######
###############################

set.seed(123)
nfolds = 5

ctrl <- trainControl(method = "cv", number = nfolds, classProbs = TRUE, summaryFunction = twoClassSummary)
formula <- "ML_class ~ Prob_malign_by_nodule + PRSscaled + mega1 + PC1 + PC2 + PC3 + PC4 + PC5"

logit.CV <- train(as.formula(formula), data = dat, method = "glm", trControl = ctrl, family = "binomial", metric = "ROC")
# logit.CV <- train(as.formula(formula), data = dat, method = "cforest", trControl = ctrl, metric = "ROC")
# logit.CV <- train(as.formula(formula), data = dat, method = "AdaBag", trControl = ctrl, metric = "ROC")
# logit.CV <- train(as.formula(formula), data = dat, method = "adaboost", trControl = ctrl, metric = "ROC")

dat$pred_class <- predict(logit.CV, newdata = dat, "prob")$MALIGNANT


# distribution of probabilities of the DL classifier and DL+PRS classifier are different 
# to be able to compare confusion matrix metrics, for each model selecting threshold that achieved the sensitivity of ~ 0.95
# matching that of the FNA. AUROC calculation are not affected by the threshold choice

rocobj <- plot.roc(dat$ML_class, dat$pred_class,
                   main="",
                   percent=TRUE,
                   col="#1c61b6")

rocobj$auc

thresholds <- calcConfusionMatrix(rocobj, 95)

dat$biopsy <- NA
dat$biopsy[dat$pred_class >= thresholds$threshold]  <- "MALIGNANT"
dat$biopsy[dat$pred_class < thresholds$threshold]   <- "BENIGN"
bootstrapConfusionMatrix(df = data.frame(prediction = factor(dat$biopsy, levels = c("MALIGNANT", "BENIGN")),
                                         reference = factor(dat$ML_class, levels = c("MALIGNANT", "BENIGN"))))




# CNN alone or in combination with PRS plus PCs
pdf("/media/nikitos/A820EB5420EB2850/Users/nikit/OneDrive - The University of Colorado Denver/Thyroid_AI/manuscript/CNN_PRS_draft/figures/CNN_plus_PRS_PCs_EUR.pdf",
    width = 5, height = 5)
rocobj1 <- plot.roc(dat$class, dat$Prob_malign_by_nodule,
                    main="",
                    percent=TRUE,
                    col="#1c61b6")



rocobj2 <- lines.roc(dat$class, dat$pred_class, 
                     percent=TRUE, 
                     col="#FF0000")

cm1 <- calcConfusionMatrix(rocobj1, 95)
cm2 <- calcConfusionMatrix(rocobj2, 95)

segments(100, cm2$cm$sensitivity, cm2$cm$specificity, cm2$cm$sensitivity, col = "green")
segments(cm2$cm$specificity, cm2$cm$sensitivity, cm2$cm$specificity, 0, col = "green")


testobj <- roc.test(rocobj1, rocobj2, method = "delong", paired = T, alternative = "less")
text(80, 40, labels=paste("DeLong p-value =", format.pval(testobj$p.value, digits = 1)), adj=c(0, .5))
text(80, 50, labels=paste("AUROC CNN =", format(rocobj1$auc/100, digits = 3)), adj=c(0, .5))
text(80, 45, labels=paste("AUROC CNN + PRS + PCs =", format(rocobj2$auc/100, digits = 3)), adj=c(0, .5))
legend("bottomright", legend=c("CNN", "CNN + PRS + PCs"), col=c("#1c61b6", "#FF0000"), lwd=2)
dev.off()




#######################################
######  CNN + PRS + COVARIATES   ######
#######################################

set.seed(123)
nfolds = 5

ctrl <- trainControl(method = "cv", number = nfolds, classProbs = TRUE, summaryFunction = twoClassSummary)
formula <- "ML_class ~ Prob_malign_by_nodule + PRSscaled + H_AP + W_Transverse + L_Longitudinal + mega1 + Sex + age + PC1 + PC2 + PC3 + PC4 + PC5"

logit.CV <- train(as.formula(formula), data = dat, method = "glm", trControl = ctrl, family = "binomial", metric = "ROC")
# logit.CV <- train(as.formula(formula), data = dat, method = "cforest", trControl = ctrl, metric = "ROC")
# logit.CV <- train(as.formula(formula), data = dat, method = "AdaBag", trControl = ctrl, metric = "ROC")
# logit.CV <- train(as.formula(formula), data = dat, method = "adaboost", trControl = ctrl, metric = "ROC")

dat$pred_class <- predict(logit.CV, newdata = dat, "prob")$MALIGNANT


# distribution of probabilities of the DL classifier and DL+PRS classifier are different 
# to be able to compare confusion matrix metrics, for each model selecting threshold that achieved the sensitivity of ~ 0.95
# matching that of the FNA. AUROC calculation are not affected by the threshold choice

rocobj <- plot.roc(dat$ML_class, dat$pred_class,
                   main="",
                   percent=TRUE,
                   col="#1c61b6")

rocobj$auc

thresholds <- calcConfusionMatrix(rocobj, 95)

dat$biopsy <- NA
dat$biopsy[dat$pred_class >= thresholds$threshold]  <- "MALIGNANT"
dat$biopsy[dat$pred_class < thresholds$threshold]   <- "BENIGN"
bootstrapConfusionMatrix(df = data.frame(prediction = factor(dat$biopsy, levels = c("MALIGNANT", "BENIGN")),
                                         reference = factor(dat$ML_class, levels = c("MALIGNANT", "BENIGN"))))




# CNN alone or in combination with PRS plus PCs
pdf("/media/nikitos/A820EB5420EB2850/Users/nikit/OneDrive - The University of Colorado Denver/Thyroid_AI/manuscript/CNN_PRS_draft/figures/CNN_plus_PRS_Covariates_EUR.pdf",
    width = 5, height = 5)
rocobj1 <- plot.roc(dat$class, dat$Prob_malign_by_nodule,
                    main="",
                    percent=TRUE,
                    col="#1c61b6")



rocobj2 <- lines.roc(dat$class, dat$pred_class, 
                     percent=TRUE, 
                     col="#FF0000")

cm1 <- calcConfusionMatrix(rocobj1, 95)
cm2 <- calcConfusionMatrix(rocobj2, 95)

segments(100, cm2$cm$sensitivity, cm2$cm$specificity, cm2$cm$sensitivity, col = "green")
segments(cm2$cm$specificity, cm2$cm$sensitivity, cm2$cm$specificity, 0, col = "green")


testobj <- roc.test(rocobj1, rocobj2, method = "delong", paired = T, alternative = "less")
text(80, 40, labels=paste("DeLong p-value =", format.pval(testobj$p.value, digits = 1)), adj=c(0, .5))
text(80, 50, labels=paste("AUROC CNN =", format(rocobj1$auc/100, digits = 3)), adj=c(0, .5))
text(80, 45, labels=paste("AUROC CNN + PRS + cov =", format(rocobj2$auc/100, digits = 3)), adj=c(0, .5))
legend("bottomright", legend=c("CNN", "CNN + PRS + cov"), col=c("#1c61b6", "#FF0000"), lwd=2)
dev.off()
















