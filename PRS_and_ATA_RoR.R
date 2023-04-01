# Analysing PRS association with thyroid cancer subtypes. 
# Analyzing association of PRS with ATA RoR
# Nikita Pozdeyev and Martin Barrio, Nov-Dec 2022

library(MASS)
library(readxl)


ror <- read_excel('/media/nikitos/A820EB5420EB2850/Users/nikit/OneDrive - The University of Colorado Denver/Thyroid_AI/manuscript/CNN_PRS_draft/Supplementary_Tables.xlsx',
                  sheet = 4)

##########################################################
# Studying association of PRS with ATA RoR
#########################################################

# excluding benign and MTC
nobenign <- ror[!(ror$diagnosis %in% c("BENIGN", "MTC")), ]
nobenign <- nobenign[!(nobenign$ATA_RoR== "N/A"), ]

nobenign$RoR <- factor(nobenign$ATA_RoR, levels = c("Low Risk", "Intermediate Risk", "High Risk"))

table(nobenign$RoR)


m <- polr(RoR ~ PRS, data = nobenign, Hess=TRUE)
summary(m)

ctable <- coef(summary(m))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
ctable <- cbind(ctable, "p value" = p)

ci <- confint(m)


# Building a model with more covariates

nobenign$`Race/Ethnicity` <- factor(nobenign$`Race/Ethnicity`, levels = c("White or Caucasian", "Non-White", "Hispanic"))
table(nobenign$`Race/Ethnicity`)
nobenign$Race <- nobenign$`Race/Ethnicity`

m <- polr(RoR ~ PRS + Sex + Race + Age + mega1, data = nobenign, Hess=TRUE)
summary(m)

ctable <- coef(summary(m))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
ctable <- cbind(ctable, "p value" = p)
ci <- confint(m)