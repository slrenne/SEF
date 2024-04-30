# libraries
library(survminer)
library(survival)

# descriptive statistics
table(db$DIAGNOSIS)
round(proportions(table(db$DIAGNOSIS)),2)
sum(!is.na(db$DIAGNOSIS))
quantile(db$AGE, c(0.25,0.5,0.75))
sum(!is.na(db$AGE))
table(db$SEX)
sum(!is.na(db$SEX))
round(proportions(table(db$SEX)),2)
db$`SIZE_(cm)`[18] <- 9.8 
db$`SIZE_(cm)`[96] <- 15.1 
db$`SIZE_(cm)`[143] <- 5.1 
db$`SIZE_(cm)`[144] <- 5.1 
quantile(as.numeric(db$`SIZE_(cm)`), c(0.25,0.5,0.75), na.rm = TRUE)
sum(!is.na(db$`SIZE_(cm)`))
table(db$SITE)
#SITE needs to be semplified
# round(proportions(table(db$SITE)),2)
db$`MITOSIS/10HPF`[db$`MITOSIS/10HPF`%in% c("Absent","N","<1")] <- 0
db$`MITOSIS/10HPF`[146] <- 5 
db$`MITOSIS/10HPF`[149] <- 1 
db$`MITOSIS/10HPF`[128] <- 10 
quantile(as.numeric(db$`MITOSIS/10HPF`), c(0.25,0.5,0.75), na.rm = TRUE)
sum(!is.na(db$`MITOSIS/10HPF`))
table(db$NECROSIS_ptg)
round(proportions(table(db$NECROSIS_ptg)),2)
sum(!is.na(db$NECROSIS_ptg))
table(db$MUC4)
round(proportions(table(db$MUC4)),2)
sum(!is.na(db$MUC4))

for (i in 13:28) {
  print(paste(colnames(db)[i],"tot =", sum(!is.na(db[,i]))));
  print(table(db[,i]));
  print(round(proportions(table(db[,i])),2));
}

table(db$`Fusion_1canonical_2YAP1-KMT2A`)
round(proportions(table(db$`Fusion_1canonical_2YAP1-KMT2A`)),2)
sum(!is.na(db$`Fusion_1canonical_2YAP1-KMT2A`))
   
event <- as.integer(db$LFUP=="DOD")
time <- as.numeric(db$`LFUP_(Months)`)


fit <- survfit(Surv(
  time = time,
  event = event) ~ 1, 
  data = db)
ggsurvplot(fit ,  
           title = "Overall Survival", 
           xlab = "months", 
           legend = "none",
           conf.int = FALSE,
           data = db)

fit <- survfit(Surv(
  time = time,
  event = event) ~ db$`Fusion_1canonical_2YAP1-KMT2A`, 
  data = db)

ggsurvplot(fit ,  
           title = "Overall Survival", 
           xlab = "months", data = db, 
           fun = "pct",
           risk.table = TRUE,
           conf.int = TRUE,
           legend.labs = c('Canonical','YAP1-KAMT2A'),
           legend.title = "Fusion",
           pval = TRUE)
