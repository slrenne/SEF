# libraries
library(survminer)
library(survival)
library(rethinking)

# descriptive statistics
table(db$DIAGNOSIS)
round(proportions(table(db$DIAGNOSIS)),2)
sum(!is.na(db$DIAGNOSIS))
quantile(db$AGE, c(0.25,0.5,0.75))
sum(!is.na(db$AGE))
table(db$SEX)
sum(!is.na(db$SEX))
round(proportions(table(db$SEX)),2)
db[18,8] <- 9.8 
db[96,8] <- 15.1 
db[143,8] <- 5.1 
db[144,8] <- 5.1 
quantile(as.numeric(db[,8]), c(0.25,0.5,0.75), na.rm = TRUE) # size
sum(!is.na(db[,8]))
round(proportions(table(db$SEF_TYPE)),2)
#SITE needs to be semplified
# round(proportions(table(db$SITE)),2)
db$`MITOSIS/10HPF`[db$`MITOSIS/10HPF`%in% c("Absent","N","<1")] <- 0
db[,11][db[,11] %in% c("Absent","N","<1")] <- 0
db[146,11] <- 5 
db[149,11] <- 1 
db[128,11] <- 10 
db[58,50] <- 12
quantile(as.numeric(db[,11]), c(0.25,0.5,0.75), na.rm = TRUE)
sum(!is.na(db[,11]))
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

table(db[,29])
round(proportions(table(db[,29])),2)
sum(!is.na(db[,29]))
   
event <- as.integer(db$LFUP=="DOD")
time <- as.numeric(db[,50])


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
  event = event) ~ db[,29], 
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

# local recurrence survival analysis

event <- as.integer(db$LOCAL_RECURRENCE=="Y")
time <- ifelse(event==1, 
               as.numeric(db[,42]),
               as.numeric(db[,50]))


fit <- survfit(Surv(
  time = time,
  event = event) ~ 1, 
  data = db)
ggsurvplot(fit ,  
           title = "Local Recurrence", 
           xlab = "months", 
           legend = "none",
           conf.int = FALSE,
           data = db)

fit <- survfit(Surv(
  time = time,
  event = event) ~ db[,29], 
  data = db)

ggsurvplot(fit ,  
           title = "Local Recurrence", 
           xlab = "months", data = db, 
           fun = "pct",
           risk.table = TRUE,
           conf.int = TRUE,
           legend.labs = c('Canonical','YAP1-KAMT2A'),
           legend.title = "Fusion",
           pval = TRUE)

# distant metastasis survival analysis

event <- as.integer(db$DISTANT__METASTASIS_!="N")
time <- ifelse(event==1, 
               as.numeric(db[,44]),
               as.numeric(db[,50]))


fit <- survfit(Surv(
  time = time,
  event = event) ~ 1, 
  data = db)
ggsurvplot(fit ,  
           title = "Distant Metastasis", 
           xlab = "months", 
           legend = "none",
           conf.int = FALSE,
           data = db)

fit <- survfit(Surv(
  time = time,
  event = event) ~ db[,29], 
  data = db)

ggsurvplot(fit ,  
           title = "Distant Metastasis", 
           xlab = "months", data = db, 
           fun = "pct",
           risk.table = TRUE,
           conf.int = TRUE,
           legend.labs = c('Canonical','YAP1-KAMT2A'),
           legend.title = "Fusion",
           pval = TRUE)

db$PMID[167:174] <- "institutional"

# multilevel modeling
# age
!is.na(db$AGE)
dat <- list(
  A = standardize(db$AGE),
  S = as.integer(as.factor(db$PMID)), 
  F = as.integer(db$Fusion)
)

m.age <- ulam(
  alist(
    A ~ dnorm( mu, sigma_a),
    mu <- a[S] + b[F],
    # define effects using other parameters
    save> vector[38]:a <<- abar + za*sigma,
    save> vector[2]:b <<- bbar + zb*tau,
    # z-scored effects
    vector[38]:za ~ normal(0,1),
    vector[2]:zb ~ normal(0,1),
    # hyper-priors
    c(abar,bbar) ~ normal(0,1),
    c(sigma,tau,sigma_a) ~ exponential(1)
  ), data = dat, chains = 4, cores = 4, iter = 2000)

dashboard(m.age)
plot(precis(m.age, 2, pars = c("a","b", "abar")))
post <- extract.samples(m.age)
s.center <- attributes(dat$A)[[1]]
s.scale <-  attributes(dat$A)[[2]]
de_stand <- function(x){  s.center + x * s.scale }
par(mar= c(5, 7, 4, 2) + 0.1)

plot(NULL, xlim = de_stand( c(-1.7,1.5) ), ylim = c(0.0,38.5),
     xlab = 'Age', ylab = '', yaxt = 'null', 
     main = 'Age in SEF')
abline(v=de_stand(0),lty=2)
abline(h= 0.5)
idx <- order(colMeans(post$a), decreasing = FALSE)
for ( i in 1:38 ) {
  x0g <- de_stand(quantile(post$a[,idx[i]], 0.05))
  x1g <- de_stand(quantile(post$a[,idx[i]], 0.95))
  segments(x0=x0g, y0=i, x1=x1g)
  x0g <- de_stand(quantile(post$a[,idx[i]], 0.25))
  x1g <- de_stand(quantile(post$a[,idx[i]], 0.75))
  segments(x0=x0g, y0=i, x1=x1g, lwd = 5, col = 2)
}
label <- levels(as.factor(db$PMID))
text(y = 38:1,  de_stand(-2),
     labels = label[idx],  
     srt=0,  adj=1,    xpd=TRUE)
legend('topleft', lwd = c(1,5), legend = c('95CI','50CI'), col = 1:2)
x0g <- de_stand(quantile(post$abar+post$bbar, 0.05))
x1g <- de_stand(quantile(post$abar+post$bbar, 0.95))
segments(x0=x0g, y0=-0.5, x1=x1g, lwd= 3)
x0g <- de_stand(quantile(post$abar+post$bbar, 0.25))
x1g <- de_stand(quantile(post$abar+post$bbar, 0.75))
segments(x0=x0g, y0=-0.5, x1=x1g, lwd = 6, col = 2)
text(y = -0.5,  de_stand(-2),
     labels = "Population age",  
     srt=0,  adj=1,    xpd=TRUE)
par(mfrow = c(2,1))
plot(NULL, xlim = de_stand( c(-1.7,1.5) ), ylim = c(0 , 3),
     xlab = 'Age', ylab = '', yaxt = 'null', 
     main = 'Age in SEF')
abline(v=de_stand(0),lty=2)
for ( i in 1:2 ) {
  x0g <- de_stand(quantile(post$b[,i], 0.05))
  x1g <- de_stand(quantile(post$b[,i], 0.95))
  segments(x0=x0g, y0=i, x1=x1g)
  x0g <- de_stand(quantile(post$b[,i], 0.25))
  x1g <- de_stand(quantile(post$b[,i], 0.75))
  segments(x0=x0g, y0=i, x1=x1g, lwd = 5, col = 4)
}
text(y = 2:1,  de_stand(-2),
     labels = c("Canonical"," YAP1-KAMT2A"),  
     srt=0,  adj=1,    xpd=TRUE)
F1 <- rnorm(4000, post$abar + post$b[,1], post$sigma_a)
F2 <- rnorm(4000, post$abar + post$b[,2], post$sigma_a)
F_contrast <- F2 - F1
dens(F_contrast, adj = 1 , lwd = 3, col = 4, show.zero = TRUE, 
     main = "Contrast by fusion", xlab = "posterior age contrast")


# multilevel modeling
# size
db$SIZE_.cm. <- as.numeric(db$SIZE_.cm.)

i<- !is.na(db$SIZE_.cm.)

dat <- list(
  D = standardize(db$SIZE_.cm.[i]),
  S = as.integer(as.factor(db$PMID[i])), 
  F = as.integer(db$Fusion[i])
)

m.size <- ulam(
  alist(
    D ~ dnorm( mu, sigma_a),
    mu <- a[S] + b[F],
    # define effects using other parameters
    save> vector[38]:a <<- abar + za*sigma,
    save> vector[2]:b <<- bbar + zb*tau,
    # z-scored effects
    vector[38]:za ~ normal(0,1),
    vector[2]:zb ~ normal(0,1),
    # hyper-priors
    c(abar,bbar) ~ normal(0,1),
    c(sigma,tau,sigma_a) ~ exponential(1)
  ), data = dat, chains = 4, cores = 4, iter = 2000)

dashboard(m.size)
plot(precis(m.size, 2, pars = c("a","b", "abar")))
post <- extract.samples(m.size)
s.center <- attributes(dat$D)[[1]]
s.scale <-  attributes(dat$D)[[2]]
de_stand <- function(x){  s.center + x * s.scale }
par(mar= c(5, 7, 4, 2) + 0.1)
par(mfrow=c(1,1))
plot(NULL, xlim = de_stand( c(-1.7,1.5) ), ylim = c(0.0,38.5),
     xlab = 'Size (cm)', ylab = '', yaxt = 'null', 
     main = 'Size in SEF')
abline(v=de_stand(0),lty=2)
abline(h= 0.5)
idx <- order(colMeans(post$a), decreasing = FALSE)
for ( i in 1:38 ) {
  x0g <- de_stand(quantile(post$a[,idx[i]], 0.05))
  x1g <- de_stand(quantile(post$a[,idx[i]], 0.95))
  segments(x0=x0g, y0=i, x1=x1g)
  x0g <- de_stand(quantile(post$a[,idx[i]], 0.25))
  x1g <- de_stand(quantile(post$a[,idx[i]], 0.75))
  segments(x0=x0g, y0=i, x1=x1g, lwd = 5, col = 2)
}
label <- levels(as.factor(db$PMID))
text(y = 38:1,  de_stand(-2),
     labels = label[idx],  
     srt=0,  adj=1,    xpd=TRUE)
legend('topleft', lwd = c(1,5), legend = c('95CI','50CI'), col = 1:2)
x0g <- de_stand(quantile(post$abar+post$bbar, 0.05))
x1g <- de_stand(quantile(post$abar+post$bbar, 0.95))
segments(x0=x0g, y0=-0.5, x1=x1g, lwd= 3)
x0g <- de_stand(quantile(post$abar+post$bbar, 0.25))
x1g <- de_stand(quantile(post$abar+post$bbar, 0.75))
segments(x0=x0g, y0=-0.5, x1=x1g, lwd = 6, col = 2)
text(y = -0.5,  de_stand(-2),
     labels = "Population size",  
     srt=0,  adj=1,    xpd=TRUE)

par(mfrow = c(2,1))
plot(NULL, xlim = de_stand( c(-1.7,1.5) ), ylim = c(0 , 3),
     xlab = 'Size in cm', ylab = '', yaxt = 'null', 
     main = 'Size in SEF')
abline(v=de_stand(0),lty=2)
for ( i in 1:2 ) {
  x0g <- de_stand(quantile(post$b[,i], 0.05))
  x1g <- de_stand(quantile(post$b[,i], 0.95))
  segments(x0=x0g, y0=i, x1=x1g)
  x0g <- de_stand(quantile(post$b[,i], 0.25))
  x1g <- de_stand(quantile(post$b[,i], 0.75))
  segments(x0=x0g, y0=i, x1=x1g, lwd = 5, col = 4)
}
text(y = 2:1,  de_stand(-2),
     labels = c("Canonical"," YAP1-KAMT2A"),  
     srt=0,  adj=1,    xpd=TRUE)
F1 <- rnorm(4000, post$abar + post$b[,1], post$sigma_a)
F2 <- rnorm(4000, post$abar + post$b[,2], post$sigma_a)
F_contrast <- F2 - F1
dens(F_contrast, adj = 1 , lwd = 3, col = 4, show.zero = TRUE, 
     main = "Contrast by fusion", xlab = "posterior size contrast")


event <- as.integer(db$LFUP=="DOD")
time <- as.numeric(db[,50])


#survival analysis
#Metastasis free survival
i<- !is.na(db[,50]) & db[,50]!=0
dat <- list(
  E = as.integer(db$LFUP[i]=="DOD"),
  F = as.integer(db$Fusion[i]),
  A = standardize(db$AGE[i]),
  d = as.integer(db$DIAGNOSIS[i] != "SEF") + 1L, # 1 SEF, 2 hybrid
  months_to_event = as.numeric(db[i,50])
)
dat$L <- length(dat[[1]])


m.os <-rstan::stan(
  file = "os.stan",  # Stan program
  data = dat,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
  refresh = 1             # no progress shown
)

dashboard(m.os)

post <- extract.samples(m.os) 

png("output/figures/OS.png", units="in", width=6, height=8, res=300)
# Set plot layout
layout(mat = matrix(c(1, 2,0, 0), 
                    nrow = 2, 
                    ncol = 2),
       heights = c(3, 1),    # Heights of the two rows
       widths = c(3, 0))     # Widths of the two columns
par(mar = c(5, 4, 4, 2)+0.1)
plot(NULL, bty= "L", ylim = c(0,1), xlim = c(0,200), 
       xlab = "Months", ylab="Overall Survival",
       main=paste("Kaplan-Meier of SEF OS"))

N <- 1:200
x <- seq(0,200, by=0.01)
R_col <- c(1,4)

for(i in 1:1) {
  for(n in 1:2) {
    lambda <- 1/exp(post$a[,i,n])
    hdpi_lambda <- HPDI(lambda)
    lambda_m <- mean(lambda)
    y1 <- dexp(x,hdpi_lambda[2])/hdpi_lambda[2]
    y2 <- dexp(x,hdpi_lambda[1])/hdpi_lambda[1]
    polygon(c(x,rev(x)),c(y2,rev(y1)),col=col.alpha(R_col[n]), border = NA)
    lines(dexp(N,lambda_m)/lambda_m, col = R_col[n], lty = n)}}
legend( "bottomleft", legend = c("Canonical","YAP1-KAMT2A"), lty = 1:2, col = c(1, 4))

med_s <- post$a
med_s_diff <- matrix(NA, nrow = 4000, ncol = 1)
for(i in 1)  med_s_diff[,i] <-  exp(med_s[,i,2]) / exp(med_s[,i,1])


par(mar = c(5, 4, 4, 2)+0.1)
plot(NULL, xlim=c(0,8), ylim = c(0,0.5),
     xlab = 'HR of Canonical Fusion', ylab = 'Density', #yaxt = 'null', 
     main = 'Hazard Ratio')
abline(v=1, lty = 2)

dens(med_s_diff, xlim = c(0,9), add = TRUE, col = 4, lwd = 3)
dev.off()

