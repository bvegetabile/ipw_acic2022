######################################################
## Example code from American Causal Inference Conference Workshop
## 
## Last update: 05/23/2022
##
## Note: Code may not replicate results in slides
##       due to differences in underlying data,
##       and/or changes to underlying algorithms


## clear worksapce
rm(list=ls())

## load packages
# install.packages(c("twang","CBPS","ebal"))
library(twang)
library(CBPS)
library(ebal)
library(cobalt)
library(entbal)


## This is only a subset of the data included in the session 5 slides.
#  - Note a subset of this data is available in twang: data("AOD")

# data("AOD", package = 'twang')
AOD = read.csv("data/aod_big.csv")
AOD <- subset(AOD, trtvar %in% c("ATM", "EAT"))

# drop NA
AOD = na.omit(AOD)

# change race to a factor
AOD$race4g = as.factor(AOD$race4g)

## fit gbm and extract propensity score weights
ps.atm <- ps(atm ~ age + female + race4g + sfs + sps +
               sds + ias + ces + eps + imds + bcs +
               prmhtx , data=AOD, 
             estimand="ATT",
             n.trees=10000,
             shrinkage = 0.001,
             stop.method="es.max")

ps.tmp <- ps(atm ~ age + female + race4g + sfs + sps +
               sds + ias + ces + eps + imds + bcs +
               prmhtx , data=AOD, 
             estimand="ATT",
             n.trees=1000,
             shrinkage = 0.1,
             stop.method="es.max")

# diagnostic plot - number of iterations
plot(ps.atm, plots = "optimize")
plot(ps.tmp, plots = "optimize")

# overview of balance
round(summary(ps.atm), 4)

# marginal relationship from between a covariate and the treatment
plot(ps.atm$gbm.obj, i.var="ias",
     n.trees=ps.atm$desc$es.max.ATT$n.trees)

# relative influence of each covariate
summary(ps.atm$gbm.obj,
        n.trees = ps.atm$desc$es.max.ATT$n.trees, 
        las = 2)

# summarize balance for each covariate
bal.table(ps.atm)
do.call(cbind,lapply(bal.table(ps.atm) , function(x) x[,"std.eff.sz",drop=F] ))


# other diagnostic plots
# boxplot of propensity score by treatment group 
plot(ps.atm, plots = 2) 

# Let's also plot weights by propensity scores 
plot(unlist(ps.atm$ps), unlist(ps.atm$w), pch = 19, col = rgb(0,0,0,0.5))

# standardized effect size before and after weighting
plot(ps.atm, plots = 3) 

# QQ plot of t-test pvalues
plot(ps.atm, plots = 4) 

# QQ plot of KS p-values
plot(ps.atm, plots = 5)


##################################
## perform the outcome analysis ##

# use survey package - install if not already installed
# install.packages("survey")
library(survey)

# extract the propensity score weights
AOD$w <- get.weights(ps.atm, estimand = "ATT",
                     stop.method = "es.max")

# use svyglm to incorporate weights
design.ps <- svydesign(ids = ~1, weights = ~w, data = AOD)
glm1 <- svyglm(sfs8p12 ~ atm, design = design.ps)
summary(glm1)




############################################
## use logistic regression instead of GBM ##

## fit logistic regression and derive propensity score weights for ATT
m.logit <- glm(atm ~ age + female + race4g + sfs + sps +
                 sds + ias + ces + eps + imds + bcs +
                 prmhtx, family=binomial , data=AOD)



# for ATT, the weights are 1 for treated
# and the odds of treatment, or exp(log-odds), for control
AOD$w.logit <- ifelse(AOD$atm==1,  1, exp(predict(m.logit)))

# dxwts can assess balance for any weights
dx.logit <- dx.wts(x=AOD$w.logit, data=AOD,
                   estimand="ATT",
                   vars=c("age","female","race4g","sfs","sps","sds","ias","ces","eps","imds","bcs","prmhtx"),
                   treat.var="atm",
                   x.as.weights=TRUE)

# summary of balance
dx.logit$summary

# full balance tables
bal.table(dx.logit)

# estimate effect using logistic based weights
design.logit <- svydesign(ids = ~1, weights = ~w.logit, data = AOD)
glm6 <- svyglm(sfs8p12 ~ atm, design = design.logit) 
summary(glm6)

################################################################################
# CBPS Example

plog <- glm(atm ~ age + female + race4g + sfs + sps + sds + ias + ces + eps + imds + bcs +
              prmhtx , family=binomial , data=AOD)
AOD$ps1 <- ifelse(AOD$atm==1, 1, exp(predict(plog)))


b1 <- dx.wts(AOD$ps1, data=AOD,
             vars=c("age", "female", "race4g", "sfs", "sps", "sds", "ias", "ces", "eps", "imds", "bcs", "prmhtx"),
             treat.var="atm", estimand="ATT", x.as.weights=TRUE, sampw=NULL, perm.test.iters=0)

pcbps <- CBPS(atm ~ age + female + race4g + sfs + sps + sds + ias + ces + eps + imds +
                bcs + prmhtx, 
              method = 'over',
              data=AOD , ATT=TRUE)
AOD$ps2 <- pcbps$weights


b2 <- dx.wts(AOD$ps2, data=AOD ,
             vars=c("age", "female", "race4g", "sfs", "sps", "sds", "ias", "ces", "eps", "imds", "bcs", "prmhtx"),
             treat.var="atm", estimand="ATT", x.as.weights=TRUE, sampw=NULL, perm.test.iters=0)

pgbm <- ps(atm ~
             age + female + race4g + sfs + sps + sds + ias + ces + eps + imds + bcs + prmhtx , 
           data=AOD ,
           estimand="ATT", n.trees=10000, stop.method="es.max")
AOD$ps3 <- unlist(pgbm$w)

b3 <- dx.wts(AOD$ps3, data=AOD ,
             vars=c("age", "female", "race4g", "sfs", "sps", "sds", "ias", "ces", "eps", "imds", "bcs", "prmhtx"),
             treat.var="atm", estimand="ATT", x.as.weights=TRUE, sampw=NULL, perm.test.iters=0)


rbind(b1$summary[2,],b2$summary[2,],b3$summary[2,])


# Compare ATT Estiamtes
library(survey)
d1 <- svydesign(id=~1, weights=~ps1, data=AOD) 
d2 <- svydesign(id=~1, weights=~ps2, data=AOD) 
d3 <- svydesign(id=~1, weights=~ps3, data=AOD)
f1 <- svyglm(sfs8p12 ~ atm, design=d1) 
f2 <- svyglm(sfs8p12 ~ atm, design=d2) 
f3 <- svyglm(sfs8p12 ~ atm, design=d3)
res <- list(logit=summary(f1)$coef, 
            cbps=summary(f2)$coef,
            gbm=summary(f3)$coef)
print(res)


# Entropy balancing 
library(entbal)
library(ebal)

ebpars <- ebpars_default_binary(estimand = 'ATT')

system.time(entbal_wts <- entbal(atm ~ age + female + race4g + sfs + sps + sds + ias + ces + eps + imds + bcs +
                                   prmhtx, data=AOD, eb_pars = ebpars))

AOD$ps4 <- entbal_wts$wts

b4 <- dx.wts(AOD$ps4, data=AOD,
             vars=c("age", "female", "race4g", "sfs", "sps", "sds", "ias", "ces", "eps", "imds", "bcs", "prmhtx"),
             treat.var="atm", estimand="ATT", x.as.weights=TRUE, sampw=NULL, perm.test.iters=0)

tmpx <- model.matrix(~ -1 
                     + age + female + race4g + sfs + sps 
                     + sds + ias + ces + eps + imds + bcs 
                     + prmhtx, 
                     data = AOD)
colnames(tmpx)
tmpx <- tmpx[,-3]
eb_wts <- ebalance(AOD$atm, tmpx)

AOD$ps5 <- 1 
AOD$ps5[AOD$atm==0] <- eb_wts$w

b5 <- dx.wts(AOD$ps5, data=AOD,
             vars=c("age", "female", "race4g", "sfs", "sps", "sds", "ias", "ces", "eps", "imds", "bcs", "prmhtx"),
             treat.var="atm", estimand="ATT", x.as.weights=TRUE, sampw=NULL, perm.test.iters=0)


# CBPS exact

system.time(pcbps_exact <- CBPS(atm ~ age + female + race4g + sfs + sps + sds + ias + ces + eps + imds +
                                  bcs + prmhtx, 
                                data=AOD, 
                                method = 'exact', 
                                ATT=TRUE))
AOD$ps6 <- pcbps_exact$weights
b6 <- dx.wts(AOD$ps6, data=AOD,
             vars=c("age", "female", "race4g", "sfs", "sps", "sds", "ias", "ces", "eps", "imds", "bcs", "prmhtx"),
             treat.var="atm", estimand="ATT", x.as.weights=TRUE, sampw=NULL, perm.test.iters=0)

bal_res <- rbind(b1$summary[2,],
                 b2$summary[2,],
                 b3$summary[2,],
                 b4$summary[2,],
                 b5$summary[2,], 
                 b6$summary[2,])
bal_res$type <- c('GLM', 'CBPS', 'GBM', 'ENTBAL', 'EBAL', "CBPS2")
bal_res

d1 <- svydesign(id=~1, weights=~ps1, data=AOD) 
d2 <- svydesign(id=~1, weights=~ps2, data=AOD) 
d3 <- svydesign(id=~1, weights=~ps3, data=AOD)
d4 <- svydesign(id=~1, weights=~ps4, data=AOD)
d5 <- svydesign(id=~1, weights=~ps5, data=AOD)
d6 <- svydesign(id=~1, weights=~ps6, data=AOD)

f1 <- svyglm(sfs8p12 ~ atm, design=d1) 
f2 <- svyglm(sfs8p12 ~ atm, design=d2) 
f3 <- svyglm(sfs8p12 ~ atm, design=d3)
f4 <- svyglm(sfs8p12 ~ atm, design=d4) 
f5 <- svyglm(sfs8p12 ~ atm, design=d5)
f6 <- svyglm(sfs8p12 ~ atm, design=d6)

res <- list(logit=summary(f1)$coef, 
            cbps=summary(f2)$coef,
            gbm=summary(f3)$coef,
            entbal=summary(f4)$coef, 
            ebal=summary(f5)$coef,
            cbps_exact = summary(f6)$coef)
print(res)
