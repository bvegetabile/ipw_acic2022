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
library(viridis)
library(survey)

## This is only a subset of the data included in the session 5 slides.
#  - Note a subset of this data is available in twang: data("AOD")
#  - To run this locally with the data in twang set demo <- FALSE

demo <- TRUE

if(demo){
  AOD = read.csv("data/aod_big.csv")
  AOD <- subset(AOD, trtvar %in% c("ATM", "EAT"))
  
  # drop NA
  AOD = na.omit(AOD)
  
  # change race to a factor
  AOD$race4g = as.factor(AOD$race4g)
  
  aod_frmla <- atm ~ age + female + race4g + sfs + sps + sds + ias + ces + eps + imds + bcs + prmhtx 
  aod_frmla2 <- atm ~ (age + female + race4g + sfs + sps + sds + ias + ces + eps + imds + bcs + prmhtx)^2
  
  bal_vars=c("age","female","race4g","sfs","sps","sds","ias","ces","eps","imds","bcs","prmhtx")
} else {
  data("AOD", package = 'twang')
  AOD <- AOD[AOD$treat!= 'scy', ]
  AOD$atm <- ifelse(AOD$treat == 'metcbt5', 1, 0)
  aod_frmla <- atm ~ illact + crimjust + subprob + subdep + white
  aod_frmla2 <- atm ~ (illact + crimjust + subprob + subdep + white)^2
  bal_vars=c('illact', 'crimjust', 'subprob', 'subdep', 'white')
  AOD$sfs8p12 <- AOD$suf12
}




## fit gbm and extract propensity score weights
ps.atm <- ps(aod_frmla, data=AOD, 
             estimand="ATT",
             n.trees=10000,
             shrinkage = 0.001,
             stop.method="es.max")

# extract the propensity score weights
AOD$w_twang <- get.weights(ps.atm, estimand = "ATT", stop.method = "es.max")

# Let's use cobalt's love plots for the rest of the examples
cobalt::love.plot(aod_frmla, data=AOD, weights = AOD$w_twang, 
                  stats = c("mean.diffs", "ks.statistics"), 
                  threshold = c(m = .1, ks = .1), 
                  binary = "std",
                  shapes = c("circle", "triangle"),
                  colors = c("red", "blue"))


############################################
## use logistic regression instead of GBM ##

## fit logistic regression and derive propensity score weights for ATT
m.logit <- glm(aod_frmla, family=binomial , data=AOD)

# for ATT, the weights are 1 for treated
# and the odds of treatment, or exp(log-odds), for control
AOD$w_logit <- ifelse(AOD$atm==1,  1, exp(predict(m.logit)))

# dxwts can assess balance for any weights
dx.logit <- dx.wts(x=AOD$w_logit, data=AOD,
                   estimand="ATT",
                   vars=bal_vars,
                   treat.var="atm",
                   x.as.weights=TRUE)

# summary of balance
dx.logit$summary

# full balance tables
bal.table(dx.logit)

# estimate effect using logistic based weights
design.logit <- svydesign(ids = ~1, weights = ~w_logit, data = AOD)
glm6 <- svyglm(sfs8p12 ~ atm, design = design.logit) 
summary(glm6)

# Let's use cobalt's love plots for the rest of the examples
cobalt::love.plot(aod_frmla, data=AOD, weights = list('TWANG' = AOD$w_twang, 'GLM' = AOD$w_logit), 
                  stats = c("mean.diffs", "ks.statistics"), 
                  threshold = c(m = .1, ks = .1), 
                  binary = "std",
                  colors = viridis(3))


################################################################################
# ALL OTHER MODELS AT ONCE

# CBPS Example
pcbps <- CBPS(aod_frmla, 
              method = 'over',
              data=AOD , ATT=TRUE)
AOD$w_cbps <- pcbps$weights

# Covariate balancing higher moments
pcbps_exact <- CBPS(aod_frmla2, 
                    data=AOD, 
                    method = 'exact', 
                    ATT=TRUE)
AOD$w_cbps_2 <- pcbps_exact$weights


# Entropy balancing 
ebpars <- ebpars_default_binary(estimand = 'ATT', which_z = 1)
ebpars$n_moments <- 1
entbal_wts <- entbal(aod_frmla, data=AOD, eb_pars = ebpars)
AOD$w_eb1 <- entbal_wts$wts


# Entropy balancing 
ebpars <- ebpars_default_binary(estimand = 'ATT', which_z = 1)
entbal_wts <- entbal(aod_frmla, data=AOD, eb_pars = ebpars)
AOD$w_eb3 <- entbal_wts$wts

# Let's use cobalt's love plots for the rest of the examples
cobalt::love.plot(aod_frmla, data=AOD, 
                  weights = list('TWANG' = AOD$w_twang, 
                                 'GLM' = AOD$w_logit,
                                 'CBPS-OVER' = AOD$w_cbps,
                                 'CBPS-EXACT-2'=AOD$w_cbps_2,
                                 'EB-1mom' = AOD$w_eb1,
                                 'EB-3mom' = AOD$w_eb3), 
                  stats = c("mean.diffs", 'variance.ratio', "ks.statistics"), 
                  shapes = c(15:20,15),
                  threshold = c(m = .1, ks = .1), 
                  binary = "std",
                  colors = viridis(n=7, alpha = 0.1))

# ------------------------------------------------------------------------------
# OUTCOME ESTIMATION 

d1 <- svydesign(id=~1, weights=~w_twang, data=AOD) 
d2 <- svydesign(id=~1, weights=~w_logit, data=AOD) 
d3 <- svydesign(id=~1, weights=~w_cbps, data=AOD)
d4 <- svydesign(id=~1, weights=~w_cbps_2, data=AOD)
d5 <- svydesign(id=~1, weights=~w_eb1, data=AOD)
d6 <- svydesign(id=~1, weights=~w_eb3, data=AOD)

f1 <- svyglm(sfs8p12 ~ atm, design=d1) 
f2 <- svyglm(sfs8p12 ~ atm, design=d2) 
f3 <- svyglm(sfs8p12 ~ atm, design=d3)
f4 <- svyglm(sfs8p12 ~ atm, design=d4) 
f5 <- svyglm(sfs8p12 ~ atm, design=d5)
f6 <- svyglm(sfs8p12 ~ atm, design=d6)

res <- list(twang=round(summary(f1)$coef,3), 
            logit=round(summary(f2)$coef,3),
            cbps=round(summary(f3)$coef,3),
            cbps2=round(summary(f4)$coef,3), 
            eb1=round(summary(f5)$coef,3),
            eb3 = round(summary(f6)$coef,3))
print(res)


