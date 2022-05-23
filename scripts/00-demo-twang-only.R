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
#  - To run this locally with the data in twang set demo <- FALSE

demo <- TRUE

if(demo){
  AOD = read.csv("data/aod_big.csv")
  AOD <- subset(AOD, trtvar %in% c("ATM", "EAT"))
  
  # drop NA
  AOD = na.omit(AOD)
  
  # change race to a factor
  AOD$race4g = as.factor(AOD$race4g)
  
  aod_frmla <- atm ~ age + female + race4g + sfs + sps +
    sds + ias + ces + eps + imds + bcs +
    prmhtx 
  aod_frmla <- (atm ~ age + female + race4g + sfs + sps +
                  sds + ias + ces + eps + imds + bcs +
                  prmhtx)^2
  
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

ps.tmp <- ps(aod_frmla, data=AOD, 
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
if(demo){
  plot(ps.atm$gbm.obj, i.var="ias",
       n.trees=ps.atm$desc$es.max.ATT$n.trees)
} else{
  plot(ps.atm$gbm.obj, i.var="crimjust",
       n.trees=ps.atm$desc$es.max.ATT$n.trees)
}


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
AOD$w <- get.weights(ps.atm, estimand = "ATT", stop.method = "es.max")

# use svyglm to incorporate weights
design.ps <- svydesign(ids = ~1, weights = ~w, data = AOD)
glm1 <- svyglm(sfs8p12 ~ atm, design = design.ps)
summary(glm1)


