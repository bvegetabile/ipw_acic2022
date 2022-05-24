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
  AOD$trtvar <- AOD$treat
}

## fit gbm and extract propensity score weights
ps.atm <- ps(aod_frmla, data=AOD, 
             estimand="ATT",
             n.trees=10000,
             shrinkage = 0.001,
             stop.method="es.max")

# extract the propensity score weights
AOD$w <- get.weights(ps.atm, estimand = "ATT", stop.method = "es.max")
AOD$ps <- unlist(ps.atm$ps)

# ------------------------------------------------------------------------------
# Balance tables from COBALT 

# Simple balance table
cobalt::bal.tab(aod_frmla, data=AOD, weights = 'w', un = TRUE)

# Balance table with mean balance
cobalt::bal.tab(aod_frmla, data=AOD, weights = 'w', 
                disp = c("means"), un = TRUE, 
                stats = c("mean.diffs"))

# Full balance table
cobalt::bal.tab(aod_frmla, data=AOD, weights = 'w', 
                disp = c("means", "sds"), un = TRUE, 
                stats = c("mean.diffs", "variance.ratios"))

# ------------------------------------------------------------------------------
# Balance plots

bal.plot(trtvar ~ ps, data=AOD, weights = AOD$w,
         var.name = "ps", which = "both",
         type = "histogram", mirror = TRUE)
if(demo){
  cobalt::bal.plot(aod_frmla, data=AOD, weights = AOD$w, which = 'both', var.name = 'sfs')  
} else{
  cobalt::bal.plot(aod_frmla, data=AOD, weights = AOD$w, which = 'both', var.name = 'crimjust')
}

cobalt::love.plot(aod_frmla, data=AOD, weights = AOD$w, 
                  stats = c("mean.diffs", "ks.statistics"), 
                  threshold = c(m = .1, ks = .1), 
                  binary = "std",
                  shapes = c("circle", "triangle"),
                  colors = c("red", "blue"))

# ------------------------------------------------------------------------------

