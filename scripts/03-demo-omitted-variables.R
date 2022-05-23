######################################################
## Example code from American Causal Inference Conference Workshop
## 
## Last update: 05/23/2022
##
## Note: Code may not replicate results in slides
##       due to differences in underlying data,
##       and/or changes to underlying algorithms

## Adapted from the vignette for OVtool
## https://cran.r-project.org/web/packages/OVtool/vignettes/OVtool.html


## clear worksapce
rm(list=ls())

## load packages
# install.packages(c("twang","CBPS","ebal"))
library(twang)
library(cobalt)
library(viridis)
library(survey)
library(OVtool)

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
  outcome_form <- sfs8p12 ~ atm + age + female + race4g + sfs + sps + sds + ias + ces + eps + imds + bcs + prmhtx 
  
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
# ------------------------------------------------------------------------------
# Original Outcome model
des <- svydesign(id=~1, weights=~w_twang, data=AOD) 
mod <- svyglm(outcome_form, design=des) 


# ------------------------------------------------------------------------------
# Sensitivity Analysis for TWANG
ovmod <- outcome_model(ps_object = ps.atm, 
              stop.method = "es.max", 
              data = sud,
              weights = "w_twang",
              treatment = "atm",
              outcome = "sfs8p12",
              model_covariates = bal_vars,
              estimand = "ATT")
saveRDS(ovmod, 'models/sensitivity_outcomemodel.rds')

sens <- ov_sim(ovmod, plot_covariates = bal_vars)
saveRDS(sens, 'models/sensitivity_sensitivity_output.rds')


plot.ov(sens, print_graphic = 3, col = 'color')

# ------------------------------------------------------------------------------
# Interpretation: 

# The solid black contours represent the treatment effect contour lines and the 
# red lines (sometimes dashed) represent the p-value thresholds. The key on the 
# right side of the graphic shows where various p-value cutoff lines are, 
# including p = 0.05. The blue points on the plot represent the observed 
# covariate correlations with the outcome (y-axis) and effect size associations 
# with the treatment indicator (x-axis).

# If an unobserved confounder had a relationship with the outcome and the 
# treatment indicator that was equivalent to those in the dotted region

# Then the researcher would conclude that their results would likely be 
# sensitive to an unobserved confounder at a p-value threshold of 0.05.

# If the blue points all existed in contours greater than the 0.05 p-value 
# contour, then unobserved confounders with similar associations would retain 
# the significant effect and allow the user to conclude that the results are 
# reasonably robust.
