# Example R code implementing maximum-likelihood parametric G-computation
# with Cox regression as the outcome model and survival outcomes

library("survival") # to fit Cox proportional hazards regression
library("copula") # for simulating BC covariates from Gaussian coupla
library("boot") # for non-parametric bootstrap

set.seed(555) # set seed for reproducibility

# setwd("C:/Users/Antonio/Desktop/Gcomp_indirect_comparisons_simstudy/Example/Survival")

### Parametric G-computation with maximum-likelihood estimation ###

rm(list=ls())
AC.IPD <- read.csv("AC_IPD_survival.csv") # load AC patient-level data
BC.ALD <- read.csv("BC_ALD_survival.csv") # load BC aggregate-level data

# matrix of pairwise correlations between IPD covariates
rho <- cor(AC.IPD[,c("X1","X2","X3","X4")])
#  covariate simulation for BC trial using copula package
cop <- normalCopula(param=c(rho[1,2],rho[1,3],rho[1,4],rho[2,3],
                            rho[2,4],rho[3,4]),
                    dim=4, dispstr="un") # AC IPD pairwise correlations
# sample covariates from approximate joint distribution using copula
mvd <- mvdc(copula=cop, margins=c("norm", "norm", # Gaussian marginals
                                  "norm", "norm"),
            # BC covariate means and standard deviations
            paramMargins=list(list(mean=BC.ALD$mean.X1, sd=BC.ALD$sd.X1),
                              list(mean=BC.ALD$mean.X2, sd=BC.ALD$sd.X2),
                              list(mean=BC.ALD$mean.X3, sd=BC.ALD$sd.X3),
                              list(mean=BC.ALD$mean.X4, sd=BC.ALD$sd.X4)))
# simulated BC pseudo-population of size 1000
x_star <- as.data.frame(rMvdc(1000, mvd))
colnames(x_star) <- c("X1", "X2", "X3", "X4")

# function to be resampled by non-parametric bootstrap
gcomp.ml <- function(data, indices) {
  dat = data[indices,]
  # outcome Cox regression model fitted to IPD using maximum likelihood
  outcome.model <- coxph(Surv(time, status)~trt*X1+trt*X2+X3+X4, data=dat)
  # event time selected for unit 50 (random selection) 
  unit.time <- 50
  # estimated cumulative baseline hazard
  hat.H0 <- basehaz(outcome.model)[unit.time,1] 
  # counterfactual datasets (two hypothetical worlds)
  data.trtA <- data.trtC <- x_star
  # intervene on treatment while keeping set covariates fixed
  data.trtA$trt <- 1 # dataset where everyone receives treatment A
  data.trtC$trt <- 0 # dataset where all observations receive C
  # linear predictor where everyone receives treatment A
  LP.A <- with(outcome.model, x_star$X1*(coefficients["X1"] + coefficients["trt:X1"]) + 
                 x_star$X2*(coefficients["X2"] + coefficients["trt:X2"]) + 
                 x_star$X3*coefficients["X3"] + x_star$X4*coefficients["X4"] +
                 coefficients["trt"])
  # linear predictor where all observations receive treatment C
  LP.C <- with(outcome.model, x_star$X1*coefficients["X1"] + x_star$X2*coefficients["X2"] +
                 x_star$X3*coefficients["X3"] + x_star$X4*coefficients["X4"])
  # predict individual survival probabilities, conditional on treatment/covariates
  hat.S.A.i <- exp(-hat.H0)^exp(LP.A) 
  hat.S.C.i <- exp(-hat.H0)^exp(LP.C)
  # mean survival probability prediction under each treatment
  hat.P.A <- mean(hat.S.A.i)
  hat.P.C <- mean(hat.S.C.i)
  # estimate marginal A vs. B log hazard ratio (mean difference in expected log hazard)
  # by transforming from survival probability to linear predictor scale 
  hat.Delta.AC <- log(-log(hat.P.A)) - log(-log(hat.P.C))
  return(hat.Delta.AC)
} 

# non-parametric bootstrap with 1000 resamples (ignore warnings)
boot.object <- boot::boot(data=AC.IPD, statistic=gcomp.ml, R=1000)
# bootstrap mean of marginal A vs. C treatment effect estimate
hat.Delta.AC <- mean(boot.object$t)
# bootstrap variance of A vs. C treatment effect estimate
hat.var.Delta.AC <- var(boot.object$t)
# marginal log hazard ratio for B vs. C reported in BC article
hat.Delta.BC <- BC.ALD$logHR_B
# variance of B vs. C in aggregate outcomes in published article
hat.var.Delta.BC <- BC.ALD$var_logHR_B
# marginal treatment effect for A vs. B
hat.Delta.AB <- hat.Delta.AC - hat.Delta.BC
# variance for A vs. B
hat.var.Delta.AB <- hat.var.Delta.AC + hat.var.Delta.BC
# construct Wald-type normal distribution-based confidence interval
uci.Delta.AB <- hat.Delta.AB + qnorm(0.975)*sqrt(hat.var.Delta.AB)
lci.Delta.AB <- hat.Delta.AB + qnorm(0.025)*sqrt(hat.var.Delta.AB)