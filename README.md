# Parametric G-computation for Compatible Indirect Treatment Comparisons with Limited Individual Patient Data: Code

### Antonio Remiro-Azócar, Anna Heath, Gianluca Baio
### *remiroantonio@gmail.com*
### *2021*

This repository contains the R code used for my paper [Parametric G-computation for Compatible Indirect Treatment Comparisons with Limited Individual Patient Data][1], co-authored with [Prof. Gianluca Baio][2] and [Prof. Anna Heath][3]. 

## Utilizing the Scripts

In order to use this repository, the user must first download a copy to their local machine. The user must set the working directory to the location where the download was made. To run the pipeline, the user should then open and run the scripts in the following order:

|          Script           | Explanation                                                  |
| :-----------------------: | ------------------------------------------------------------ |
|       `gen_data.R`        | Specifies the settings of the simulation study and saves them in `"./binary_settings.RData"`. Generates the data for the simulation study according to the settings (saving the data to the `"./Data/"` subdirectory) |
| `population_adjustment.R` | Performs the population-adjusted indirect comparison methods on the simulated data (saving the point estimates and variances to the `"./Results/"` subdirectory) |
|       `analysis.R`        | Processes the results of the simulation study and computes and graphs the relevant performance metrics (the analyses are saved to the `"./Analysis/"` subdirectory) |

In addition, the `functions.R` script contains a user-defined function for weight estimation in MAIC and functions to evaluate the performance measures of interest. The file `./Analysis/scenarios.csv` records the parameter values or settings for each scenario and the key performance measures/results associated with each (as presented in the paper). 

The `./Example` subdirectory features example `R` code implementing matching-adjusted indirect comparison (MAIC), conventional simulated treatment comparison (STC), maximum-likelihood parametric G-computation and Bayesian parametric G-computation, as per Appendix D in the Supplementary Material of the article. 

The `doSNOW` package is used to parallelize the performance of the methods, distributing the tasks to different cores of the computer. 

The code presented here was prepared in `RStudio` using `R` version `3.6.3` in a Windows architecture, with 64-bit operating system. The following packages and version were used:

* `boot 1.3-24` required for use of the non-parametric bootstrap in MAIC and maximum-likelihood parametric G-computation
* `copula 0.999.20` simulates covariates from a multivariate Gaussian copula for the covariate simulation step of maximum-likelihood parametric G-computation and Bayesian parametric G-computation. 
* `doSNOW 1.0.18` used in combination with `foreach()` to start up local clusters that distribute parallel tasks to different cores
* `dplyr 1.0.2` for data manipulation
* `ggplot2 3.3.2` to plot the simulation study results (Figures 2-4 in the article)
* `gridExtra 2.3` to plot the simulation study results (Figures 2-4 in the article)
* `MASS 7.3-51.5` to simulate correlated covariates in the data-generating process for the simulation study, drawing from a multivariate normal distribution 
* `parallel 3.6.3` to detect the number of CPU cores
* `rstanarm 2.21.1` for fitting outcome regressions and drawing outcomes from their posterior predictive distribution in Bayesian parametric G-computation.  

[1]: https://doi.org/10.1002/jrsm.1565
[2]: http://www.statistica.it/gianluca/
[3]: https://sites.google.com/site/annaheathstats/
