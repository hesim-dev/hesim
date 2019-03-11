[![Travis-CI Build Status](https://travis-ci.org/InnovationValueInitiative/hesim.svg?branch=master)](https://travis-ci.org/InnovationValueInitiative/hesim)
[![Coverage Status](https://codecov.io/gh/InnovationValueInitiative/hesim/branch/master/graph/badge.svg)](https://codecov.io/gh/InnovationValueInitiative/hesim)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/hesim)](https://cran.r-project.org/package=hesim)

# Overview
`hesim` is an R package for health economic simulation modeling and decision analysis that provides a general framework for integrating statistical analyses with economic evaluation. The package currently supports N-state partitioned survival models (PSMs) and individual-level continuous time state transition models (CTSTMs), as well as summarizing the output of probabilistic sensitivity analysis (PSA). It is designed for high performance simulation modeling and heavily utilizes `Rcpp` and `data.table`. `hesim` is being actively developed and subsequent versions will provide support for cohort-level state transition models. 

Features of the current version include:

* N-state PSMs and individual-level CTSTMs
* Modeling patient heterogeneity 
* Parameter estimates from either an R based model or from an external source
* Sampling parameters of a statistical model fit with R using Monte Carlo methods or bootstrapping
* Parameter uncertainty propagated with PSA
* Performing cost-effectiveness analyses and representing decision uncertainty from PSAs
* Simulation code written in `C++` to boost performance

# Installation
```r
# Install the most up to date development version from GitHub:
# install.packages("devtools")
devtools::install_github("InnovationValueInitiative/hesim")

# Install v0.1 (without economic modeling) from CRAN:
install.packages("hesim")


```