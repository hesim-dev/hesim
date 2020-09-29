[![Travis-CI Build Status](https://travis-ci.org/hesim-dev/hesim.svg?branch=master)](https://travis-ci.org/hesim-dev/hesim)
[![Coverage Status](https://codecov.io/gh/hesim-dev/hesim/branch/master/graph/badge.svg)](https://codecov.io/gh/hesim-dev/hesim)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/hesim)](https://cran.r-project.org/package=hesim)

# Overview
`hesim` is a modular and computationally efficient R package for health economic simulation modeling and decision analysis that provides a general framework for integrating statistical analyses with economic evaluation. The package supports cohort discrete time state transition models (DTSTMs), N-state partitioned survival models (PSMs), and individual-level continuous time state transition models (CTSTMs), encompassing both Markov (time-homogeneous and time-inhomogeneous) and semi-Markov processes. It heavily utilizes `Rcpp` and `data.table`, making individual-patient simulation, probabilistic sensitivity analysis (PSA), and incorporation of patient heterogeneity fast. 

Features of the current version can be summarized as follows:

* Cohort DTSTMs, N-state PSMs, and individual-level CTSTMs that encompass Markov and semi-Markov processes
* Options to build models via mathematical expressions using nonstandard evaluation or directly from fitted statistical models
* Parameter estimates from either an R based model or from an external source
* Convenience functions for sampling model parameters from parametric distributions or via bootstrapping
* Parameter uncertainty propagated with PSA
* Modeling patient heterogeneity 
* Performing cost-effectiveness analyses and representing decision uncertainty from PSAs
* Simulation code written in `C++` to boost performance

# Installation
```r
# Install from CRAN:
install.packages("hesim")

# Install the most up to date development version from GitHub:
# install.packages("devtools")
devtools::install_github("hesim-dev/hesim")
```