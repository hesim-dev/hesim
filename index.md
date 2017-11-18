[![Travis-CI Build Status](https://travis-ci.org/InnovationValueInitiative/hesim.svg?branch=master)](https://travis-ci.org/InnovationValueInitiative/hesim)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/hesim)](https://cran.r-project.org/package=hesim)

# Overview
`hesim` is an R package for health economic simulation modeling and decision analysis. The package can help facilitate computationally intensive simulation modeling and be used to analyze the output of simulation models. Current functionality includes:

* Individualized cost-effectiveness analysis
* Random sampling for probabilistic sensitivity analysis (PSA) and individual patient simulation (IPS)

To ensure that simulations can be run (and analyzed) in a reasonable amount of time, most functions are written in C++ using `Rcpp` and data manipulations are performed using the `data.table` package. `hesim` is therefore well suited for IPS, PSA, and quantifying structural uncertainty.

# Installation
`hesim` can be installed from GitHub using `devtools`:

```r
install.packages("devtools")
library(devtools)
devtools::install_github("InnovationValueInitiative/hesim")
```

It can then be loaded into `R`:

```r
library(hesim)
```