# Overview
`hesim` is an R package for health economic simulation modeling and decision analysis. The package has two primary aims: 

* Analyze the output of discrete and continuous time simulation models.
* Provide functions to randomly sample from probability distributions commonly used in health-economic simulation modeling.

To ensure that simulations can be run (and analyzed) in a reasonable amount of time, most functions are written in C++ using `Rcpp` and data manipulations are performed using the `data.table` package. `hesim` is therefore well suited for individual patient simulation, probabilistic sensitivity analysis, and quantifying structural uncertainty.

# Installation
`hesim` can be installed from GitHub using `devtools`:

```r
install.packages("devtools")
library(devtools)
devtools::install_github("dincerti/hesim")
```

It can then be loaded into `R`:

```r
library(hesim)
```