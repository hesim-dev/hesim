`hesim` is an R package for health economic simulation modeling and decision analysis. The package can help facilitate computationally intensive simulation modeling and be used to analyze the output of simulation models. Current functionality includes:

* Individualized cost-effectiveness analysis
* Random sampling from probability distributions commonly used in health-economic simulation modeling

To ensure that simulations can be run (and analyzed) in a reasonable amount of time, most functions are written in C++ using `Rcpp` and data manipulations are performed using the `data.table` package. `hesim` is therefore well suited for individual patient simulation, probabilistic sensitivity analysis, and quantifying structural uncertainty.

# Documentation
Documentation and package vignettes are available on the package [website](https://innovationvalueinitiative.github.io/hesim/).

# Installation
`hesim` can be installed from GitHub using `devtools`:

```r
install.packages("devtools")
library(devtools)
devtools::install_github("InnovationValueInitiative/hesim")
```