# About

Simulation models are commonly used in health economics for cost-effectiveness analysis, modeling disease burden, and other purposes. The `hesim` package is designed for computationally intensive simulation and has the following three goals:

* Run very fast patient-level and cohort simulation models.
* Analyze the outcomes of both discrete and continuous time simulation models.
* Help users perform Bayesian cost-effectiveness analysis (e.g. probabilistic sensitivity analyses) at the subgroup level.

Since both PSA and patient-level models are slow using base R, key functions are written in C++ using `Rcpp` and data manipulations are performed using the `data.table` package for speed. `hesim` is a natural complement to existing packages for survival analysis and multi-state modeling (`survival`, `flexsurv`, `msm`, `mstate`), which are needed for estimating the parameters of disease models.

`hesim` is under development. It currently provides a number of functions for cost-effectiveness analysis at the subgroup level, but new functionality for simulating and analyzing patient-level and cohort health-economic models is being added.

# Documentation

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