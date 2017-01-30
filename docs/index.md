# hesim

## About

Simulation models are commonly used in health economics for cost-effectiveness analysis, modeling disease burden, and other purposes. The `hesim` package is designed for computationally intensive simulation and has the following three goals:

* Help users conduct decision analysis and probabilistic sensitivity analysis at the subgroup level.
* Run very fast patient-level and cohort simulation models.
* Analyze the outcomes of both discrete and continuous time simulation models.

Since both PSA and patient-level models are computationally intensive, key functions are written in C++ using `Rcpp` and data manipulations are performed using the `data.table` package for speed. `hesim` is a natural complement to existing packages for survival analysis and multi-state modeling (`survival`, `flexsurv`, `msm`, `mstate`), which are needed for estimating the parameters of disease models.

`hesim` is under development. It currently provides a number of functions for decision analysis and PSA at the subgroup level, but new functionality for simulating and analyzing patient-level and cohort health-economic models will be added shortly.

## Installation
`hesim` can be installed from GitHub using `devtools`

```r
install.packages("devtools")
library(devtools)
devtools::install_github("dincerti/hesim")
```

The package can then be loaded into `R`

```r
library(hesim)
```