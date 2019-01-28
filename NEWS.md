## hesim 0.1.0.9000
The current development version of `hesim`.

### Highlights
`hesim` now provides a general framework for integrating statistical models with economic evaluation. Users can 
build a decision model by specifying a model structure, which consists of a set of statistical models for disease progression, utilities, and costs. Each statistical model is used to simulate outcomes as a function of estimated parameters and input data. N-state partitioned survival models (PSMs) and individual-level continuous time state transition models (iCTSTMs) are now supported. 

### API changes
* The argument `sim` was renamed `sample` in `icea()`, `icea_pw()`, and `incr_effect()`.
* Custom functions and variables are no longer supported in `icea` and `icea_pw()`.

### New features
Economic models---which combine the disease, utility, and cost models---can be constructed with the following classes:
* `Psm()` for PSMs
* `IndivCtstm()` for iCTSTMs

Disease models are constructed using the classes:
* `PsmCurves()` to simulate survival curves for each endpoint of interest
* `IndivCtstmTrans()` to simulate health state transitions with a iCTSTM

Utility and cost models are constructed with the `StateVals` class. 

See the [introduction](file:///Users/devin.incerti/IVI/hesim/docs/articles/intro.html) for more details.

## hesim 0.1.0
The initial CRAN submission containing support for individualized cost-effectiveness analysis but not for model development.