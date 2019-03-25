## hesim 0.2.0
### Highlights
`hesim` now provides a general framework for integrating statistical models with economic evaluation. Users 
build a decision model by specifying a model structure, which consists of a set of statistical models for disease progression, utilities, and costs. Each statistical model is used to simulate outcomes as a function of estimated parameters and input data. N-state partitioned survival models (PSMs) and individual-level continuous time state transition models (iCTSTMs) are now supported. 

### API changes
* The argument `sim` was renamed `sample` in `icea()`, `icea_pw()`, and `incr_effect()`.
* Custom functions and variables are no longer supported in `icea` and `icea_pw()`.

### New features
Economic models---which combine the disease, utility, and cost models---are constructed with the following classes:
* `Psm()` for PSMs
* `IndivCtstm()` for iCTSTMs

Disease models are constructed using the classes:
* `PsmCurves()` to simulate survival curves for each endpoint of interest
* `IndivCtstmTrans()` to simulate health state transitions with a iCTSTM

Utility and cost models are constructed with the `StateVals()` class. 

The economic models are used to simulate disease progression (`$sim_disease()`, `sim_stateprobs()`), quality-adjusted life-years (QALYs) (`sim_qalys()`), and costs (`sim_costs()`). Parameter uncertainty is propagated to model outcomes using probabilistic sensitivity analysis. Summaries of the simulated costs and QALYs are used to perform model-based cost-effectiveness analyses (CEAs) and represent decision uncertainty with `icea.ce()` and `icea_pw.ce()`.

## hesim 0.1.0
The initial CRAN submission containing support for CEA but not for model development. Decision uncertainty is represented using cost-effectiveness planes, cost-effectiveness acceptability curves, cost-effectiveness acceptability frontiers, and the expected value of perfect information. CEAs by subgroup (i.e., individualized CEAs) are performed with `icea()` and `icea_pw()`.