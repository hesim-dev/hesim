## hesim 0.1.0.9000
The current development version of `hesim`.

### Highlights
`hesim` now provides a general framework for integrating statistical models with economic evaluation. Users can 
build a decision model by specifying a model structure, which consists of a set of statistical models for disease progression, utility values, and costs. Each statistical model is used to simulate outcomes as a function of estimated parameters. N-state partitioned survival analysis is now supported. 

### New features
Fitted statistical models
* There are functions to create objects that store collections of fitted statistical models or `formula` objects. These include `partsurvfit()`, `formula_list()`, and `flexsurvreg_list()`.  

Parameters
* Functions prefixed by `params_` create objects storing samples of parameters of fitted statistical models for probabilistic sensitivity analysis.
* `form_params()` is a generic function for constructing parameter objects from a fitted statistical model or a `formula` object. Parameters can be sampled using Monte Carlo multivariate normal approximations or via bootstrapping. 
* Current support for flexible survival modeling (`params_surv()`, `params_surv_list()`) and linear regression (`params_lm()`). Splines and parametric distributions (exponential, Weibull, Gompertz, gamma, lognormal log-logistic, generalized gamma) are supported for survival modeling.

Input data
* `hesim_data()` creates an object of class `hesim_data` for storing a collection of data tables or data frames for simulation modeling.
* `expand_hesim_data()` combines some or all of the data tables or data frames in `hesim_data()` into a single long dataset.
* `input_data()` creates an object of class `input_data`, which contains data for predicting or simulating values with a statistical model.
* `form_input_data()` creates an object of class `input_data` from a fitted statistical model or a `formula` object.

Health state values
* The R6 class `StateVals` simulates the costs or utilities associated with health states.
* `form_StateVals()` creates a `StateVals` object from fitted statistical models or `formula` objects.

Partitioned survival models
* The R6 class `Psm` simulates outcomes from an N-state partitioned survival model. 
* A `Psm` object is instantiated with a set of survival models (the R6 class `PsmCurves`) and models for costs and utility (the R6 class `StateVals`).
* `form_PsmCurves()` creates a `PsmCurves` object from fitted statistical models or `formula` objects.

Datasets
* `psm4_exdata` provides example datasets for parameterizing a partitioned survival model. 

## hesim 0.1.0
The initial CRAN submission.