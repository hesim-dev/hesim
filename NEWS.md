## hesim 0.4.1
Minor updates to the documentation and fixes to small problems in the `C++` code identified with the CRAN package checks.

## hesim 0.4.0
### New features

* `IndivCtstmTrans` objects can be constructed from a `params_surv_list` using `create_IndivCtstmTrans.params_surv_list()`. 

* Survival models can randomly sample from piecewise exponential and proportional hazards Weibull distributions. A `fixed` distribution has also been added so that survival times can be set to a single constant value. Random number generation from truncated versions of these distributions is also supported. Note that functionality beyond random number generation (e.g., hazard functions, cumulative hazard functions, cumulative density functions) is not yet complete. See `?params_surv`. 

* A new vignette incorporates the two bullets above and shows how a time-inhomogeneous Markov model can be simulated using individual patient simulation. 

* Disease progression (i.e., a trajectory through a multi-state model) can be simulated using the `sim_disease()` method of the `hesim::IndivCtstmTrans` class.

* A more computationally efficient approach to simulation of time-inhomogeneous Markov cohort models has been added to the corresponding vignette. This was aided by the new `tpmatrix_id()` and `tparams_transprobs.tpmatrix()` functions.

* The "Articles" on the package website have been reorganized so that they align more closely with the different types of economic models. 

### API changes

* `icea()` and `icea_pw()` have been deprecated and replaced by `cea()` and `cea_pw()`. 

### Bug fixes

* The `lys` argument for the `$sim_qalys()` method of `hesim::Psm` and `hesim::CohortDtstm` classes now works. 

* The `$sim_stateprobs()` argument for the `hesim::Psm` class now properly returns the `patient_wt` column. 

## hesim 0.3.1
Fixes a small bug in the `C++` code identified with the CRAN package checks.

## hesim 0.3.0
### Highlights
`hesim` now supports discrete time state transition models via `hesim::CohortDtstm` objects. Users can build a model by either fitting multinomial logistic regressions with `nnet::multinom()` or with a mathematical expression using `define_model()`. Furthermore, `$summarize()` methods now have a `by_grp` option to facilitate subgroup analyses. 

### New features

* The `hesim::CohortDtstm` class simulates cohort discrete time state transition models (cDTSTM). State transitions in a cDTSTM are simulated using the `hesim::CohortDtstmTrans` class, which can be constructed from a `multinom_list()` object or using `define_model()`.

* `$summarize()` methods now have a `by_grp = "TRUE"` option to facilitate subgroup analyses. If `by_grp = "FALSE"`, then estimates are aggregated across groups. A new `patient_wt` argument in the `patients` table of `hesim_data()` can be used to weight groups during the aggregation.

* `hesim::tparams` objects can now be used to store transformed parameters used to simulate outcomes such as means (i.e., `tparams_mean()`) that have already been predicted as a function of covariates.

* General cumulative hazard functions can now be simulated using a discrete time approximation where the probability of an event during each time period is simulated from a Bernoulli distribution. This is more efficient than the previous method based on a `C++` version of the `sample()` function. See the `random_method = "discrete"` option in `params_surv()`.

* `rdirichlet_mat()` has a new argument `output` so that multiple object types can be returned. 

### API changes

* The auxiliary argument `random_method = "sample"` in `params_surv()` is deprecated and `random_method = "discrete"` should be used instead.

* `stateval_tbl` now contains a `grp_var` column used to assign state values to "groups" of patients. This is distinct from `grp_id` in `hesim_data()`, which is used to define groups for subgroup analyses. 

* The public field `input_mats` has been renamed `input_data` in `R6` classes for disease progression and state values. This is a more generic name and will allow for potential feature enhancements in which `input_data` is a data frame rather than a matrix.

* `rdirichlet_mat()` has been modified to better facilitate sampling from transition matrices within the context of `define_model()` and `tparams_transprobs()`. One implication is that the number of rows in `alpha` must now be less than or equal to the number of columns and that the number of columns can be greater than the number of rows. 

## hesim 0.2.3
Remove a documented `...` that was not used in `weibullNMA()`.

## hesim 0.2.2
No longer use deprecated `C++` function `bind2nd()`.

## hesim 0.2.1
The `input_mats` class now contains an element `time_reset`. If `TRUE`, then time intervals reset each time a patient enters a new health state. In other words, state values can depend on time since entering a health state.

To illustrate, consider an oncology application with three health states (stable disease, progressed disease, and death). In these models it is common to assume that patients begin second line treatment after disease progression. Suppose the second line treatment is a chemotherapy that patients take for 12 cycles (or approximately 1 year). Then drug costs would accrue for the first year but not afterwards.

State values like this can be specified by setting `time_reset = TRUE` in `create_StateVals.stateval_tbl()`.

```r
hesim_dat <- hesim_data(strategies = data.frame(strategy_id = c(1, 2)),
                        patients = data.frame(patient_id = seq(1, 3)),
                        states = data.frame(state_id = c(1, 2)))
drugcosts <- stateval_tbl(tbl = data.frame(state_id = rep(c(1, 2), each = 2),
                                           time_start = c(0, 1, 0, 1),
                                           est = c(10000, 0, 12500, 0)),
                                  dist = "fixed",
                                  hesim_data = hesim_dat)  
drugcostsmod <- create_StateVals(drugcosts, time_reset = TRUE) 
```

## hesim 0.2.0
### Highlights
`hesim` now provides a general framework for integrating statistical models with economic evaluation. Users 
build a decision model by specifying a model structure, which consists of a set of statistical models for disease progression, utilities, and costs. Each statistical model is used to simulate outcomes as a function of estimated parameters and input data. N-state partitioned survival models (PSMs) and individual-level continuous time state transition models (iCTSTMs) are now supported. 

### New features
Economic models---which combine the disease, utility, and cost models---are constructed with the following classes:
* `hesim::Psm()` for PSMs
* `hesim::IndivCtstm()` for iCTSTMs

Disease models are constructed using the classes:
* `hesim::PsmCurves` to simulate survival curves for each endpoint of interest
* `hesim::IndivCtstmTrans` to simulate health state transitions with a iCTSTM

Utility and cost models are constructed with the `hesim::StateVals` class. 

The economic models are used to simulate disease progression (`$sim_disease()`, `$sim_stateprobs()`), quality-adjusted life-years (QALYs) (`$sim_qalys()`), and costs (`$sim_costs()`). Parameter uncertainty is propagated to model outcomes using probabilistic sensitivity analysis. Summaries of the simulated costs and QALYs are used to perform model-based cost-effectiveness analyses (CEAs) and represent decision uncertainty with `icea.ce()` and `icea_pw.ce()`.

### API changes
* The argument `sim` was renamed `sample` in `icea()`, `icea_pw()`, and `incr_effect()`.
* Custom functions and variables are no longer supported in `icea()` and `icea_pw()`.

## hesim 0.1.0
The initial CRAN submission containing support for CEA but not for model development. Decision uncertainty is represented using cost-effectiveness planes, cost-effectiveness acceptability curves, cost-effectiveness acceptability frontiers, and the expected value of perfect information. CEAs by subgroup (i.e., individualized CEAs) are performed with `icea()` and `icea_pw()`.