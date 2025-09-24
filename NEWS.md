## hesim 0.5.6
Ensure compatibility with `ggplot v4.0.0`. 

## hesim 0.5.5
Minor updates to the `.Rd` files and tests to fix problem identified with the CRAN package checks.

## hesim 0.5.4
* For individual CTSTMs, adds experimental support for `clock="mixt"` option with a `transition_types` argument for transition-specific clocks (#106).
* Fixes several CRAN checks that were producing NOTEs related to package
documentation (#114), Rd files (#118), and C++ system requirements (#119).

## hesim 0.5.3
Minor updates to the `.Rd` files to fix problems with the HTML version of the manual identified with the CRAN package checks.

## hesim 0.5.2

### API change
* The default value of the `check` argument of `eval_rng()` has been changed to `TRUE`.

### Bug fixes
* The Markov inhomogeneous individual-level model vignette has been corrected so that `eval_rng()` can be summarized using `summary.eval_rng()`.  

* `print.params_surv()` has been fixed for the piecewise exponential distribution.

* `custom()` works properly when `n=1`.

* Modifies `expmat()` to eliminate `length > 1 in coercion to logical` warnings in R-devel caused by `msm::MatrixExp()`.

* Fixed a bug introduced in the previous version in which the `clock` argument of `create_IndivCtstmTrans()` was not being passed to `IndivCtstmTrans$new()`. 

## hesim 0.5.1

### New features

* A `survival` object can now be constructed manually with `survival()` and simulated state probabilities can be computed from survival curves with `sim_stateprobs.survival()`. These features are useful for partitioned survival analyses when a user would like to use a survival model not supported by `hesim` to predict survival curves (#49).

* `tpmatrix()` is more flexible: 
    + There is now a `complement` argument where users can define transitions that are "complements". This is particularly helpful for creating large transition probability matrices since there is no longer a need to manually enter `"C"` in the `...` portion of `tpmatrix()`. In other words, users can pass a single object to `tpmatrix()` and use the `complement` argument to ensure probabilities sum to 1 in each row. One reasonable workflow with `define_tparams()` would be to (i) create a single matrix of initial values (say zeros), (ii) modify the transitions that differ from the initial values, and (iii) pass the resulting object to `tpmatrix()` while using the `complement` argument.
    
    + The `states`, `prefix`, and `sep` arguments can be used to name the columns (i.e., the transitions) and the states. The named states are, in turn, passed to the `$value` element of a `tparams_transprobs` object with `tparams_transprobs.tpmatrix()`. 

* An `eval_rng()` object can now be passed directly to `define_model()`, meaning that users can sample parameters prior to defining a model. Previously users could only pass an `rng_def` object to `define_model()`, which meant that sampling could only occur while evaluating a `model_def` object. 

* There are now summary and print methods for parameter, transformed parameter, and `eval_rng` objects. See `summary.params()`, `summary.tparams_transprobs()`, `summary.tparams_mean()`, and `summary.eval_rng()`.

* An `input_mats` object can be converted to a `data.table` with `as.data.table.input_mats()` and printed to the console in a less verbose way than in prior versions with `print.input_mats()`.

* `sim_ev()`, `sim_costs()`, and `sim_qalys()` are now exported functions that give users additional flexibility in their modeling pipelines and provide improved documentation for computation of expected values in cohort models. `sim_ev()` is particularly useful for computing outcomes that depend on time in state other than costs or quality-adjusted life-years (QALYs). 

* Multiple absorbing states (or none at all) are possible in `hesim::CohortDtstm` and `hesim::IndivCtstm` models (previously the final health state was always assumed to be a death state). In cohort models, the absorbing states can be set manually using the `absorbing` field in the `hesim::CohortDtstmTrans` class; if not, they are set automatically based on the transition probabilities. The number of health states in state value models (class `hesim::StateVals`) must equal the number of health states in the transition models less the number of absorbing states. 

* A new `create_CohortDtstmTrans.params_mlogit_list()` method allows the transition component of a cohort discrete time state transition model (cDTSTM) to be created directly from multinomial logistic regression parameter objects.

* The coefficient elements of parameter objects can be constructed from any object (e.g., data frame) than can be passed to `as.matrix()` (rather than only from matrices as in previous versions). See, for instance, `params_surv()`. 

### Bug fixes

* `sim_stateprobs.survival()` handles scenarios where survival curves cross better, ensuring that state probabilities sum to 1 (#56).

## hesim 0.5.0
### New features
* A transition intensity matrix can be created from a multi-state model fit using `msm::msm()` with the  `qmatrix.msm()` method. Similarly, a `CohortDtstmTrans` object can be created with `create_CohortDtstmTrans.msm()`.

* The `...` argument to `create_PsmCurves()` now passes arguments to `create_params.partsurvfit()` when `object` is of class `flexsurvreg_list`. This allows more control over bootstrapping (i.e., use of the `max_errors` and `silent` arguments). 

* `summary.ce()` is a new summary method for a `hesim::ce` object that produces confidence intervals for QALYs and each cost category; `format.summary.ce()` formats the output for pretty printing. 

* `icer()` generates a tidy table of incremental cost-effectiveness ratios (ICERs) given output from `cea_pw()`; `format.icer()` formats the output for pretty printing.

* `plot_ceplane()`, `plot_ceac()`, `plot_ceaf()`, and `plot_evpi()` plot the cost-effectiveness plane, cost-effectiveness acceptability curve (CEAC), cost-effectiveness acceptability frontier (CEAF), and expected value of perfect information (EVPI), respectively. 

* `autoplot.survival()` and `autoplot.stateprobs()` plot survival curves and state probabilities, respectively.

### API changes
* The first column of each matrix listed in the `coef` element returned by `create_params.flexsurvreg()` is now named "(Intercept)" instead of the name of the corresponding parameter.

* The `create_params()` methods now use the argument `uncertainty` to draw parameters  and the old arguments `point_estimate` and `bootstrap` are deprecated. This also affects `create_CohortDtstmTrans()`, `create_IndivCtstmTrans()`, and `create_PsmCurves()`.

* `icer_tbl()` has been deprecated in favor of `icer()`.

* The column `trans` in the data table returned by the `$hazard()` and `$cumhazard()` methods from the `hesim::CtstmTrans` class has been renamed `transition_id`. 

### Documentation

* Performance benchmarks are now provided [here](../articles/benchmarks.html).

### Bug fixes
* The `$summarize()` method of `hesim::Psm` now contains the `by_grp` argument so that subgroup analyses can be performed.

## hesim 0.4.2
### New features
* There are new functions to construct (and debug the construction of) the multiple transition probability matrices stored in `tparams_transprobs()` objects and used for cDTSTMs. These can either stored as 3-dimensional arrays or as 2-dimensional tabular objects (i.e., `data.table`, `data.frame`, `matrix`).

    + `as_array3()` and `as_tbl2()` lets users convert 2-dimensional tabular objects where each row stores a flattened square matrix to a 3-dimensional array of square matrices and vice versa.
    
    + `qmatrix()` lets you store transition intensity matrices which can be used to construct transition probability matrices with the matrix exponential via `expmat()`. The latter is a simple wrapper around `msm::MatrixExp()` that computes the matrix exponential for all matrices in an array rather than just a single matrix.
    
    + `apply_rr()` applies relative risks (stored in a 2-dimensional tabular object) to (potentially multiple elements) of transition probability matrices stored in an array. This function is vectorized so it can be performed very quickly even for large arrays.
    
    + `as.data.table.tparams_transprobs()` converts the array of transition probability matrices stored in a `tparams_transprobs` object into a `data.table` which can be helpful for debugging to ensure that the right transition probability matrices correspond to the right observations (i.e., treatment strategies, patients, etc.).
    
    + The `tpmatrix` element of `define_tparams()` can now be a 3-dimensional array in addition to the output of `tpmatrix()` to increase flexibility for the user.
    
* 2-dimensional tabular objects (in addition to vectors) can now be passed to `...` with `tpmatrix()`. See the new examples. 
    
* A new dataset `hesim::onc3` was added as an example multi-state dataset for an oncology use case with 3 health states (Stable, Progression, Death) and 3 possible transitions (Stable -> Progression, Stable -> Death, and Progression -> Death). This is similar to `hesim::mstate3_exdata` but does not allow for reversible transitions and does not contain cost or utility data.

* The function `as_pfs_os()` can convert a multi-state dataset in the same format as `hesim::onc3` into a dataset with one row per patients containing time to event information for progression free survival (PFS) and overall survival (OS).

### Bug fixes
* The `cycle_length` field in `CohortDtstmTrans` was fixed so that it corresponds to a model cycle in terms of years (e.g., a value of 2 means a model cycle is 2 years long and that state probabilities are computed every 2 years with `$sim_stateprobs()`).

* The simulated dataset `multinom3_exdata` was fixed by removing a bug where some patients were simulated to have died more than once. 

* Fixed bug where the `$sim_costs()` method of `IndivCtstm` was erroneously returning a life-years column in addition to the costs column.

* Modification to creation of input matrices from a `flexsurvreg` object to properly capture levels of factor variables.  

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
`hesim` now supports cDTSTM via `hesim::CohortDtstm` objects. Users can build a model by either fitting multinomial logistic regressions with `nnet::multinom()` or with a mathematical expression using `define_model()`. Furthermore, `$summarize()` methods now have a `by_grp` option to facilitate subgroup analyses. 

### New features

* The `hesim::CohortDtstm` class simulates cDTSTMs. State transitions in a cDTSTM are simulated using the `hesim::CohortDtstmTrans` class, which can be constructed from a `multinom_list()` object or using `define_model()`.

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
