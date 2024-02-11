## Release summary
The primary purpose of this release is to fix problems identified by the 
CRAN checks. Additional details of the release can be found at:
https://hesim-dev.github.io/hesim/news/index.html#hesim-054

This is a re-submission to fix additional NOTEs caught during the submission:

```
Flavor: r-devel-linux-x86_64-debian-gcc, r-devel-windows-x86_64
Check: S3 generic/method consistency, Result: NOTE
  Mismatches for apparent methods not registered:
  create_StateVals:
    function(object, ...)
  create_StateVals.eval_model:
    function(object, cost, name, init_args)

  check:
    function(object, ...)
  check.params_lm:
    function(object)

  check:
    function(object, ...)
  check.params_surv:
    function(object)

  check:
    function(object, ...)
  check.eval_model:
    function(x)

  check:
    function(object, ...)
  check.hesim_data:
    function(x)

  check:
    function(object, ...)
  check.input_mats:
    function(object)

  check:
    function(object, ...)
  check.tparams_transprobs:
    function(object)

  check:
    function(object, ...)
  check.id_attributes:
    function(object)

  check:
    function(object, ...)
  check.params_mlogit:
    function(object)

  check:
    function(object, ...)
  check.eval_rng:
    function(object)

  check:
    function(object, ...)
  check.coeflist:
    function(coefs)

  check:
    function(object, ...)
  check.params_surv_list:
    function(object)

  check:
    function(object, ...)
  check.model_def:
    function(x)

  new_tparams_transprobs:
    function(object, ...)
  new_tparams_transprobs.data.table:
    function(object)

  new_tparams_transprobs:
    function(object, ...)
  new_tparams_transprobs.array:
    function(object, tpmatrix_id, times, grp_id, patient_wt)

  Apparent methods for exported generics not registered:
    create_StateVals.eval_model plot_ceac.default sim_ev.NULL
    tparams_transprobs.eval_model
```


## Test environments
* Local OS X, R 4.2.2
* Ubuntu 20.04.6 (on GitHub actions), R-devel, R 4.3.2
* Microsoft Windows Server 2022 10.0.20348 (on GitHub actions) R 4.3.2
* win-builder (devel, release)
* R-hub builder

## Win Builder R-devel results
0 errors | 0 warnings | 0 notes
