#ifndef RCPP_hesim_H_
#define RCPP_hesim_H_

#include "hesim_RcppExports.h"
#include <hesim/statmods.h>
#include <hesim/time_fun.h>

/** \mainpage 

# Overall Description
This site documents the @c hesim C++ library, which is desisgned for use with 
the @c R package by the same name.
 */

/** @defgroup stats Statistics
 *  Probability distributions, random number generation, and statistical functions.
 */

/** @defgroup statmods Statistical models
 *  Prediction and random sampling from different statistical models.
 */

/** @defgroup statevals Health state values
 *  Simulate costs and QALYs based on time spent in health states.
 */

/** @defgroup test Tests
 *  Tests of C++ code. These functions are exported to @c R using the
 *  @c R package @c Rcpp and the @c R package @c testthat.
 */

#endif // RCPP_hesim_H_GEN_
