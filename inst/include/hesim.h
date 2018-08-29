#ifndef hesim_H_
#define hesim_H_

#include <hesim/statmods/statmods.h>
#include <hesim/time_fun.h>

/** \mainpage 

# Overall Description
This site documents the @c hesim C++ library, which is designed for use with 
the @c R package by the same name.
 */

/** @defgroup general General
 *  General classes and functions to facilitate model development. 
 */

/** @defgroup math Mathematics
 *  Functions and classes for numerical computing including integration and
 *  root finding.
 */

/** @defgroup stats Statistics
 *  Probability distributions, random number generation, and statistical functions.
 */

/** @defgroup statmods Statistical models
 *  Prediction and random sampling from different statistical models.
 */

/** @defgroup psm Partitioned survival models
 *  Classes and functions for simulating partitioned survival models.
 */

/** @defgroup ctstm Continuous time state transiton models
 *  Classes and functions for simulating continuous time state transiton models.
 */

/** @defgroup Rbase R base
 *  Functions from the @c R standard library. Some functions are adapted for use
 *  with @c C++.
 */

/** @defgroup test Tests
 *  Tests of C++ code. These functions are exported to @c R using the
 *  @c R package @c Rcpp and the @c R package @c testthat.
 */

#endif // hesim_H_
