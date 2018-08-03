#ifndef hesim_H_
#define hesim_H_

#include <hesim/statmods.h>
#include <hesim/time_fun.h>

/** \mainpage 

# Overall Description
This site documents the @c hesim C++ library, which is desisgned for use with 
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

/** @defgroup test Tests
 *  Tests of C++ code. These functions are exported to @c R using the
 *  @c R package @c Rcpp and the @c R package @c testthat.
 */

#endif // hesim_H_
