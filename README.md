
# Health economic simulation modeling <img src="man/figures/logo.png" align="right" width="90" />

<!-- badges: start -->

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/hesim)](https://cran.r-project.org/package=hesim)
[![R build
status](https://github.com/hesim-dev/hesim/workflows/R-CMD-check/badge.svg)](https://github.com/hesim-dev/hesim/actions)
[![Coverage
Status](https://codecov.io/gh/hesim-dev/hesim/branch/master/graph/badge.svg)](https://codecov.io/gh/hesim-dev/hesim)
<!-- badges: end -->

## Overview

`hesim` is a modular and computationally efficient R package for health
economic simulation modeling and decision analysis that provides a
general framework for integrating statistical analyses with economic
evaluation. The package supports cohort discrete time state transition
models (DTSTMs), N-state partitioned survival models (PSMs), and
individual-level continuous time state transition models (CTSTMs),
encompassing both Markov (time-homogeneous and time-inhomogeneous) and
semi-Markov processes. It heavily utilizes `Rcpp` and `data.table`,
making individual-level simulation, probabilistic sensitivity analysis
(PSA), and incorporation of patient heterogeneity fast.

Features of the current version can be summarized as follows:

  - Cohort DTSTMs, individual-level CTSTMs, and N-state PSMs that
    encompass Markov and semi-Markov processes
  - Options to build models directly from fitted statistical models or
    by defining them in terms of expressions
  - Parameter estimates from either an `R` based model or from an
    external source
  - Convenience functions for sampling model parameters from parametric
    distributions or via bootstrapping
  - Parameter uncertainty propagated with PSA
  - Modeling patient heterogeneity
  - Performing cost-effectiveness analyses and representing decision
    uncertainty from PSAs
  - Simulation code written in `C++` to boost performance

## Installation

You can install the [current
release](https://hesim-dev.github.io/hesim/) from CRAN or the most up to
date [development version](https://hesim-dev.github.io/hesim/dev/) from
GitHub.

``` r
# Install from CRAN:
install.packages("hesim")

# Install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("hesim-dev/hesim")
```

## Getting started

There are two good places to start:

1.  The [Introduction to
    `hesim`](https://hesim-dev.github.io/hesim/articles/intro.html)
    article provides a quick introduction.

2.  Our [preprint](https://arxiv.org/abs/2102.09437) describes the
    package (including mathematical details) more thoroughly.

You might also want to explore our example analyses which can be found
in the preprint and web articles. They are summarized in the table
below, with some drawn from the [Decision Modeling for Health Economic
Evaluation](https://www.herc.ox.ac.uk/downloads/decision-modelling-for-health-economic-evaluation)
textbook. Key areas of focus are the (i) statistical models of disease
progression (in terms of the baseline risk and relative treatment
effects) and (ii) the available data (either individual patient data
(IPD) or aggregate-level data).

<table class="table" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="empty-cells: hide;border-bottom:hidden;" colspan="1">

</th>

<th style="empty-cells: hide;border-bottom:hidden;" colspan="1">

</th>

<th style="empty-cells: hide;border-bottom:hidden;" colspan="1">

</th>

<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Baseline risk

</div>

</th>

<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Treatment effect

</div>

</th>

<th style="empty-cells: hide;border-bottom:hidden;" colspan="1">

</th>

</tr>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:left;">

Name

</th>

<th style="text-align:left;">

Model

</th>

<th style="text-align:left;">

Disease model

</th>

<th style="text-align:left;">

Disease data

</th>

<th style="text-align:left;">

Disease model

</th>

<th style="text-align:left;">

Disease data

</th>

<th style="text-align:left;">

Application

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

1

</td>

<td style="text-align:left;width: 15em; ">

<a href="https://arxiv.org/pdf/2102.09437.pdf" style="     ">Preprint
4.1</a>

</td>

<td style="text-align:left;">

iCTSTM

</td>

<td style="text-align:left;">

Multi-state model

</td>

<td style="text-align:left;">

IPD

</td>

<td style="text-align:left;">

Coefficient (AFT)

</td>

<td style="text-align:left;">

IPD

</td>

<td style="text-align:left;">

Oncology

</td>

</tr>

<tr>

<td style="text-align:left;">

2

</td>

<td style="text-align:left;width: 15em; ">

<a href="https://arxiv.org/pdf/2102.09437.pdf" style="     ">Preprint
4.2</a>

</td>

<td style="text-align:left;">

PSM

</td>

<td style="text-align:left;">

Survival models

</td>

<td style="text-align:left;">

IPD

</td>

<td style="text-align:left;">

Coefficient (AFT)

</td>

<td style="text-align:left;">

Aggregate

</td>

<td style="text-align:left;">

Oncology

</td>

</tr>

<tr>

<td style="text-align:left;">

3

</td>

<td style="text-align:left;width: 15em; ">

<a href="https://arxiv.org/pdf/2102.09437.pdf" style="     ">Preprint
4.3</a>

</td>

<td style="text-align:left;">

cDTSTM

</td>

<td style="text-align:left;">

Multi-state model (panel data)

</td>

<td style="text-align:left;">

IPD

</td>

<td style="text-align:left;">

RR

</td>

<td style="text-align:left;">

Aggregate

</td>

<td style="text-align:left;">

Oncology

</td>

</tr>

<tr>

<td style="text-align:left;">

4

</td>

<td style="text-align:left;width: 15em; ">

<a href="https://hesim-dev.github.io/hesim/articles/markov-cohort.html" style="     ">Simple
Markov cohort</a>

</td>

<td style="text-align:left;">

cDTSTM

</td>

<td style="text-align:left;">

Multinomial

</td>

<td style="text-align:left;">

Aggregate

</td>

<td style="text-align:left;">

RR

</td>

<td style="text-align:left;">

Aggregate

</td>

<td style="text-align:left;">

HIV

</td>

</tr>

<tr>

<td style="text-align:left;">

5

</td>

<td style="text-align:left;width: 15em; ">

<a href="https://hesim-dev.github.io/hesim/articles/markov-inhomogeneous-cohort.html" style="     ">Time
inhomogeneous Markov (cohort)</a>

</td>

<td style="text-align:left;">

cDTSTM

</td>

<td style="text-align:left;">

Custom

</td>

<td style="text-align:left;">

Aggregate

</td>

<td style="text-align:left;">

Coefficient (HR)

</td>

<td style="text-align:left;">

Aggregate

</td>

<td style="text-align:left;">

Hip replacement

</td>

</tr>

<tr>

<td style="text-align:left;">

6

</td>

<td style="text-align:left;width: 15em; ">

<a href="https://hesim-dev.github.io/hesim/articles/mlogit.html" style="     ">Multinomial
logit</a>

</td>

<td style="text-align:left;">

cDTSTM

</td>

<td style="text-align:left;">

Multinomial logit

</td>

<td style="text-align:left;">

IPD

</td>

<td style="text-align:left;">

Coefficient (OR)

</td>

<td style="text-align:left;">

IPD

</td>

<td style="text-align:left;">

Generic

</td>

</tr>

<tr>

<td style="text-align:left;">

7

</td>

<td style="text-align:left;width: 15em; ">

<a href="https://hesim-dev.github.io/hesim/articles/markov-inhomogeneous-indiv.html" style="     ">Time
inhomogeneous Markov (individual)</a>

</td>

<td style="text-align:left;">

iCTSTM

</td>

<td style="text-align:left;">

Custom

</td>

<td style="text-align:left;">

Aggregate

</td>

<td style="text-align:left;">

Coefficient (HR)

</td>

<td style="text-align:left;">

Aggregate

</td>

<td style="text-align:left;">

Hip replacement

</td>

</tr>

<tr>

<td style="text-align:left;">

8

</td>

<td style="text-align:left;width: 15em; ">

<a href="https://hesim-dev.github.io/hesim/articles/mstate.html" style="     ">Semi-Markov
multi-state</a>

</td>

<td style="text-align:left;">

iCTSTM

</td>

<td style="text-align:left;">

Multi-state model

</td>

<td style="text-align:left;">

IPD

</td>

<td style="text-align:left;">

Coefficient (AFT)

</td>

<td style="text-align:left;">

IPD

</td>

<td style="text-align:left;">

Generic

</td>

</tr>

<tr>

<td style="text-align:left;">

9

</td>

<td style="text-align:left;width: 15em; ">

<a href="https://hesim-dev.github.io/hesim/articles/psm.html" style="     ">4-state
PSM</a>

</td>

<td style="text-align:left;">

PSM

</td>

<td style="text-align:left;">

Survival models

</td>

<td style="text-align:left;">

IPD

</td>

<td style="text-align:left;">

Coefficient (AFT)

</td>

<td style="text-align:left;">

IPD

</td>

<td style="text-align:left;">

Oncology

</td>

</tr>

</tbody>

<tfoot>

<tr>

<td style="padding: 0; border: 0;" colspan="100%">

<span style="font-style: italic;">Note: </span> <sup></sup> iCTSTM =
Individual-level continuous time state transition model; PSM =
partitioned survival model; cDTSTM = Cohort discrete time state
transition model. AFT = accelerated failure time; RR = relative risk; HR
= hazard ratio; OR = odds ratio. IPD = individual patient data.

</td>

</tr>

</tfoot>

</table>

## Citing hesim

If you use `hesim`, please cite as follows:

``` 

  Devin Incerti and Jeroen P Jansen (2021). hesim: Health Economic
  Simulation Modeling and Decision Analysis. arXiv:2102.09437
  [stat.AP], URL https://arxiv.org/abs/2102.09437.

A BibTeX entry for LaTeX users is

  @Misc{incerti2021hesim,
    author = {Devin Incerti and Jeroen P. Jansen},
    title = {hesim: Health Economic Simulation Modeling and Decision Analysis},
    year = {2021},
    eprint = {2102.09437},
    archiveprefix = {arXiv},
    primaryclass = {stat.AP},
    url = {https://arxiv.org/abs/2102.09437},
  }
```
