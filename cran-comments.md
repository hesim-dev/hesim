## Resubmission
This is a resubmission. 

In the original submission, it was noted that there was an error in one of the package tests on a "noLD" platform. This was due to a rounding error and corrected by adding tolerance to the problematic test. The updated package has been checked using R-Hub builder on a platform with long doubles disabled and returned no errors.

In addition, there was an incorrect local URL in the `psm.Rmd` vignette. This has been corrected. 

## Release summary
This release contains a minor update as outlined on the package 
website [here](https://hesim-dev.github.io/hesim/news/index.html).

## Test environments
* Local OS X, R 3.6.0
* Ubuntu 14.04 (on travis-ci), R 3.6.0
* win-builder (devel and release)

## R CMD check results
0 ERRORs | 0 WARNINGs | 1 NOTE

* checking installed package size ...
    installed size is  7.2Mb
    sub-directories of 1Mb or more:
      doc    4.6Mb
      libs   1.7Mb
      
  This is compiled code in the libs/ directory and documentation in the docs/ directory. 
