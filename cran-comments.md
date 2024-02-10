## Release summary
The primary purpose of this release is to fix problems identified by the 
CRAN checks. Additional details of the release can be found at:
https://hesim-dev.github.io/hesim/news/index.html#hesim-054

## Test environments
* Local OS X, R 4.2.2
* Ubuntu 20.04.6 (on GitHub actions), R-devel, R 4.3.2
* Microsoft Windows Server 2022 10.0.20348 (on GitHub actions) R 4.3.2
* win-builder (devel, release)
* R-hub builder

## Local R CMD check results
0 errors | 0 warnings | 2 notes

* checking installed package size ... NOTE
    installed size is  5.9Mb
    sub-directories of 1Mb or more:
      doc    2.2Mb
      libs   1.8Mb

* checking dependencies in R code ... NOTE
  Namespace in Imports field not imported from: ‘R6’
    All declared Imports should be used.
