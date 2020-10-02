## Release summary
This patch fixes a bug in the compiled code identified with the CRAN package check results and noted in "Additional isssues". It also makes a few minor updates to the documentation.

## Test environments
* Local OS X, R 4.0.0
* Ubuntu 16.04 (on travis-ci), R 4.0.2
* Ubuntu 18.04 (with valgrind), R 4.0.0
* R-hub builder
* win-builder (devel, release)

## R CMD check results
0 errors | 0 warnings | 1 note

* checking installed package size ... NOTE
    installed size is  5.4Mb
    sub-directories of 1Mb or more:
      doc    2.0Mb
      libs   1.9Mb
