## Release summary
This release contains a minor update to address an email from Kurt Hornik to fix a warning about  'Documented arguments not in \\usage' in the r-devel checks from recent bug fix (PR#16223).

## Test environments
* Local OS X, R 3.6.3
* Ubuntu 16.04 (on travis-ci), R 3.6.2
* win-builder (devel and release)

## R CMD check results
0 ERRORs | 0 WARNINGs | 1 NOTE

* checking installed package size ...
    installed size is  7.2Mb
    sub-directories of 1Mb or more:
      doc    4.6Mb
      libs   1.7Mb
      
  This is compiled code in the libs/ directory and documentation in the docs/ directory. 
