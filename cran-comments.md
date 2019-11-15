## Release summary
This release contains a minor update to address an email from Professor Brian Ripley that the `auto_ptr`, `bind1st`, `bind2nd`, `mem_fun_ref`, and `ptr_fun` are now deprecated in`C++` and should be removed from `R` packages. `hesim` contained a function that used `bind2nd`, but it is now replaced with an alternative implementation.

## Test environments
* Local OS X, R 3.6.1
* Ubuntu 16.04 (on travis-ci), R 3.6.1
* win-builder (devel and release)

## R CMD check results
0 ERRORs | 0 WARNINGs | 1 NOTE

* checking installed package size ...
    installed size is  7.2Mb
    sub-directories of 1Mb or more:
      doc    4.6Mb
      libs   1.7Mb
      
  This is compiled code in the libs/ directory and documentation in the docs/ directory. 
