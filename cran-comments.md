## Resubmission
This is a resubmission. The original submission generated the error:

*checking re-building of vignette outputs ... [3s] WARNING
Error in re-building vignettes:
  ...
Error: processing vignette 'cea-vignette.Rmd' failed with diagnostics:
DLL 'stringi' not found: maybe not installed for this architecture?
Execution halted*

It was suggested that this was a likely a race condition and I was told to resubmit. 

In addition, I have made the following changes:

* Converted the DESCRIPTION title to title case.

* Removed 'An R package' from the description and elaborated on the functionality 
provided by the package.

* Added references for the methods described in the Description' field

* Added more small executable examples in my RD files.

## Test environments
* OS X, R 3.4.1
* CentOS release 6.7 R 3.2.3
* Ubuntu 14.04 (on travis-ci), R 3.4.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs, WARNING or NOTEs on OS X, CENTOS or Ubuntu

There were no ERRORs or WARNING on win-builider. There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Devin Incerti <devin.incerti@gmail.com>'

This note is occuring because this is a new submission.
