## Release summary
The purpose of this release is to fix memory access issues identified by 
`clang-ASAN` and `gcc-ASAN` during the CRAN checks. 

The error was first reproduced on R-hub `clang-asan` and `gcc-asan`. After
the fix, it was then confirmed that `clang-asan` and `gcc-asan` did not 
produce any errors.

## R CMD check results
0 errors | 0 warnings | 0 notes
