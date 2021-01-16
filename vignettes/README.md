Vignettes for the package are stored in this directory. Note that the `benchmarks` vignette is not currently included in the CRAN release and must be precompiled with:

```{r}
knitr::knit("vignettes/benchmarks.Rmd.orig", output = "vignettes/benchmarks.Rmd")
```