Vignettes for the package are stored in this directory. Note that the `benchmarks` vignette is not currently included in the CRAN release and must be precompiled. The following code with generate the `R` Markdown file and build a `pkgdown` webpage. 

```{r}
setwd("vignettes")
knitr::knit("benchmarks.Rmd.orig", output = "benchmarks.Rmd")
setwd("..")
pkgdown::build_article("benchmarks")
```