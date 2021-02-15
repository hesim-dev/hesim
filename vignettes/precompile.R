setwd("vignettes")
knitr::knit("benchmarks.Rmd.orig", output = "benchmarks.Rmd")
setwd("..")
pkgdown::build_article("benchmarks")