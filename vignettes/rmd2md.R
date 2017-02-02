library("knitr")
opts_chunk$set(fig.path = "figs/", fig.width = 8, fig.height = 5)
opts_knit$set(base.url = "", root.dir = getwd())
knit("cea-vignette.Rmd", "cea-vignette.md")