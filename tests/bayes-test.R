rm(list = ls())
library("foreign")
library("MCMCpack")
library("nnet")

# data
# dv = apply: likelihood of applying to graduate school
dat <- read.dta("http://www.ats.ucla.edu/stat/data/ologit.dta")
newdat <- data.frame(
  pared = rep(0:1, 200),
  public = rep(0:1, each = 200),
  gpa = rep(seq(from = 1.9, to = 4, length.out = 100), 4))

# predict_MCMCmnl -------------------------------------------------------------
mlogit <- multinom(apply ~ pared + public + gpa, data = dat)
bmlogit <- MCMCmnl(apply ~ pared + public + gpa, data = dat, mcmc = 5000,
                   mcmc.method = c("IndMH"), baseline = "unlikely")
x <- cbind(1, as.matrix(newdat))
mlogit_prob(coef(mlogit), x[1, ], 3)
predict(mlogit, newdata = newdat[1, ], "prob")
mlogit_prob(matrix(as.matrix(bmlogit)[1, ], nrow = 2, ncol = 4), x[1, ], 3)
head(predict_MCMCmnl(as.matrix(bmlogit), x, 3)[, ,1])

