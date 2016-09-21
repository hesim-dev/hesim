library("foreign")
library("nnet")

# matrixC ----------------------------------------------------------------------
matrixC(c(1, 2, 3, 4, 5, 6), 2, 3)

# mlogit_prob ------------------------------------------------------------------
dat <- read.dta("http://www.ats.ucla.edu/stat/data/hsbdemo.dta")
dat$prog2 <- relevel(dat$prog, ref = "academic")
mlogit <- multinom(prog2 ~ ses + write, data = ml)
beta <- c(coef(mlogit))
predict(mlogit, newdata = dat[1, c("ses", "write")], "probs")
x <- model.matrix(mlogit)[1, ]
mlogit_prob(x, beta, 3)

# mlogit_transprob -------------------------------------------------------------
beta3 <- matrix(rep(beta, 3), nrow = 3, byrow = T)
pmat <- mlogit_transprob(x, beta3, 3)

# transprob_addmort ------------------------------------------------------------
pmat <- transprob_addmort(pmat, c(.2, .3, .4), 3)
rowSums(pmat)
