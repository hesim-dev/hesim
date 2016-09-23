library("foreign")
library("nnet")

# matrixC ----------------------------------------------------------------------
matrixC(c(1, 2, 3, 4, 5, 6), 2, 3)

# mlogit_prob ------------------------------------------------------------------
dat <- read.dta("http://www.ats.ucla.edu/stat/data/hsbdemo.dta")
dat$prog2 <- relevel(dat$prog, ref = "academic")
mlogit <- multinom(prog2 ~ ses + write, data = dat)
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

# markov_mlogit ----------------------------------------------------------------
beta.array <- array(NA, dim = c(2, length(beta), 3))
for (i in 1:dim(beta.array)[3]){
  for (j in 1:dim(beta.array)[1]){
    beta.array[j,, i] <- beta
  }
}
markov_mlogitC(t(as.matrix(x)), beta.array, c(100, 0, 0), 20, 100)

# markov_pmat ------------------------------------------------------------------
p <- array(NA, dim = c(1, 4, 3))
p[,, 1] <- c(.2, .8, .1, .9)
p[,, 2] <- c(.3, .7, .5, .5)
p[,, 3] <- c(.4, .6, .6, .4)
p.index <- c(0, 1, 1, 2)
z0 <- matrix(rep(c(1000, 0), 4), nrow = 4, ncol = 2, byrow = T)
markov_pmatC(p, p.index, z0, 19, 100, 2)
