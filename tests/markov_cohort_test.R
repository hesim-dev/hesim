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
pmatadj <- transprob_addmort(pmat, c(.2, .3, .4), 3, c(0, 0, 0), 1)
rowSums(pmatadj)

# markov_mlogit ----------------------------------------------------------------
beta.array <- array(NA, dim = c(2, length(beta), 3))
for (i in 1:dim(beta.array)[3]){
  for (j in 1:dim(beta.array)[1]){
    beta.array[j,, i] <- beta
  }
}
markov_mlogitC(t(as.matrix(x)), beta.array, c(100, 0, 0), 20, 100)

# markov_pmat ------------------------------------------------------------------
p <- array(NA, dim = c(2, 2, 3))
p[,, 1] <- matrix(c(.2, .8, .1, .9), 2, 2, byrow = T)
p[,, 2] <- matrix(c(.3, .7, .5, .5), 2, 2, byrow = T)
p[,, 3] <- matrix(c(.4, .6, .6, .4), 2, 2, byrow = T)
p.index <- c(0, 1, 1, 2)
z0 <- matrix(rep(c(1000, 0, 0), 4), nrow = 4, ncol = 3, byrow = T)
mort.prob <- array(NA, dim = c(nrow(lifetable_male), 2, 1))
mort.prob[, , 1] <- matrix(lifetable_male$qx, nrow = dim(mort.prob)[1], ncol = 2)
markov_cohort_trans(z0 = z0, ncycles = rep(19, 4), pmat = p,
                    pmat_index = p.index, mortadj = TRUE,
                    mortprob = mort.prob, mortprob_index = rep(0, 4))
markov_cohort_trans(z0 = z0[, 1:2], ncycles = rep(19, 4), pmat = p,
                    pmat_index = p.index, mortadj = FALSE,
                    mortprob = mort.prob, mortprob_index = rep(0, 4))
