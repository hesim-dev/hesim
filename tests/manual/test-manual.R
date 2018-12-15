context("Manual ctstm.R unit tests")
library("flexsurv")
library("mstate")
library("data.table")
library("Rcpp")
library("ggplot2")
rm(list = ls())

# Simulate flexsurv fit and compare state probabilities ------------------------
n_pats <- 10000
tmat <- rbind(c(NA, 1, 2), 
              c(NA, NA, 3), 
              c(NA, NA, NA))
transitions <- create_trans_dt(tmat)
transitions[, trans := factor(transition_id)]
hesim_dat <- hesim_data(strategies <- data.table(strategy_id = 1),
                        patients <- data.table(patient_id = 1:n_pats),
                        transitions = transitions)
hesim_edat <- expand(hesim_dat, by = c("strategies", "patients", "transitions"))

# Exponential model
fit <- flexsurvreg(Surv(years, status) ~ trans, data = bosms3, dist = "weibull")

## Simulate flexsurv
flexsurv_stateprobs <- pmatrix.simfs(fit, M = n_pats, t = 5, trans = tmat)

## Simulate hesim
transmod <- create_IndivCtstmTrans(fit, 
                                   data = hesim_edat, trans_mat = tmat,
                                   point_estimate = TRUE)
hesim_stateprobs <- transmod$sim_stateprobs(t = 5)

## Compare
stateprobs <- cbind(flexsurv_stateprobs[1, ],
                    hesim_stateprobs$prob)
colnames(stateprobs) <- c("flexsurv", "hesim")
print(stateprobs)


# Simulate survival from arbitrary cumulative hazards --------------------------
module <- Rcpp::Module('distributions', PACKAGE = "hesim")

compute_surv <- function(step, lower, upper, hazfun){
  time <- seq(lower, upper, by = step)
  cumhaz <- rep(NA, length(time))
  cumhaz[1] <- 0
  for (i in 2:length(time)){
    cumhaz[i] <- (step * do.call("hazfun", list(time[i]))) + cumhaz[i - 1]
  }
  surv <- exp(-cumhaz)
  dat <- data.frame(time = time, surv = surv, lab = "Analytical") 
  return(dat)
}

# Test #1 = fractional polynomial from [0, inf)
FracPoly <- module$fracpoly
gamma = c(-1.2, -.567, 1.15)
powers = c(1, 0)
fp <- new(FracPoly, gamma = gamma, powers = powers)
fp$max_x_ <- 40

## Random sample with hesim
r1 <- replicate(10000, fp$random())
fun <- ecdf(r1)
esurv <- 1 - fun(time)
dat1 <- data.frame(time = time, surv = esurv, lab = "Random")

## Analytically compute survival
dat2 <- compute_surv(step = 1/12, lower = 0, upper = fp$max_x_, 
                     hazfun = fp$hazard)

## Compare
dat <- rbind(dat1, dat2)
ggplot(dat, aes(x = time, y = surv, col = lab)) + geom_line()

# Test #2 = truncated exponential distribution
Exponential <- module$exponential
exp <- new(Exponential, rate = 1.5)
lower <- 5; upper <- 10
step <- 1/12

## Sample using inverse CDF method
r1 <- replicate(1000, exp$trandom(lower, upper))

## Sample from arbitrary cumulative hazard with R
surv_df <- compute_surv(step = step, lower = lower, upper = upper,
                        hazfun = exp$hazard)
r2 <- replicate(1000,
                hesim:::C_test_rsurv(time = surv_df$time, est = surv_df$surv, 
                           type = "surv", time_inf = FALSE))

## Compare
time <- seq(lower, upper, step)
fun1 <- ecdf(r1)
fun2 <- ecdf(r2)
esurv1 <- 1 - fun1(time)
esurv2 <- 1 - fun2(time)
dat1 <- data.frame(time = time, surv = esurv1, lab = "Empirical CDF")
dat2 <- data.frame(time = time, surv = esurv2, lab = "Empirical hazard")
dat <- rbind(dat1, dat2)
ggplot(dat, aes(x = time, y = surv, col = lab)) + geom_line() 
