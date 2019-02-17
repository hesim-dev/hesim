context("Manual ctstm.R unit tests")
library("flexsurv")
library("mstate")
library("data.table")
library("Rcpp")
library("ggplot2")
# setwd("tests/manual") # tests should be run of this directory
source("../testthat/helpers.R")
rm(list = ls())

# State probabilities: hesim vs. mstate ----------------------------------------
# Model structure
tmat <- rbind(c(NA, 1, 2), 
              c(3, NA, 4), 
              c(NA, NA, NA))

# Fit survival models
fit_data <- data.table(ctstm3_exdata$transitions)
fit_data[, trans := factor(trans)]

fit_models <- function(clock = c("reset", "forward")){
  clock <- match.arg(clock)
  
  # Parametric models
  dists <- c("exponential", "weibull", "gompertz", "gamma", "lnorm", "llogis",
             "gengamma")
  fits <- vector(mode = "list", length = length(dists))
  names(fits) <- dists
  for (i in 1:length(fits)){
    if (clock == "reset"){
      fits[[i]] <-  flexsurvreg(Surv(years, status) ~ trans, data = fit_data, 
                                dist = dists[i]) 
    } else{
      fits[[i]] <-  flexsurvreg(Surv(Tstart, Tstop, status) ~ trans, data = fit_data, 
                                dist = dists[i]) 
    }
    print(paste0("Fit ", dists[i], " model"))
  } 
  
  # Spline model
  if (clock == "reset"){
    fits$spline <- flexsurvspline(Surv(years, status) ~ trans, data = fit_data) 
  } else{
    fits$spline <- flexsurvspline(Surv(Tstart, Tstop, status) ~ trans, data = fit_data) 
  }
  print(paste0("Fit spline model"))
  
  # Return
  return(fits)
}
fits_reset <- fit_models(clock = "reset")
fits_forward <- fit_models(clock = "forward")

# Compute cumulative hazards
hesim_msfit <- function(fit, tmat, t_grid){
  input_dat <- data.frame(strategy_id = 1, patient_id = 1, 
                          transition_id = 1, 
                          trans = factor(1:max(tmat, na.rm = TRUE)))
  setattr(input_dat, "class", c("expanded_hesim_data", "data.table", "data.frame"))
  attr(input_dat, "id_vars") <- c('strategy_id', "patient_id", "transition_id")
  transmod <- create_IndivCtstmTrans(fit, 
                                     input_data = input_dat, trans_mat = tmat,
                                     point_estimate = TRUE)
  cumhaz <- transmod$cumhazard(t = t_grid)[, .(t, cumhazard, trans)]
  setnames(cumhaz, c("t", "cumhazard", "trans"), c("time", "Haz", "trans"))
  return(cumhaz)
}

# Simulate 
## mstate
sim_stprobs_mstate <- function(fit, n_patients, t_grid,
                               clock = c("reset", "forward")){
  clock <- match.arg(clock)
  cumhaz <- hesim_msfit(fit, tmat, t_grid)
  if (clock == "reset"){
    stprobs <- mstate::mssample(Haz = cumhaz, trans = tmat, tvec = t_grid, 
                                clock = clock, M = n_patients) 
  } else{
    msfit <- list(Haz = cumhaz,
                  trans = tmat)
    class(msfit) <- "msfit"
    stprobs <- probtrans(msfit, predt = 0, variance = FALSE)[[1]]
  }
  stprobs <- data.table(stprobs)
  stprobs <- melt(stprobs, id.vars = "time", 
                  variable.name = "state", value.name = "prob")
  stprobs[, state := factor(state,
                            levels = paste0("pstate", 1:3),
                            labels = paste0("State ", 1:3))] 
  stprobs[, lab := "mstate"]
  return(stprobs)
}

## hesim
sim_stprobs_hesim <- function(fit, n_patients, t_grid, 
                              clock = c("reset", "forward")){
  clock <- match.arg(clock)
  # Input data
  transitions <- create_trans_dt(tmat)
  transitions[, trans := factor(transition_id)]
  hesim_dat <- hesim_data(strategies <- data.table(strategy_id = 1),
                          patients <- data.table(patient_id = 1:n_patients),
                          transitions = transitions)
  hesim_edat <- expand(hesim_dat, by = c("strategies", "patients", "transitions")) 
  
  # Simulate
  ## hesim
  transmod <- create_IndivCtstmTrans(fit, 
                                     input_data = hesim_edat, trans_mat = tmat,
                                     point_estimate = TRUE,
                                     clock = clock) 
  stprobs <- transmod$sim_stateprobs(t = t_grid)  
  stprobs[, state := factor(state_id,
                                   levels = 1:3,
                                   labels = paste0("State ", 1:3))]
  stprobs[, c("sample", "strategy_id", "state_id") := NULL]
  stprobs[, lab := "hesim"]
  setnames(stprobs, "t", "time")
  return(stprobs)
}

## comparison plot
plot_comparison1 <- function(fit, n_patients = 1000, clock){
  t_grid <- seq(0, max(fit_data$Tstop), .01)
  mstate_stprobs <- sim_stprobs_mstate(fit, n_patients, t_grid, clock)
  hesim_stprobs <- sim_stprobs_hesim(fit, n_patients, t_grid, clock)
  
  # plot
  pdat <- rbind(mstate_stprobs, hesim_stprobs)
  p <- ggplot(pdat, aes(x = time, y = prob, col = lab)) +
       geom_line() + 
       facet_wrap(~state) +
    xlab("Years") + ylab("Probability in health state") +
    scale_color_discrete(name = "") + theme_minimal() +
    theme(legend.position = "bottom") 
  return(p)
}

plot_comparisons <- function(fits, n_patients, clock,
                             filename){
  pdf(paste0("figs/", filename))
  for (i in 1:length(fits)){
    p <- plot_comparison1(fits[[i]], n_patients, clock)
    p <- p + labs(title = names(fits)[[i]])
    print(p)
    print(paste0("Completed plot for ", names(fits)[[i]], " model."))
  }
  dev.off()
}
plot_comparisons(fits_reset, 1000, "reset", 
                 "stprobs-ictstm-reset-hesim-vs-mstate.pdf")
plot_comparisons(fits_forward, 10000, "forward", 
                 "stprobs-ictstm-forward-hesim-vs-mstate.pdf") 

# State probabilities: fractional polynomial vs. weibull -----------------------
t_grid <- seq(0, max(fit_data$Tstop), .01)

# Weibull NMA fit
weiNMA_fit <- flexsurvreg(Surv(years, status) ~ trans, data = fit_data, 
                         dist = hesim_survdists$weibullNMA)
weiNMA_params <- create_params(weiNMA_fit, point_estimate = TRUE)

# Equivalent fractional polynomial parameters
fp_params <- weiNMA_params  
fp_params$dist <- "fracpoly"
fp_params$aux <- list(powers = c(0, 0),
                      cumhaz_method = "riemann", 
                      step = .02,
                      random_method = "sample")
names(fp_params$coefs) <- c("gamma0", "gamma1")
colnames(fp_params$coefs$gamma0)[1] <- "gamma0"
colnames(fp_params$coefs$gamma1)[1] <- "gamma1"

# Simulate
sim_stprobs_fp <- function(obj, param_names, n_patients, mod_name){
  transitions <- create_trans_dt(tmat)
  hesim_dat <- hesim_data(strategies <- data.table(strategy_id = 1),
                        patients <- data.table(patient_id = 1:n_patients),
                        transitions = transitions)
  input_dat <- expand(hesim_dat, by = c("strategies", "patients", "transitions")) 
  input_dat[, trans2 := 1 * (transition_id == 2)]
  input_dat[, trans3 := 1 * (transition_id == 3)]
  input_dat[, trans4 := 1 * (transition_id == 4)]
  input_dat[, (param_names) := 1]
  transmod <- create_IndivCtstmTrans(obj, 
                                    input_data = input_dat, trans_mat = tmat,
                                    clock = "reset") 
  stprobs <- transmod$sim_stateprobs(t = t_grid)
  stprobs[, lab := mod_name]
  return(stprobs)
}
weiNMA_stprobs <- sim_stprobs_fp(weiNMA_params, c("a0", "a1"), 1000, 
                                 "Weibull")
fp_stprobs <- sim_stprobs_fp(fp_params, c("gamma0", "gamma1"), 1000, 
                             "Fractional polynomial")
pdat <- rbind(weiNMA_stprobs, fp_stprobs)
p <- ggplot(pdat, aes(x = t, y = prob, col = lab)) +
            geom_line() + 
            facet_wrap(~state_id) +
            xlab("Years") + ylab("Probability in health state") +
            scale_color_discrete(name = "") + theme_minimal() +
            theme(legend.position = "bottom") 
ggsave("figs/stprobs-reset-fracpoly.pdf", p, width = 5, height = 7)

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
fp <- new(FracPoly, gamma = gamma, powers = powers,
          cumhaz_method = "quad", step = .01, random_method = "riemann")
fp$max_x_ <- 40
lower <- 0
upper <- fp$max_x_
step <- 1/12
time <- seq(lower, upper, step)

## Random sample with hesim
r1 <- replicate(10000, fp$random())
fun <- ecdf(r1)
esurv <- 1 - fun(time)
dat1 <- data.frame(time = time, surv = esurv, lab = "Random")

## Analytically compute survival
dat2 <- compute_surv(step = step, lower = lower, upper = upper, 
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
