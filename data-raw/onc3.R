# Data for a 3-state (Stable, Progression, Death) oncology model
rm(list = ls())

# Simulate multi-state dataset -------------------------------------------------
sim_onc3_data <- function(n = 2500, seed = NULL){
  
  if (!is.null(seed)) set.seed(seed)
  
  # Data 
  age_mu <- 60
  data <- data.table(
    intercept = 1,
    strategy_id = 1,
    strategy_name = sample(c("SOC", "New 1", "New 2"), n, replace = TRUE,
                             prob = c(1/3, 1/3, 1/3)),
    patient_id = 1:n,
    female = rbinom(n, 1, .5),
    age = rnorm(n, mean = age_mu, sd = 5.5)
  )
  data[, `:=` (new1 = ifelse(strategy_name == "New 1", 1, 0),
               new2 = ifelse(strategy_name == "New 2", 1, 0))]
  attr(data, "id_vars") <- c("strategy_id", "patient_id")
  
  # Transition matrix
  tmat <- rbind(
    c(NA, 1,  2),
    c(NA,  NA, 3),
    c(NA, NA, NA)
  )
  trans_dt <- create_trans_dt(tmat)
  
  # Parameters for each transition
  get_scale <- function(shape, mean) {
    scale <- mean/(gamma(1 + 1/shape))
    scale_ph <- scale^{-shape}
    return(scale_ph)
  }
  
  matrixv <- function(v) {
    x <- matrix(v); colnames(x) <- "intercept"
    return(x)
  }
  
  params_wei <- function(shape, mean,
                         beta_new1 = log(1), 
                         beta_new2 = log(1),
                         beta_age, beta_female){
    log_shape <- matrixv(log(shape))
    scale = get_scale(shape, mean)
    beta_intercept <- log(scale) - mean(data$age) * beta_age
    scale_coefs <-  matrix(c(beta_intercept, beta_new1, beta_new2, 
                             beta_age, beta_female), 
                           ncol = 5)
    colnames(scale_coefs) <- c("intercept", "new1", "new2", "age", "female")
    params_surv(coefs = list(shape = log_shape,
                             scale = scale_coefs),
                dist = "weibullPH")
  }
  
  mstate_params <- params_surv_list(
    
    # 1. S -> P
    params_wei(shape = 2, mean = 6.25, 
               beta_new1 = log(.7), beta_new2 = log(.6),
               beta_female = log(1.4), beta_age = log(1.03)),
    
    # 2. S -> D
    params_wei(shape = 2.5, mean = 10,
               beta_new1 = log(.85), beta_new2 = log(.75),
               beta_female = log(1.2), beta_age = log(1.02)),
    
    # 3. P -> D
    params_wei(shape = 3.5, mean = 8, beta_new1 = log(1),
               beta_female = log(1.3), beta_age = log(1.02))
  )
  
  # Create multi-state model
  mstatemod <- create_IndivCtstmTrans(mstate_params, 
                                      input_data = data,
                                      trans_mat = tmat,
                                      clock = "reset",
                                      start_age = data$age)
  
  # Simulate data
  ## Observed "latent" transitions
  sim <- mstatemod$sim_disease(max_age = 100)
  sim[, c("sample", "grp_id", "strategy_id") := NULL]
  sim <- cbind(
    data[match(sim$patient_id, data$patient_id)][, patient_id := NULL],
    sim
  )
  sim[, ":=" (intercept = NULL, strategy_id = NULL, status = 1, added = 0)]
  
  ## Add all possible states for each transition
  ### Observed 1->2 add 1->3
  sim_13 <- sim[from == 1 & to == 2]
  sim_13[, ":=" (to = 3, status = 0, final = 0,  added = 1)]
  sim <- rbind(sim, sim_13)
  
  ### Observed 1->3 add 1->2
  sim_12 <- sim[from == 1 & to == 3 & added == 0]
  sim_12[, ":=" (to = 2, status = 0, final = 0, added = 1)]
  sim <- rbind(sim, sim_12)
  
  ### Sort and clean
  sim <- merge(sim, trans_dt, by = c("from", "to")) # Add transition ID
  setorderv(sim, c("patient_id", "from", "to"))
  sim[, added := NULL]
  
  ## Add right censoring
  rc <- data.table(patient_id = 1:n,
                   time = stats::rexp(n, rate = 1/15))
  sim[, time_rc := rc[match(sim$patient_id, rc$patient_id)]$time]
  sim[, status := ifelse(time_stop < 15 & time_stop < time_rc, status, 0)]
  sim[, time_stop := pmin(time_stop, 15, time_rc)]
  sim <- sim[time_start <= pmin(time_rc, 15)]
  
  ## Final data cleaning
  sim[, strategy_id := fcase(
    strategy_name == "SOC", 1L,
    strategy_name == "New 1", 2L,
    strategy_name == "New 2", 3L
  )]
  sim[, strategy_name := factor(strategy_id, 
                                levels = c(1, 2, 3),
                                labels = c("SOC", "New 1", "New 2"))]
  label_states <- function (x) {
    fcase(
      x == 1, "Stable",
      x == 2, "Progression",
      x == 3, "Death"
    )
  }
  sim[, from := label_states(from)]
  sim[, to := label_states(to)]
  sim[, c("new1", "new2", "final", "time_rc") := NULL]
  
  # Return
  sim[, time := time_stop - time_start]
  return(sim[, ])
}
onc3 <- sim_onc3_data(n = 3000, seed = 102)

# Check that coefficient estimates are consistent with "truth"
fit_weibull <- function(i) {
  flexsurvreg(Surv(time, status) ~ strategy_name + female + age,
              data = onc3, subset = (transition_id == i), dist = "weibullPH")
}
fit_weibull(1)
fit_weibull(2)
fit_weibull(3)

# Save
save(onc3, file = "../data/onc3.rda", compress = "bzip2")