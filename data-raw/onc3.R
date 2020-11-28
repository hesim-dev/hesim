# Data for a 3-state (Stable, Progression, Death) oncology model
rm(list = ls())

# Simulate multi-state dataset -------------------------------------------------
sim_onc3_data <- function(n = 2000, seed = NULL){
  
  if (!is.null(seed)) set.seed(seed)
  
  # Data  
  data <- data.table(
    intercept = 1,
    strategy_id = 1,
    patient_id = 1:n,
    female = rbinom(n, 1, .5),
    new = rbinom(n, 1, .5),
    age = rnorm(n, mean = 60, sd = 5.5)
  )
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
    return(mean/(gamma(1 + 1/shape)))
  }
  
  matrixv <- function(v) {
    x <- matrix(v); colnames(x) <- "intercept"
    return(x)
  }
  
  params_wei <- function(shape, mean, 
                         beta_new = log(.6), 
                         beta_female = log(1.4)){
    log_shape <- matrixv(log(shape))
    scale = get_scale(shape, mean)
    beta_intercept <- log(scale) - beta_new
    scale_coefs <-  matrix(c(beta_intercept, beta_new, beta_female), 
                           ncol = 3)
    colnames(scale_coefs) <- c("intercept", "new", "female")
    params_surv(coefs = list(shape = log_shape,
                             scale = scale_coefs),
                dist = "weibull")
  }
  
  mstate_params <- params_surv_list(
    
    # 1. H -> S
    params_wei(shape = 2, mean = 1/.16),
    
    # 2. H -> D
    params_wei(shape = 3, mean = 10),
    
    # 3. S -> D
    params_wei(shape = 3.5, mean = 1/.12, beta_new = log(1))
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
  rc <- data.table(patient_id = 1:2000,
                   time = stats::rexp(n, rate = 1/15))
  sim[, time_rc := rc[match(sim$patient_id, rc$patient_id)]$time]
  sim[, status := ifelse(time_stop < 15 & time_stop < time_rc, status, 0)]
  sim[, time_stop := pmin(time_stop, 15, time_rc)]
  sim <- sim[time_start <= pmin(time_rc, 15)]
  
  ## Final data cleaning
  sim[, strategy_id := ifelse(new == 0, 1, 2)]
  sim[, strategy_name := factor(strategy_id, 
                                levels = c(1, 2),
                                labels = c("SOC", "New"))]
  label_states <- function (x) {
    fcase(
      x == 1, "Stable",
      x == 2, "Progression",
      x == 3, "Death"
    )
  }
  sim[, from := label_states(from)]
  sim[, to := label_states(to)]
  sim[, new := NULL]
  sim[, final := NULL]
  sim[, time_rc := NULL]
  
  # Return
  sim[, time := time_stop - time_start]
  return(sim[, ])
}
onc3 <- sim_onc3_data(seed = 101)

# Save
save(onc3, file = "../data/onc3.rda", compress = "bzip2")