# Simulate data for 4-state partitioned survival models
n_strategies <- 3

# Survival curves for 4-state partitioned survival model  ----------------------
sim_surv_data <-  function(X, beta){
  n.obs <- nrow(X)
  lograte <- X %*% beta
  latent.surv.times <- rexp(n.obs, exp(lograte))
  censoring.times <- rexp(n.obs, exp(-2))
  
  time <- pmin(latent.surv.times, censoring.times)
  status <- as.numeric(latent.surv.times <= censoring.times)
  return(data.frame(time = time, status = status))
}

sim_part_surv4_curves <- function(n_strategies = 3, n_patients){ # n_strategies
  set.seed(101)
  n_obs <- n_strategies * n_patients
  
  # Parameters
  beta1 <- c(0, .02, .01, -.1, 0)
  beta2 <- c(-.4, .02, .01, -.3, -.4)
  beta3 <- c(-.8, .02, .01, -.4, -.3)
  beta <- matrix(c(beta1, beta2, beta3), nrow = 3, byrow = TRUE) 
  
  # Data
  female <- rbinom(n_obs, 1, .5)
  age <- rlnorm(n_obs, meanlog = 4, sdlog = .25)
  strategy_id <- factor(rep(seq_len(n_strategies), n_patients))
  sim_data <- data.frame(female, age, strategy_id = strategy_id)
  
  # Simulate
  n_curves <- nrow(beta)
  X <- model.matrix(~female + age + strategy_id, sim_data)
  colnames(X)[colnames(X) == "(Intercept)"] <- "intercept"
  surv_data <- data.frame(female = female, age = age, strategy_id = strategy_id)
  censoring_times <- rexp(n_obs, exp(-2))
  prior_times <- rep(0, n_obs)
  
  for (i in 1:n_curves){
    latent_surv_times <- rexp(n_obs, exp(X %*% beta[i, ]))
    status_str <- paste0("endpoint", i, "_status")
    time_str <- paste0("endpoint", i, "_time")
    surv_data[[time_str]] <- pmin(latent_surv_times + prior_times, censoring_times) 
    surv_data[[status_str]] <- as.numeric(latent_surv_times + prior_times <= censoring_times)
  }
  return(surv_data)
}

survival_simdata <- sim_part_surv4_curves(n_patients = 500)
  

# Costs in 3 health states  ----------------------------------------------------
# Medical costs
sim_cost3_medical <- function(n_patients){
  set.seed(101)
  
  # Parameters
  beta_intercept <- 30000
  beta_female <- 1000
  beta_state2 <- -5000
  beta_state3 <- 10000
  beta <- matrix(c(beta_intercept, beta_female, beta_state2, beta_state3),
                 ncol = 1)
  
  # Data
  patients_df <- data.frame(patient_id = seq_len(n_patients),
                            female = rbinom(n_patients, 1, .5))
  state <- data.frame(state_name = factor(paste0("state", seq(1, 3))))
  data <- merge(patients_df, state)
  
  # Simulate
  X <- model.matrix(~female + state_name, data)
  mean <- X %*% beta
  shape <- 10
  rate <- shape/mean
  y <- rgamma(nrow(X), rate = rate, shape = shape)
  data$costs <- y
  
  # Return
  return(data)
}

cost.medical.simdata <- sim_cost3_medical(n_patients = 30)

# Drug costs
sim_cost3_drugs <- function(){
  costs <- c(120000, 130000, 200000)
  costs.df <- data.frame(strategy_id = seq(1, 3),
                      costs = costs)
  return(costs.df)
}

# Combine costs
costs_simdata <- list(medical = cost.medical.simdata,
                     drugs = sim_cost3_drugs())

# Utility in 3 health states  --------------------------------------------------
utility <- data.frame(state = paste0("state", seq(1, 3)),
                      lower = c(.8, .7, .6),
                      upper = c(.9, .8, .7))

# Save -------------------------------------------------------------------------
psm4_exdata <- list(survival = survival_simdata,
                    costs = costs_simdata,
                    utility = utility)
save(psm4_exdata, file = "../data/psm4_exdata.rda", compress = "bzip2")
