# Simulate data for 4-state partitioned survival models
n_strategies <- 3

# Survival curves for 4-state partitioned survival model  ----------------------
# Simulate a 4-state partitioned survival model using a mutli-state framework
# The following transitions are possible
# 1 -> 2, 1 -> 4
# 2 -> 3, 2 -> 4
# 3 -> 4

sim_psm4_survival <- function(n_patients = 500) {
  set.seed(101)
  n_strategies <- 3
  n_obs <- n_strategies * n_patients

  # Parameters
  x_names <- c("intercept", "strategy_id2", "strategy_id3", "female", "age")
  beta <- rbind(
    c(log(1 / 2), log(.8), log(.65), log(.9), log(1.01)), # 1 -> 2
    c(log(1 / 10), log(.9), log(.8), log(.9), log(1.01)), # 1 -> 4
    c(log(1 / 2), log(.8), log(.65), log(.9), log(1.01)), # 2 -> 3
    c(log(1 / 6), log(.85), log(.75), log(.9), log(1.01)), # 2 -> 4
    c(log(1 / 2), log(.85), log(.75), log(.9), log(1.01)) # 3 -> 4
  )
  colnames(beta) <- x_names

  # Data
  ## Dataset
  female <- rbinom(n_obs, 1, .5)
  age <- rlnorm(n_obs, meanlog = 4, sdlog = .25)
  strategy_id <- factor(rep(seq_len(n_strategies), n_patients))
  sim_data <- data.frame(female, age, strategy_id = strategy_id)

  ## Design matrix
  n_trans <- nrow(beta)
  x <- model.matrix(~ strategy_id + female + age, sim_data)

  # Simulate
  censored_time <- rexp(n_obs, 1 / 10)

  ## General function to simulate competing risks
  sim_crisk <- function(x, beta, died = rep(FALSE, nrow(x)),
                        start_status = rep(1, nrow(x)), start_time = 0) {
    latent_time_k <- matrix(0, nrow = nrow(x), ncol = nrow(beta))
    for (i in 1:ncol(latent_time_k)) {
      latent_time_k[, i] <- rexp(nrow(x), exp(x %*% beta[i, ]))
    }
    latent_time <- ifelse(died | start_status == 0,
      start_time,
      apply(latent_time_k, 1, min) + start_time
    )
    time <- pmin(latent_time, censored_time)
    status <- ifelse(latent_time <= censored_time & start_status == 1, 1, 0)
    to <- apply(latent_time_k, 1, which.min)
    to[status == 0] <- NA
    to[died] <- ncol(latent_time_k)

    return(data.frame(time, status, to))
  }


  ## Transitions from state 1
  t1 <- sim_crisk(x, beta[1:2, ])
  sim_data$endpoint1_time <- t1$time
  sim_data$endpoint1_status <- t1$status

  ## Transitions from state 2
  t2 <- sim_crisk(x, beta[3:4, ],
    died = t1$to == 2,
    start_status = t1$status,
    start_time = t1$time
  )
  sim_data$endpoint2_time <- t2$time
  sim_data$endpoint2_status <- t2$status

  ## Transitions from state 3
  t3 <- sim_crisk(x, beta[5, , drop = FALSE],
    died = t2$to == 2,
    start_status = t2$status,
    start_time = t2$time
  )
  sim_data$endpoint3_time <- t3$time
  sim_data$endpoint3_status <- t3$status

  # Return
  return(sim_data)
}

survival_simdata <- sim_psm4_survival(n_patients = 500)

# Costs in 3 health states  ----------------------------------------------------
# Medical costs
sim_cost3_medical <- function(n_patients) {
  set.seed(101)

  # Parameters
  beta_intercept <- 30000
  beta_female <- 1000
  beta_state2 <- -5000
  beta_state3 <- 10000
  beta <- matrix(c(beta_intercept, beta_female, beta_state2, beta_state3),
    ncol = 1
  )

  # Data
  patients_df <- data.frame(
    patient_id = seq_len(n_patients),
    female = rbinom(n_patients, 1, .5)
  )
  state <- data.frame(state_name = factor(paste0("state", seq(1, 3))))
  data <- merge(patients_df, state)

  # Simulate
  X <- model.matrix(~ female + state_name, data)
  mean <- X %*% beta
  shape <- 10
  rate <- shape / mean
  y <- rgamma(nrow(X), rate = rate, shape = shape)
  data$costs <- y

  # Return
  return(data)
}

cost.medical.simdata <- sim_cost3_medical(n_patients = 30)

# Drug costs
sim_cost3_drugs <- function() {
  costs <- c(120000, 130000, 200000)
  costs.df <- data.frame(
    strategy_id = seq(1, 3),
    costs = costs
  )
  return(costs.df)
}

# Combine costs
costs_simdata <- list(
  medical = cost.medical.simdata,
  drugs = sim_cost3_drugs()
)

# Utility in 3 health states  --------------------------------------------------
utility <- data.frame(
  state = paste0("state", seq(1, 3)),
  lower = c(.8, .7, .6),
  upper = c(.9, .8, .7)
)

# Save -------------------------------------------------------------------------
psm4_exdata <- list(
  survival = survival_simdata,
  costs = costs_simdata,
  utility = utility
)
save(psm4_exdata, file = "../data/psm4_exdata.rda", compress = "bzip2")
