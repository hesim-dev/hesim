# Simulate data for partitioned survival analysis with 3 treatment strategies
n.strategies <- 3

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
  n.obs <- n.strategies * n_patients
  
  # Parameters
  beta1 <- c(0, .02, .01, -.1, 0)
  beta2 <- c(-.4, .02, .01, -.3, -.4)
  beta3 <- c(-.8, .02, .01, -.4, -.3)
  beta <- matrix(c(beta1, beta2, beta3), nrow = 3, byrow = TRUE) 
  
  # Data
  female <- rbinom(n.obs, 1, .5)
  age <- rlnorm(n.obs, meanlog = 4, sdlog = .25)
  strategy.id <- factor(rep(seq_len(n_strategies), n_patients))
  sim.data <- data.frame(female, age, strategy_id = strategy.id)
  
  # Simulate
  n.curves <- nrow(beta)
  X <- model.matrix(~female + age + strategy_id, sim.data)
  colnames(X)[colnames(X) == "(Intercept)"] <- "intercept"
  surv.data <- vector(mode = "list", length = n.curves)
  for (i in 1:n.curves){
    surv.data[[i]] <- sim_surv_data(X, beta = beta[i, ])
    colnames(surv.data[[i]]) <- c(paste0("endpoint", i, "_time"),
                                  paste0("endpoint", i, "_status"))
  }
  sim.data <- cbind(sim.data, do.call("cbind", surv.data))
  return(sim.data)
}

survival.simdata <- sim_part_surv4_curves(n_patients = 500)
  

# Costs in 3 health states  ----------------------------------------------------
# Medical costs
sim_cost3_medical <- function(n_patients){
  set.seed(101)
  
  # Parameters
  beta.intercept <- 30000
  beta.female <- 1000
  beta.state2 <- -5000
  beta.state3 <- 10000
  beta <- matrix(c(beta.intercept, beta.female, beta.state2, beta.state3),
                 ncol = 1)
  
  # Data
  patients.df <- data.frame(patient_id = seq_len(n_patients),
                            female = rbinom(n_patients, 1, .5))
  state <- data.frame(state_name = factor(paste0("state", seq(1, 3))))
  data <- merge(patients.df, state)
  
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
costs.simdata <- list(medical = cost.medical.simdata,
                     drugs = sim_cost3_drugs())

# Sample Weibull NMA data  -----------------------------------------------------
# part_surv4_weibull_NMA <- function(){
#   
# }

# Save -------------------------------------------------------------------------
part_surv4_simdata <- list(survival = survival.simdata,
                          costs = costs.simdata)
save(part_surv4_simdata, file = "../data/part_surv4_simdata.rda", compress = "bzip2")
