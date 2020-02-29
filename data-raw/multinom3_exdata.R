# Data for a 3-state discrete time state transition model
rm(list = ls())
library("data.table")
strategy_names <- c("Reference", "Intervention")
state_names <- c("Healthy", "Sick", "Dead")

# Simulate transitions based on multinomial logit model-------------------------
rlnorm_trunc <- function(n = 1000, meanlog = 0, sdlog = 1, lower = 18, 
                         upper = 100){
  u <- runif(n, 0, 1)
  cdf_lower <- plnorm(lower, meanlog, sdlog)
  cdf_upper <- plnorm(upper, meanlog, sdlog)
  v <- cdf_lower + (cdf_upper - cdf_lower) * u
  qlnorm(v, meanlog, sdlog)
}

sim_mlogit <- function(n = 10000){
  # Covariates
  age <- rlnorm_trunc(n, meanlog = log(55), sdlog = 1)
  female <- rbinom(n, 1, prob = .52)
  strategy_id <- sample.int(2, size = n, replace = TRUE)
  strategy_name <- factor(strategy_id,
                          levels = 1:2,
                          labels = strategy_names)
  age_cat <- factor(
    x = 1 * (age >= 40 &  age < 60) +
        2 * (age >= 60),
    labels = c("Age < 40",
               "40 <= Age < 60",
               "Age >= 60")
  )
  
  # Coefficients
  beta_h <- list(
    s = c(intercept = 0, treatment = log(.75),
          female = log(.85),
          age_40to60 = log(1.2), age_60plus = log(1.4),
          year3to6 = log(1.1), year7plus = log(1.2)),
    d = c(intercept = 0, treatment = log(.8),
          female = log(.80),
          age_40to60 = log(1.2), age_60plus = log(1.5),
          year3to6 = log(1.1), year7plus = log(1.2))
  )
  beta_s <- list(
    d = c(intercept = 0, treatment = log(.9),
          female = log(.75),
          age_40to60 = log(1.2), age_60plus = log(1.5),
          year3to6 = log(1.1), year7plus = log(1.2))
  )
  beta <- list(h = beta_h, s = beta_s)
  
  # Predicted probabilities
  predict_probs <- function(x, beta, state){
    z <- matrix(NA, nrow = nrow(x), ncol = 3)
    if (sum(state == "Healthy") > 0){
      x_healthy <- x[which(state == "Healthy"), ]
      z[which(state == "Healthy"), ] <- cbind(1, exp(x_healthy %*% beta$h$s), 
                                              exp(x_healthy %*% beta$h$d))
    }
    if (sum(state == "Sick") > 0){
      x_sick <- x[which(state == "Sick"), ]
      z[which(state == "Sick"), ] <- cbind(0, 1,  exp(x_sick %*% beta$s$d))
    }
    return(z/rowSums(z))
  }
  
  # Simulate
  sim_cat <- function(n, p){
    cat_id <- apply(p, 1, function (x) which(rmultinom(1, size = 1, x) == 1))
    cat <- factor(cat_id, levels = 1:3, labels = state_names)
    return(cat)
  }
  years <- seq(1, 12, 1)
  n_years <- length(years)
  sim <- vector(mode = "list", length = n_years)
  alive <- 1:n
  from_state <- rep("Healthy", n)
  for (t in 1:n_years){
    # Break out of loop if everyone has died
    if (length(alive) == 0){
      break
    }
    
    # Create data and design matrix
    ## Year dummies
    if(years[t] < 3){
      year_dummies <- c(0, 0)
    } else if (years[t] >= 3 & years[t] <= 6){
      year_dummies <- c(1, 0)
    } else{
      year_dummies <- c(0, 1)
    }
    year_dummies <- matrix(year_dummies, ncol = 2, nrow = length(alive),
                           byrow = TRUE)
    colnames(year_dummies) <- c("year3to6", "year7plus")

    ## Data
    data <- data.frame(patient_id = alive,
                       strategy_id = strategy_id[alive],
                       strategy_name = strategy_name[alive],
                       age = age[alive], age_cat = age_cat[alive],
                       female = female[alive], 
                       year = years[t],
                       year_dummies)
    x_t <- model.matrix(~strategy_id + female + age_cat + year3to6 + year7plus, data)
    
    # Simulate class categories
    p <- predict_probs(x_t, beta, from_state)
    to_states <- sim_cat(length(alive), p)
    
    # Store data
    sim[[t]] <- data.table(data, 
                           state_from = from_state, 
                           state_to = to_states)
    
    # Update variables for next loop iteration
    alive <- which(sim[[t]]$state_to != "Dead")
    from_state <- sim[[t]]$state_to[alive]
  }
  sim <- rbindlist(sim)
  
  # Return
  return(sim)
}
transitions <- sim_mlogit()
transitions[year == 1, .N, by = "state_to"]
transitions[, c("year3to6", "year7plus") := NULL]
transitions[, year_cat := factor(
  x =  1 * (year >= 3 & year <= 6) +
       2 * (year >= 7),
  labels = c("Year < 3",
             "3 <= Year <= 6",
             "Year >= 7"))]

# Costs ------------------------------------------------------------------------
costs <- list()

# Drugs
costs$drugs <- data.frame(strategy_id = c(1, 2),
                          strategy_name = strategy_names,
                          est = c(2000, 5000))

# Medical
costs$medical <- data.frame(state_id = c(1, 2),
                            state_name = state_names[1:2],
                            mean = c(3000, 1500),
                            se = sqrt(c(150, 100)))
# Utility ----------------------------------------------------------------------
utility <- data.frame(state_id = c(1, 2),
                      state_name = state_names[1:2],
                      mean = c(.90, .65),
                      se = sqrt(c(.03, .04)))

# Save -------------------------------------------------------------------------
multinom3_exdata <- list(transitions = data.frame(transitions),
                       costs = costs, 
                       utility = utility)
save(multinom3_exdata, file = "../data/multinom3_exdata.rda", compress = "bzip2")