format_run_time <- function(x) {
  fcase(
    x < 60, paste0(formatC(x, digits = 2), " seconds"),
    x > 60 & x < 3600, paste0(formatC(x/60, digits = 2), " minutes"),
    x > 3600, paste0(formatC(x/3600, digits = 2), " hours")
  )
}

# Semi-Markov model ------------------------------------------------------------
init_spline_weibull <- function(formula, data, k = 0) {
  # Parameter estimates from Weibull model
  fit_wei <- flexsurv::flexsurvreg(formula, data = data, dist = "weibullPH")
  wei_coef <- fit_wei$res.t[, "est"]
  
  # Initialize spline parameter vector
  inits <- rep(0, 2 + fit_wei$ncoveffs)
  cov_names <- names(coef(fit_wei))[fit_wei$covpars]
  names(inits) <-  c("gamma0", "gamma1", cov_names)
  
  # Set spline parameter values
  inits["gamma0"] <- wei_coef["scale"]
  inits["gamma1"] <- exp(wei_coef["shape"])
  inits[cov_names] <- wei_coef[cov_names]
  return(inits)
}

fit_semi_markov <- function(dist) {
  # Data
  tmat <- mstate::transMat(x = list(c(2, 3, 5, 6), 
                                    c(4, 5, 6), 
                                    c(4, 5, 6), 
                                    c(5, 6),
                                    c(),
                                    c()),
                           names = c("Tx", "Rec", "AE", "Rec+AE", 
                                     "Rel", "Death"))
  data("ebmt4", package = "mstate")
  msebmt <- mstate::msprep(data = ebmt4, trans = tmat, 
                           time = c(NA, "rec", "ae","recae", "rel", "srv"), 
                           status = c(NA, "rec.s", "ae.s", "recae.s", "rel.s", "srv.s"), 
                           keep = c("match", "proph", "year", "agecl"))
  
  # Model fitting
  n_trans <- max(tmat, na.rm = TRUE)
  fits <- fits_spline <- vector(mode = "list", length = n_trans)
  msebmt$years <- msebmt$time/365.25
  f <- formula(survival::Surv(years, status) ~ match + proph + year + agecl)
  for (i in 1:n_trans){
    data_i <- msebmt[msebmt$trans == i, ]
    fits[[i]] <- flexsurv::flexsurvreg(f, data = data_i, dist = dist)
    inits <- init_spline_weibull(f, data = data_i) 
    fits_spline[[i]] <- flexsurv::flexsurvspline(f, data = data_i, scale = "hazard",
                                                 k = 0, inits = inits)
  }
  return(list(fit = fits, fit_spline = fits_spline, data = msebmt))
}

make_semi_markov_input_data <- function(newdata, n_rep) {
  patients <- newdata[rep(1, n_rep), ]
  patients$patient_id <- 1:n_rep
  strategies <- data.frame(strategy_id = 1)  
  hesim_dat <- hesim_data(strategies = strategies,
                          patients = patients)
  input_data <- hesim::expand(hesim_dat, by = c("strategies", "patients")) 
  return(input_data)
}

sim_semi_markov_mstate <- function(object, newdata, n_rep = 100, n_samples = 100,
                                   uncertainty = c("normal", "none"),
                                   step = .01, yr_grid = seq(0, 10, .1)) {
  ptm <- proc.time()
  dat <- object$data
  tmat <- attr(dat, "trans")
  cumhaz_grid <- seq(0, max(dat$years), step)
  
  if (uncertainty == "none") {
    cumhaz <- flexsurv::msfit.flexsurvreg(object$fit, trans = tmat, 
                                          t = cumhaz_grid,
                                          newdata = newdata,
                                          variance = FALSE)
    sim <- mstate::mssample(Haz = cumhaz$Haz, 
                            trans = tmat,
                            tvec = yr_grid,
                            clock = "reset",
                            M = n_rep) 
    sim$sample <- 1
    sim <- data.table::as.data.table(sim)
  } else {
    # Predict cumulative hazard
    input_dat <- make_semi_markov_input_data(newdata = newdata, n_rep = 1)
    dismod <- create_IndivCtstmTrans(hesim::flexsurvreg_list(object$fit), 
                                     input_data = input_dat,
                                     trans_mat = tmat,
                                     clock = "reset",
                                     n = n_samples,
                                     uncertainty = uncertainty)
    cumhaz <- dismod$cumhazard(t = cumhaz_grid)  
    data.table::setnames(cumhaz, 
                         c("transition_id", "t", "cumhazard"), 
                         c("trans", "time", "Haz"))
    
    # Simulate
    sim <- vector(mode = "list", length = n_samples)
    for (s in 1:n_samples){
      sim[[s]] <- mstate::mssample(Haz = cumhaz[sample == s], 
                                   trans = tmat,
                                   tvec = yr_grid,
                                   clock = "reset",
                                   M = n_rep) 
    }
    sim <- data.table::rbindlist(sim, idcol = "sample")
  }
  
  # Return
  run_time <- proc.time() - ptm
  attr(sim, "run_time") <-  run_time
  return(sim)
}

sim_semi_markov_hesim <- function(object, newdata, n_rep = 100, n_samples = 100,
                                   uncertainty = c("normal", "none"),
                                   yr_grid = seq(0, 10, .1), spline = FALSE,
                                   step = NULL) {
  ptm <- proc.time()
  
  # Input data
  patients <- newdata[rep(1, n_rep), ]
  patients$patient_id <- 1:n_rep
  strategies <- data.frame(strategy_id = 1)  
  hesim_dat <- hesim_data(strategies = strategies,
                          patients = patients)
  input_data <- expand(hesim_dat, by = c("strategies", "patients")) 
  
  # Simulate
  fit <- if (!spline) object$fit else object$fit_spline
  dismod <- create_IndivCtstmTrans(hesim::flexsurvreg_list(fit), 
                                   input_data = input_data,
                                   trans_mat = attr(object$data, "trans"),
                                   clock = "reset",
                                   uncertainty = uncertainty,
                                   n = n_samples) 
  if (!is.null(step)) {
    dismod$params <- params_surv_list(lapply(dismod$params, function(z){
      z$aux$random_method <- "discrete"
      z$aux$cumhaz_method <- "riemann"
      z$aux$step <- step
      return(z)
    }))
  }
  sim <- dismod$sim_stateprobs(t = yr_grid)
  
  # Return
  run_time <- proc.time() - ptm
  attr(sim, "run_time") <-  run_time
  return(sim)
}

rbind_semi_markov <- function(mstate_sim, hesim_sim, hesim_sim_spline) {
  # mstate
  sim1 <- data.table::melt(mstate_sim, id.vars = c("sample", "time"), 
                           variable.name = "state_id",
                           value.name = "prob")
  sim1[, state_id := sub("pstate", "", state_id)]
  sim1[, state_id := as.integer(state_id)]
  sim1[, lab := "mstate"]
  
  # hesim
  ## Parametric model
  sim2a <- data.table::copy(hesim_sim)
  sim2a[, lab := "hesim (parametric)"]
  
  ## Spline
  sim2b <- data.table::copy(hesim_sim_spline)
  sim2b[, lab := "hesim (spline)"]
  
  # Combine
  sim2 <- rbind(sim2a, sim2b)
  sim2[, c("strategy_id", "grp_id") := NULL]
  data.table::setnames(sim2, "t", "time")
  
  # Combine
  return(rbind(sim1, sim2))
}

plot_semi_markov <- function(x, state_labels) {
   x <- x[, prob := mean(prob), by = c("time", "state_id", "lab")]
   x[, state_name := factor(state_id, labels = state_labels)]
   ggplot2::ggplot(x, ggplot2::aes(x = time, y = prob, col = lab)) +
     ggplot2::geom_line() +
     ggplot2::facet_wrap(~state_name, scales = "free_y") +
     ggplot2::xlab("Years") + 
     ggplot2::ylab("Probability in health state") +
     ggplot2::scale_x_continuous(breaks = seq(0, max(x$t), 2)) +
     ggplot2::scale_color_discrete(name = "") +
     ggplot2::theme_bw() +
     ggplot2::theme(legend.position = "bottom") 
}

benchmark_semi_markov <- function(n_patients = 100, n_samples = 1, 
                                 uncertainty = c("normal", "none"),
                                 dist = "gompertz", step = .01,
                                 spline_step = FALSE) {
  fit <- fit_semi_markov(dist = dist)
  pat2 <- data.frame(fit$data[fit$data$id == 2, 
                              c("match", "proph", "year", "agecl")][1, ])
  yr_grid = seq(0, 10, .1)
  mstate_sim <- sim_semi_markov_mstate(fit, newdata = pat2, n_rep = n_patients,
                                        n_samples = n_samples, uncertainty = uncertainty,
                                        step = step, yr_grid = yr_grid) 
  hesim_sim <- sim_semi_markov_hesim(fit, newdata = pat2, n_rep = n_patients,
                                     n_samples = n_samples, uncertainty = uncertainty,
                                     yr_grid = yr_grid)
  spline_step <-  if (spline_step) step else NULL
  hesim_sim_spline <- sim_semi_markov_hesim(fit, newdata = pat2, n_rep = n_patients,
                                            n_samples = n_samples, uncertainty = uncertainty,
                                            yr_grid = yr_grid, spline = TRUE,
                                            step = spline_step)
  sim <- rbind_semi_markov(mstate_sim = mstate_sim, hesim_sim = hesim_sim,
                           hesim_sim_spline = hesim_sim_spline)
  plot <- plot_semi_markov(sim, state_labels = rownames(attr(fit$data, "trans")))
  if (uncertainty == "none") n_samples <- 1
  return(list(
    sim = sim,
    run_time = c(
      mstate = attr(mstate_sim, "run_time")[["elapsed"]],
      hesim = attr(hesim_sim, "run_time")[["elapsed"]],
      hesim_spline = attr(hesim_sim_spline, "run_time")[["elapsed"]]
    ),
    plot = plot,
    n_patients = n_patients,
    uncertainty = uncertainty,
    n_samples = n_samples
  ))
}

semi_markov_table <- function(x) {
  run_times <- lapply(x, function (z) format_run_time(z$run_time))
  
  make_row <- function(x, y) {
    data.table(x$n_patients, x$n_samples, matrix(y, nrow = 1))
  }
  
  mapply(make_row, x, run_times, SIMPLIFY = FALSE) %>%
    rbindlist() %>%
    setnames(new = c("# of patients",  "# of PSA samples", "mstate", 
                     "hesim (parametric)", "hesim (spline)")) %>%
    kable() %>%
    kable_styling() %>%
    add_header_above(c(" " = 2, "Run time" = 3)) 
}

# Markov model -----------------------------------------------------------------
sim_markov_heemod <- function(n_samples) {
  ptm <- proc.time()
  
  # Define parameters
  param <- define_parameters(
    age_init = 60,
    sex = 0,
    
    ## age increases with cycles
    age = age_init + markov_cycle,
    
    ## operative mortality rates
    omrPTHR = .02,
    omrRTHR = .02,
    
    ## re-revision mortality rate
    rrr = .04,
    
    ## parameters for calculating primary revision rate
    cons = -5.49094,
    ageC = -.0367,
    maleC = .768536,
    lambda = exp(cons + ageC * age_init + maleC * sex),
    log_gamma = 0.3740968,
    gamma = exp(log_gamma),
    
    log_rrNP1 = -1.344473,
    rrNP1 = exp(log_rrNP1),
    
    ## revision probability of primary procedure
    standardRR = 1 - exp(lambda * ((markov_cycle - 1) ^ gamma -
                                     markov_cycle ^ gamma)),
    np1RR = 1 - exp(lambda * rrNP1 * ((markov_cycle - 1) ^ gamma - 
                                        markov_cycle ^ gamma)),
    
    ## age-related mortality rate
    sex_cat = ifelse(sex == 0, "FMLE", "MLE"),
    mr = get_who_mr(age, sex_cat, country = "GBR", local = TRUE),
    
    ## state values
    u_SuccessP = .85,
    u_RevisionTHR = .30,
    u_SuccessR = .75,
    c_RevisionTHR = 5294
  )
  
  # Define transitions 
  mat_standard <- define_transition(
    state_names = c(
      "PrimaryTHR",
      "SuccessP",
      "RevisionTHR",
      "SuccessR",
      "Death"
    ),
    0, C, 0,          0, omrPTHR,
    0, C, standardRR, 0, mr,
    0, 0, 0,          C, omrRTHR+mr,
    0, 0, rrr,        C, mr,
    0, 0, 0,          0, 1
  )
  
  mat_np1 <- define_transition(
    state_names = c(
      "PrimaryTHR",
      "SuccessP",
      "RevisionTHR",
      "SuccessR",
      "Death"
    ),
    0, C, 0,          0, omrPTHR,
    0, C, np1RR,      0, mr,
    0, 0, 0,          C, omrRTHR+mr,
    0, 0, rrr,        C, mr,
    0, 0, 0,          0, 1
  )
  
  # Define strategies
  strat_standard <- define_strategy(
    transition = mat_standard,
    PrimaryTHR = define_state(
      utility = 0,
      cost = 0
    ),
    SuccessP = define_state(
      utility = discount(u_SuccessP, .015),
      cost = 0
    ),
    RevisionTHR = define_state(
      utility = discount(u_RevisionTHR, .015),
      cost = discount(c_RevisionTHR, .06)
    ),
    SuccessR = define_state(
      utility = discount(u_SuccessR, .015),
      cost = 0
    ),
    Death = define_state(
      utility = 0,
      cost = 0
    ),
    starting_values = define_starting_values(
      cost = 394
    )
  )
  
  strat_np1 <- define_strategy(
    transition = mat_np1,
    PrimaryTHR = define_state(
      utility = 0,
      cost = 0
    ),
    SuccessP = define_state(
      utility = discount(u_SuccessP, .015),
      cost = 0
    ),
    RevisionTHR = define_state(
      utility = discount(u_RevisionTHR, .015),
      cost = discount(c_RevisionTHR, .06)
    ),
    SuccessR = define_state(
      utility = discount(u_SuccessR, .015),
      cost = 0
    ),
    Death = define_state(
      utility = 0,
      cost = 0
    ),
    starting_values = define_starting_values(
      cost = 579
    )
  )
  
  
  # Run deterministic model 
  res_mod <- run_model(
    standard = strat_standard,
    np1      = strat_np1,
    parameters = param,
    cycles = 60,
    cost = cost,
    effect = utility,
    init = c(1L, 0, 0, 0, 0)
  )
  
  # Run PSA
  rr_coef <- c(0.3740968, -5.490935, -0.0367022, 0.768536, -1.344474)
  names(rr_coef) <- c("lngamma", "cons", "age", "male", "np1")
  
  rr_vcov <- matrix(
    c(0.0474501^2, -0.005691, 0.000000028, 0.0000051, 0.000259,
      -0.005691, 0.207892^2, -0.000783, -0.007247, -0.000642,
      0.000000028, -0.000783, 0.0052112^2, 0.000033, -0.000111,
      0.0000051, -0.007247, 0.000033, 0.109066^2, 0.000184,
      0.000259, -0.000642, -0.000111, 0.000184, 0.3825815^2),
    ncol = 5, nrow = 5, byrow = TRUE
  )
  
  rsp <- define_psa(
    omrPTHR ~ beta(shape1 = 2, shape2 = 98),
    omrRTHR ~ beta(shape1 = 2, shape2 = 98),
    rrr ~ beta(shape1 = 4, shape2 = 96),
    u_SuccessP ~ beta(shape1 = .85, shape2 = .03),
    u_RevisionTHR ~ beta(shape1 = .30, shape2 = .03),
    u_SuccessR ~ beta(shape1 = .75, shape2 = .04),
    c_RevisionTHR ~ gamma(mean = 5294, sd = sqrt(1487)),
    log_gamma ~ normal(0.3740968, 0.002251512),
    cons ~ normal(-5.49094, 0.04321908),
    ageC ~ normal(-.0367, 0.00002715661),
    maleC ~ normal(.768536, 0.01189539),
    log_rrNP1 ~ normal(-1.3444740, 0.1463686)
  )
  
  pm <- run_psa(
    model = res_mod,
    psa = rsp,
    N = n_samples
  )
  
  # 
  run_time <- proc.time() - ptm
  return(list(sim = pm, run_time = run_time))
}

sim_markov_hesim_cohort <- function(n_samples) {
  
  ptm <- proc.time()
  
  # Model setup
  strategies <- data.table(strategy_id = 1:2,
                           strategy_name = c("Standard prosthesis",
                                             "New prosthesis"))
  patients <- data.table(patient_id = 1,
                         sex = "Female",
                         age = 60)
  hesim_dat <- hesim_data(strategies = strategies,
                          patients = patients)
  
  # Parameters
  ## Estimates from literature
  ### Mortality
  mort_tbl <- rbind(
    c(35, 45, .00151, .00099),
    c(45, 55, .00393, .0026),
    c(55, 65, .0109, .0067),
    c(65, 75, .0316, .0193),
    c(75, 85, .0801, .0535),
    c(85, Inf, .1879, .1548)
  )
  colnames(mort_tbl) <- c("age_lower", "age_upper", "male", "female")
  mort_tbl <- data.frame(mort_tbl)
  
  ### Revision risk
  #### Coefficients
  rr_coef <- c(0.3740968, -5.490935, -0.0367022, 0.768536, -1.344474)
  names(rr_coef) <- c("lngamma", "cons", "age", "male", "np1")
  
  #### Variance-covariance matrix
  rr_vcov <- matrix(
    c(0.0474501^2, -0.005691, 0.000000028, 0.0000051, 0.000259,
      -0.005691, 0.207892^2, -0.000783, -0.007247, -0.000642,
      0.000000028, -0.000783, 0.0052112^2, 0.000033, -0.000111,
      0.0000051, -0.007247, 0.000033, 0.109066^2, 0.000184,
      0.000259, -0.000642, -0.000111, 0.000184, 0.3825815^2),
    ncol = 5, nrow = 5, byrow = TRUE
  )
  rownames(rr_vcov) <- colnames(rr_vcov) <- names(rr_coef)
  
  #### Combine all parameters
  params <- list(
    # Transition probabilities
    ## Operative mortality following primary THR
    omrPTHR_shape1 = 2, 
    omrPTHR_shape2 = 98,
    
    ## Revision rate for prosthesis
    rr_coef = rr_coef,
    rr_vcov = rr_vcov,
    
    ## Mortality_rates
    mr = mort_tbl,
    
    ## Operative mortality following revision THR
    omrRTHR_shape1 = 2,
    omrRTHR_shape2 = 98,
    
    ## re-revision rate
    rrr_shape1 = 4,
    rrr_shape2 = 96,
    
    # Utility
    u_mean = c(PrimaryTHR = 0, SuccessP = .85, Revision = .30, SuccessR = .75),
    u_se = c(PrimaryTHR = 0, SuccessP = .03, Revision = .03, SuccessR = .04),
    
    # Costs
    c_med_mean = c(PrimaryTHR = 0, SuccessP = 0, Revision = 5294, SuccessR = 0),
    c_med_se = c(PrimaryTHR = 0, SuccessP = 0, Revision = 1487, SuccessR = 0),
    c_Standard = 394,
    c_NP1 = 579
  )
  
  ### Random number generation
  rng_def <- define_rng({
    list( 
      omrPTHR = beta_rng(shape1 = omrPTHR_shape1, shape2 = omrPTHR_shape2),
      rr_coef = multi_normal_rng(mu = rr_coef, Sigma = rr_vcov),
      mr_male = fixed(mr$male, names = mr$age_lower),
      mr_female = fixed(mr$female, names = mr$age_lower),
      omrRTHR = beta_rng(shape1 = omrRTHR_shape1, shape2 = omrRTHR_shape2),
      rrr = beta_rng(shape1 = rrr_shape1, shape2 = rrr_shape2),
      u = beta_rng(mean = u_mean, sd = u_se),
      c_med = gamma_rng(mean = c_med_mean, sd = c_med_se),
      c_Standard = c_Standard,
      c_NP1 = c_NP1
    )
  }, n = n_samples)
  
  ### Transformed parameters (transition probability matrix)
  transitions_def <- define_tparams({
    #### Regression for revision risk
    male <- ifelse(sex == "Female", 0, 1)
    np1 <- ifelse(strategy_name == "Standard prosthesis", 0, 1)
    scale <- exp(rr_coef$cons + rr_coef$age * age + rr_coef$male * male + 
                   rr_coef$np1 * np1)
    shape <- exp(rr_coef$lngamma)
    rr <- 1 - exp(scale * ((time - 1)^shape - time^shape))
    
    #### Mortality rate
    age_new <- age + time
    mr <- mr_female[["35"]] * (sex == "Female" & age_new >= 35 & age_new < 45) +
      mr_female[["45"]] * (sex == "Female" & age_new >= 45 & age_new < 55) +
      mr_female[["55"]] * (sex == "Female" & age_new >= 55 & age_new < 65) +
      mr_female[["65"]] * (sex == "Female" & age_new >= 65 & age_new < 75) +
      mr_female[["75"]] * (sex == "Female" & age_new >= 75 & age_new < 85) +
      mr_female[["85"]] * (sex == "Female" & age_new >= 85) +
      
      mr_male[["35"]] * (sex == "Male" & age_new >= 35 & age_new < 45) +
      mr_male[["45"]] * (sex == "Male" & age_new >= 45 & age_new < 55) +
      mr_male[["55"]] * (sex == "Male" & age_new >= 55 & age_new < 65) +
      mr_male[["65"]] * (sex == "Male" & age_new >= 65 & age_new < 75) +
      mr_male[["75"]] * (sex == "Male" & age_new >= 75 & age_new < 85) +
      mr_male[["85"]] * (sex == "Male" & age_new >= 85)
    
    list(
      tpmatrix = tpmatrix(
        0, C, 0,   0, omrPTHR,
        0, C, rr,  0, mr,
        0, 0, 0,   C, omrRTHR + mr,
        0, 0, rrr, C, mr,
        0, 0, 0,   0, 1)
    )
  }, times = 1:60)
  
  statevals_def <- define_tparams({
    c_prosthesis <- ifelse(strategy_name == "Standard prosthesis",
                           c_Standard,
                           c_NP1)
    list(
      utility = u,
      costs = list(
        prosthesis = c_prosthesis,
        medical = c_med
      )
    )
  })
  
  # Simulation
  ## Construct model
  mod_def <- define_model(tparams_def = list(transitions_def, 
                                             statevals_def),
                          rng_def = rng_def, 
                          params = params)
  
  cost_args <- list(
    prosthesis = list(method = "starting"),
    medical = list(method = "wlos")
  )
  input_data <- expand(hesim_dat, by = c("strategies", "patients"))
  econmod <- create_CohortDtstm(mod_def, input_data,
                                cost_args = cost_args)
  
  # Simulation
  ## Simulate outcomes
  econmod$sim_stateprobs(n_cycles = 60)
  econmod$sim_qalys(dr = .015, integrate_method = "riemann_right")
  econmod$sim_costs(dr = .06, integrate_method = "riemann_right")
  ce_sim <- econmod$summarize()
  
  # Return
  run_time <- proc.time() - ptm
  return(list(sim = ce_sim, run_time = run_time))
}

sim_markov_hesim_indiv <- function(n_samples, n_patients) {
  ptm <- proc.time()
  
  # Model setup
  ## Treatment strategies
  strategies <- data.table(
    strategy_id = 1:2,
    strategy_name = c("Standard prosthesis", "New prosthesis")
  )
  n_strategies <- nrow(strategies)
  
  ## Patients
  patients <- data.table(
    patient_id = 1:n_patients,
    gender = "Female",
    age = 60
  )
  
  ## Health states
  states <- data.table(
    state_id = 1:4,
    state_name = c("PrimaryTHR", "SuccessP", "Revision", "SuccessR")
  ) # Non-death health states
  n_states <- nrow(states)
  
  ## Transitions
  tmat <- rbind(c(NA, 1, NA, NA, 2),
                c(NA, NA, 3, NA, 4),
                c(NA, NA, NA, 5, 6),
                c(NA, NA, 7, NA, 8),
                c(NA, NA, NA, NA, NA))
  colnames(tmat) <- rownames(tmat) <- c(states$state_name, "Death")
  
  ## "hesim data"
  hesim_dat <- hesim_data(strategies = strategies,
                          patients = patients, 
                          states = states)
  
  # Parameters
  ## Transitions
  ### Estimates from literature
  #### Revision risk
  rr_coef <- c(0.3740968, -5.490935, -0.0367022, 0.768536, -1.344474)
  names(rr_coef) <- c("lngamma", "cons", "age", "male", "np1")
  
  rr_vcov <- matrix(
    c(0.0474501^2, -0.005691, 0.000000028, 0.0000051, 0.000259,
      -0.005691, 0.207892^2, -0.000783, -0.007247, -0.000642,
      0.000000028, -0.000783, 0.0052112^2, 0.000033, -0.000111,
      0.0000051, -0.007247, 0.000033, 0.109066^2, 0.000184,
      0.000259, -0.000642, -0.000111, 0.000184, 0.3825815^2),
    ncol = 5, nrow = 5, byrow = TRUE
  )
  rownames(rr_vcov) <- colnames(rr_vcov) <- names(rr_coef)
  
  #### Storing the parameters
  params <- list(
    #### Transition 1
    ttrrPTHR = 2, # Time to recovery rate implies mean time of 1/2 years
    
    #### Transition 2
    omrPTHR_shape1 = 2, # 2 out of 100 patients receiving primary THR died
    omrPTHR_shape2 = 98,
    
    #### Transition 3
    rr_coef = rr_coef,
    rr_vcov = rr_vcov,
    
    #### Transition 4
    mr = c(.0067, .0193, .0535, .1548),
    
    #### Transition 5
    ttrRTHR = 1, # There is no rate, the time is fixed
    
    #### Transition 6: omrRTHR + mr
    
    #### Transition 7
    omrRTHR_shape1 = 4, # 4 out of 100 patients with a successful revision needed another procedure
    omrRTHR_shape2 = 96
    
    #### Transition 8: same as transition 4
  )
  
  ### Multi-state model
  matrixv <- function(v, n = NULL){
    if (length(v) == 1) v <- rep(v, n_samples) 
    m <- matrix(v)
    colnames(m) <- "cons"
    return(m)
  }
  
  prob_to_rate <- function(p, t = 1){
    (-log(1 - p))/t
  }
  
  transmod_coef_def <- define_rng({
    omrPTHR <- prob_to_rate(beta_rng(shape1 = omrPTHR_shape1,
                                     shape2 = omrPTHR_shape2))
    mr <- fixed(mr)
    mr_omrPTHR <- omrPTHR + mr
    rr <- multi_normal_rng(mu = rr_coef, Sigma = rr_vcov)
    rrr <- prob_to_rate(beta_rng(shape1 = 4, shape2 = 96))
    
    list(
      log_omrPTHR = matrixv(log(omrPTHR)),
      log_mr = lapply(as.list(log(mr)), matrixv),
      log_ttrrPTHR = matrixv(log(ttrrPTHR)),
      log_mr_omrPTHR = lapply(as.list(log(mr_omrPTHR)), matrixv),
      rr_shape = matrixv(rr$lngamma),
      rr_scale = as.matrix(rr[, -1,]),
      log_rrr = matrixv(log(rrr))
    )
  }, n = n_samples, prob_to_rate = prob_to_rate, matrixv = matrixv)
  transmod_coef <- eval_rng(transmod_coef_def, params = params)
  
  transmod_params <- params_surv_list(
    # 1. Primary THR:Successful primary (1:2)
    params_surv(coefs = list(rate = transmod_coef$log_ttrrPTHR), 
                dist = "fixed"),
    
    # 2. Primary THR:Death (1:5)
    params_surv(coefs = list(rate = transmod_coef$log_omrPTHR), 
                dist = "exp"), 
    
    # 3. Successful primary:Revision THR (2:3)
    params_surv(coefs = list(shape = transmod_coef$rr_shape,
                             scale = transmod_coef$rr_scale), 
                dist = "weibullPH"), 
    
    # 4. Successful primary:Death (2:5)
    params_surv(coefs = transmod_coef$log_mr,
                aux = list(time = c(0, 5, 15, 25)),
                dist = "pwexp"),
    
    # 5. Revision THR:Successful revision (3:4)
    params_surv(coefs = list(est = matrixv(params$ttrRTHR)),
                dist = "fixed"),
    
    # 6. Revision THR:Death (3:5)
    params_surv(coefs = transmod_coef$log_mr_omrPTHR,
                aux = list(time = c(0, 5, 15, 25)),
                dist = "pwexp"),
    
    # 7. Successful revision:Revision THR (4:3)
    params_surv(coefs = list(rate = transmod_coef$log_rrr),
                dist = "exp"), 
    
    # 8. Successful revision:Death (4:5)
    params_surv(coefs = transmod_coef$log_mr,
                aux = list(time = c(0, 5, 15, 25)),
                dist = "pwexp")
  )
  
  
  ## Utility and costs
  utility_tbl <- stateval_tbl(
    data.table(state_id = states$state_id,
               mean = c(0, .85, .3, .75),
               se = c(0, .03, .03, .04)),
    dist = "beta"
  )
  
  drugcost_tbl <- stateval_tbl(
    data.table(strategy_id = rep(strategies$strategy_id, each = n_states),
               state_id = rep(states$state_id, times = n_strategies),
               est = c(394, 0, 0, 0,
                       579, 0, 0, 0)),
    dist = "fixed"
  )
  
  medcost_tbl <- stateval_tbl(
    data.table(state_id = states$state_id,
               mean = c(0, 0, 5294, 0),
               se = c(0, 0, 1487, 0)),
    dist = "gamma"
  )
  
  # Simulation
  ## Construct model
  ### Transition model
  transmod_data <- expand(hesim_dat, by = c("strategies", "patients"))
  transmod_data[, cons := 1]
  transmod_data[, male := ifelse(gender == "Male", 1, 0)]
  transmod_data[, np1 := ifelse(strategy_name == "New prosthesis", 1, 0)]
  
  transmod <- create_IndivCtstmTrans(transmod_params, 
                                     input_data = transmod_data,
                                     trans_mat = tmat,
                                     clock = "forward",
                                     start_age = patients$age)
  
  ### Utility and cost models
  utilitymod <- create_StateVals(utility_tbl, n = transmod_coef_def$n,
                                 hesim_data = hesim_dat)
  drugcostmod <- create_StateVals(drugcost_tbl, n = transmod_coef_def$n,
                                  method = "starting",
                                  hesim_data = hesim_dat)
  medcostmod <- create_StateVals(medcost_tbl, n = transmod_coef_def$n,
                                 hesim_data = hesim_dat)
  costmods <- list(Drug = drugcostmod,
                   Medical = medcostmod)
  
  ### Economic model
  econmod <- IndivCtstm$new(trans_model = transmod,
                            utility_model = utilitymod,
                            cost_models = costmods)
  
  ## Simulate outcomes
  econmod$sim_disease(max_t = 60, max_age = 120)
  econmod$sim_qalys(dr = .015)
  econmod$sim_costs(dr = .06)
  ce_sim <- econmod$summarize()
  
  # Return
  run_time <- proc.time() - ptm
  return(list(sim = ce_sim, run_time = run_time))
}

benchmark_markov <- function(n_samples, n_patients = 1000) {
  # Simulate
  heemod_sim <- sim_markov_heemod(n_samples)
  hesim_cohort_sim <- sim_markov_hesim_cohort(n_samples)
  hesim_indiv_sim <- sim_markov_hesim_indiv(n_samples, n_patients = n_patients)
  
  # return
  return(list(
    run_time = c(heemod = heemod_sim$run_time[["elapsed"]],
                 hesim_cohort = hesim_cohort_sim$run_time[["elapsed"]],
                 hesim_indiv = hesim_indiv_sim$run_time[["elapsed"]]),
    n_samples = n_samples,
    n_patients = n_patients
  ))
}

markov_table <- function(x) {
  make_row <- function(x) {
    data.table(x$n_samples, 1, x$n_patients, x$run_time[1], x$run_time[2], x$run_time[3])
  }
  
  lapply(x, make_row) %>%
    rbindlist() %>%
    setnames(new = c("# of PSA samples", "Cohort", "Individual", 
                     "heemod", "hesim_cohort", "hesim_indiv")) %>%
    .[, heemod := format_run_time(heemod)] %>%
    .[, hesim_cohort := format_run_time(hesim_cohort)] %>%
    .[, hesim_indiv := format_run_time(hesim_indiv)] %>%
    setnames(old = c("hesim_cohort", "hesim_indiv"),
             new = c("hesim (cohort)", "hesim (individual)")) %>%
    kable() %>%
    kable_styling() %>%
    add_header_above(c(" " = 1, "# of patients" = 2, "Run time" = 3)) 
}