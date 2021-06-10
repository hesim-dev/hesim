context("Model definition unit tests")

# Random number generation -----------------------------------------------------
alpha <- matrix(c(75, 25, 33, 67), byrow = TRUE, ncol = 2)
colnames(alpha) <- rownames(alpha) <- c("A", "B")
params <- list(
  alpha = alpha,
  gamma_mean = c(A = 900, B = 1500, C = 2000),
  normal_mean = c(A = 5, B = 10, C = 3),
  unif_min= 2,
  unif_max = 3
)

rng_def <- define_rng({
  beta_mean <- c(.2, .8, 0)
  beta_sd <- c(0, .1, 0)
  list( 
    dir1 = dirichlet_rng(alpha), 
    dir2 = dirichlet_rng(unname(alpha)),
    dir3 = dirichlet_rng(unname(alpha), names = alpha_names),
    gamma = gamma_rng(mean = gamma_mean,
                      sd = gamma_mean), 
    beta = beta_rng(mean = beta_mean,
                   sd = beta_sd,
                   names = beta_colnames),
    unif_vec = uniform_rng(min = unif_min, max = unif_max),
    unif_dt = uniform_rng(min = c(2, 3), max = c(4, 5)),
    normal = normal_rng(mean = normal_mean,
                       sd = rep(0, 3)),
    mvnorm = multi_normal_rng(mu = c(5, 10),
                             Sigma = matrix(c(10, 3, 3, 2) , nrow = 2, ncol = 2)),
    fixed_vec = fixed(est = 2),
    fixed_dt = fixed(est = c(2, 3)),
    custom_vec = custom(x = c(1, 2, 3)),
    custom_dt = custom(x = matrix(1:4, nrow = 2, ncol = 2))
    
  )
}, n = 2, beta_colnames = c("A", "B", "C"), 
   alpha_names = tpmatrix_names(c("A", "B"), prefix = "")) 
rng <- eval_rng(x = rng_def, params)

test_that( "eval_rng() runs without error", {
  expect_true(inherits(rng, "list"))
})
 
test_that( "eval_rng has correct number of samples", {
  n <- sapply(rng, function (z) if (length(dim(z)) < 2) length(z) else nrow(z))
  expect_true(all(n == rng_def$n))
})

test_that( "dirichlet_rng has correct column names", {
  expect_equal(colnames(rng$dir1), rng_def$alpha_names)
  expect_equal(colnames(rng$dir2), tpmatrix_names(c("s1", "s2"), prefix = ""))
  expect_equal(colnames(rng$dir3), rng_def$alpha_names)
})

test_that( "uv_rng() produces correct random variates for each column", {
  expect_equal(c(as.matrix(rng$normal)), 
               as.numeric(rep(params$normal_mean, each = rng_def$n)))
})

test_that( "uv_rng() produces correct vector or data.table", {
  expect_true(inherits(rng$unif_vec, "numeric"))
  expect_true(inherits(rng$unif_dt, "data.table"))
})

test_that( "mom_fun_rng produces fixed samples is sd == 0", {
  expect_true(all(rng$beta$A == .2))
  expect_true(all(rng$beta$C == 0))
  expect_equal(which(colnames(rng$beta) == "A"), 1)
  expect_equal(which(colnames(rng$beta) == "C"), 3)
})

test_that( "fixed() produces vector or data.table", {
  expect_true(inherits(rng$fixed_vec, "numeric"))
  expect_true(inherits(rng$fixed_dt, "data.table"))
})

test_that( "custom() produces vector or data.table", {
  expect_true(inherits(rng$custom_vec, "numeric"))
  expect_true(inherits(rng$custom_dt, "data.table"))
})

test_that( "custom() produces correct column names", {
  rng_def <- define_rng({
    custom3 <- custom(x = matrix(1:4, ncol = 2, nrow = 2))
    names(custom3) <- c("A", "B")
    list(
      custom1 = custom(x = matrix(1:4, ncol = 2, nrow = 2)),
      custom2 = custom(x = matrix(1:4, ncol = 2, nrow = 2), 
                       names = c("c1", "c2")),
      custom3 = custom3
    )
  }, n = 2)
  rng <- eval_rng(rng_def)
  expect_equal(colnames(rng$custom1), c("v1", "v2"))
  expect_equal(colnames(rng$custom2), c("c1", "c2"))
  expect_equal(colnames(rng$custom3), c("A", "B"))
})

test_that("custom() produces a warning if n > n_samples", {
  rng_def <- define_rng({
    list(
      custom = custom(x = 1)
    )
  }, n = 2)
  expect_warning(eval_rng(rng_def),
                 paste0("The number of requested draws for the probabilistic ",
                        "sensitivity analysis (PSA), 'n', is larger than the number ",
                        "of previously sampled values from the probability ",
                        "distribution of interest. Samples for the PSA have ",
                        "consequently been drawn with replacement."), fixed = TRUE)
})

test_that( "Column names for multi-parameter RNG is as expcted", {
  expect_equal(colnames(rng$gamma), names(params$gamma_mean))
  expect_equal(colnames(rng$beta), rng_def$beta_colnames)
  expect_equal(colnames(rng$mvnorm), c("v1", "v2"))
})

test_that( "multi_normal_rng() returns correct output when n = 1", {
  fun <- function(n, m = 0, V = 1){
    define_rng({ 
      x = multi_normal_rng(mu = m, Sigma = V)
    }, n = n, m = m, V = V)
  }
  expect_true(inherits(eval_rng(fun(1)), "numeric"))
  expect_true(
    inherits(eval_rng(fun(1, m = c(0, 0), V =  matrix(c(10,3,3,2),2,2))),
             "data.table")
  )
})

test_that( "define_rng() must return a list", {
  rng_def <- define_rng({
    data.frame(2)
  })
  expect_error(eval_rng(rng_def, check = TRUE),
               "define_rng() must return a list", fixed = TRUE)
})

test_that( "define_rng() has incorrect number of samples", {
  rng_def <- define_rng({
    list(
      x = gamma_rng(mean = c(10, 10),
                    sd = c(1, 1)),
      y = c(1, 2)
    )
  }, n = 3)
  expect_error(eval_rng(rng_def, check = TRUE),
               paste0("The number of samples produced by define_rng() must be ",
               "equal to n unless a scalar (of length 1) is returned."), 
               fixed = TRUE)
})

# Model definition works as expected -------------------------------------------
# Setup an example
## Data
strategies <- data.table(strategy_id = 1:2,
                         strategy_name = c("S1", "S2"))
patients <- data.table(patient_id = 1)
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients)
data <- expand(hesim_dat)

## Random number generation
rng_def <- define_rng({
  alpha <- matrix(c(1251, 350, 116, 17,
                    0, 731, 512, 15,
                    0, 0, 1312, 437,
                    0, 0, 0, 469),
                  nrow = 4, byrow = TRUE)
  rownames(alpha) <- colnames(alpha) <- c("A", "B", "C", "D")
  list(
    p_mono = dirichlet_rng(alpha),
    rr_s2 = .509,
    utility = 1,
    c_s1 = 2278,
    c_s2 = 4500,
    c_med = gamma_rng(mean = c(A = 2756, B = 3052, C = 9007),
                     sd = c(A = 2756, B = 3052, C = 9007))
  )
}, n = 2)

## Helper function
test_eval_model <- function(tparams_def, rng_def){
  model_def <- define_model(
    tparams_def = tparams_def,
    rng_def = rng_def
  )
  return(eval_model(model_def, data))
}

# Test
test_that("define_tparams() can be a list", {
  tparams_def1 <- define_tparams({
    list(
      costs = list(medical = c_med)
    )
  })
  tparams_def2 <- define_tparams({
    list(
      costs = list(drug = ifelse(strategy_name == "S1", c_s1, c_s2))
    )
  })  
  model_def <- define_model(
    tparams_def = list(tparams_def1, tparams_def2),
    rng_def = rng_def,
    n_states = 4
  )
  expect_equal(names(eval_model(model_def, data)$costs), c("medical", "drug"))
})

test_that("define_model() works with rng_def = NULL", {
  rng <- eval_rng(rng_def) 
  tparams_def <- define_tparams({
    list(tpmatrix = tpmatrix(1, 0, 0, 1))
  })
  model_def <- define_model(tparams_def, rng_def = NULL, params = rng)
  m <- eval_model(model_def, data)
  
  expect_equivalent(m$tpmatrix, tpmatrix(1, 0, 0, 1))
  rdata <- data[rep(1:nrow(data), rng_def$n)]
  expect_equal(m$id[[1]]$strategy_id, rdata$strategy_id)
  expect_equal(m$id[[1]]$patient_id, rdata$patient_id)
  expect_equal(m$id[[1]]$sample, rep(1:rng_def$n, each = nrow(data)))
})



# Model definition throws errors -----------------------------------------------
test_that( "define_rng() returns list elements of the right class", {
  error_msg <- paste0("Each element of the list returned by define_rng() must be a ",
                      "numeric vector, matrix, data.frame, or data.table.")
  
  # Error if element is a list
  rng_def <- define_rng({list(list(y = 3))})
  tparams_def <- define_tparams({
    list(z = 3)
  })
  expect_error(test_eval_model(tparams_def, rng_def),
               error_msg, fixed = TRUE)
  
  # Error if element is an array
  rng_def <- define_rng({list(y = array(1))})
  expect_error(test_eval_model(tparams_def, rng_def),
               error_msg, fixed = TRUE)
  
  # Should work with a matrix
  rng_def <- define_rng({list(y = matrix(1))})
  tparams_def <- define_tparams({
    list(tpmatrix = tpmatrix(1, 0, 0, 1))
  })
  mod <- test_eval_model(list(tparams_def), rng_def)
  expect_equivalent(mod$tpmatrix, tpmatrix(1, 0, 0, 1))
})

test_that("define_model() requires either params or rng_def to be non NULL", {
  expect_error(
    define_model(
      tparams_def = tparams_def, rng_def = NULL, params = NULL
    ),
    "'rng_def' and 'params' cannot both be NULL."
  )
})

test_that( "eval_tparams() must return a list", {
  tparams_def <- define_tparams({
    data.frame(2)
  })
  expect_error(test_eval_model(tparams_def, rng_def),
               "define_tparams() must return a list", fixed = TRUE)
})

test_that( "tpmatrix() must be square", {
  tparams_def <- define_tparams({
    list(
      tpmatrix = tpmatrix(
        C, p_mono$A_B, p_mono$A_C
      )
    )
  })
  expect_error(test_eval_model(tparams_def, rng_def),
               "tpmatrix() must be a square matrix", fixed = TRUE)
})

test_that("tpmatrix in define_tparams() must be square", {
  tparams_def <- define_tparams({
    list(
      tpmatrix = cbind(
        1 - p_mono$A_B - p_mono$A_C,
        p_mono$A_B, 
        p_mono$A_C
      )
    )
  })
  expect_error(test_eval_model(tparams_def, rng_def),
               "tpmatrix in define_tparams() must be a square matrix",
               fixed = TRUE)
})

test_that( "Number of states cannot be NULL", {
  tparams_def <- define_tparams({
    list(utility = utility)
  })
  expect_error(test_eval_model(tparams_def, rng_def),
               "'n_states' cannot be NULL.")
})

test_that( "Number of states in state value models must be correct", {
  tparams_def <- define_tparams({
    # Costs
    list(
      tpmatrix = tpmatrix(
        C, p_mono$A_B, p_mono$A_C,
        0, C, p_mono$A_D,
        0, 0, 1
      ),
      costs = list(medical = c_med),
      utility = utility
    )
  })
  expect_error(test_eval_model(tparams_def, rng_def),
               paste0("The number of columns in each element of 'costs' must ",
                      "equal 'n_states' - 1."),
               fixed = TRUE)
  
  # Utility
  tparams_def <- define_tparams({
    # Costs
    list(utility = cbind(1, 1)
    )
  })
  model_def <- define_model(
    tparams_def = tparams_def,
    rng_def = rng_def,
    n_states = 4
  )
  expect_error(eval_model(model_def, data),
               paste0("The number of columns in 'utility' must ",
                      "equal 'n_states' - 1."),
               fixed = TRUE)
  
})

test_that("costs in define_tparams() must be a list", {
  tparams_def <- define_tparams({
    list(
      costs = c_med
    )
  })
  model_def <- define_model(
    tparams_def = tparams_def,
    rng_def = rng_def,
    n_states = 3
  )
  expect_error(eval_model(model_def, data),
               "The 'costs' element returned by define_tparams() must be a list",
               fixed = TRUE)
})