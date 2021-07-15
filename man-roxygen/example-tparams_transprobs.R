hesim_dat <- hesim_data(strategies = data.frame(strategy_id = 1:2),
                        patients = data.frame(patient_id = 1:3))
input_data <- expand(hesim_dat, by = c("strategies", "patients"))    

# tpmatrix objects provide a convenient way to construct
# tparams_transprobs() objects
tpmat_id <- tpmatrix_id(input_data, n_samples = 2)      
p_12 <- runif(nrow(tpmat_id), .6, .7) + 
  .05 * (tpmat_id$strategy_id == 2)
tpmat <- tpmatrix(
  C, p_12,
  0, 1
)
tprobs <- tparams_transprobs(tpmat, tpmat_id)
names(tprobs) # Names of list elements

# Convert to data.table in wide format
as.data.table(tprobs)

# Convert to data.table in long format
as.data.table(tprobs, long = TRUE)

# Summary where each column is a vector
summary(tprobs)

# Summary where each column is a matrix
ps <- summary(tprobs, id = tpmat_id, unflatten = TRUE)
ps
ps$mean