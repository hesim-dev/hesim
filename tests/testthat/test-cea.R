context("Cost-effectiveness analysis")

# Output for testing ----------------------------------------------------------
nsims <- 1000

# cost
c <- vector(mode = "list", length = 6)
names(c) <- c("Arm 1, Grp 1", "Arm 1, Grp 2", "Arm 2, Grp 1",
              "Arm 2, Grp 2", "Arm 3, Grp 1", "Arm 3, Grp 2")
c[[1]] <- rlnorm(nsims, 2, .1)
c[[2]] <- rlnorm(nsims, 2, .1)
c[[3]] <- rlnorm(nsims, 11, .15)
c[[4]] <- rlnorm(nsims, 11, .15)
c[[5]] <- rlnorm(nsims, 11, .15)
c[[6]] <- rlnorm(nsims, 11, .15)

# effectiveness
e <- c
e[[1]] <- rnorm(nsims, 8, .2)
e[[2]] <- rnorm(nsims, 8, .2)
e[[3]] <- rnorm(nsims, 10, .8)
e[[4]] <- rnorm(nsims, 10.5, .8)
e[[5]] <- rnorm(nsims, 8.5, .6)
e[[6]] <- rnorm(nsims, 11, .6)

# cost and effectiveness by arm and simulation
ce <- data.table::data.table(sim = rep(seq(nsims), length(e)),
                 arm = rep(paste0("Arm ", seq(1, 3)), 
                           each = nsims * 2),
                 grp = rep(rep(c("Group 1", "Group 2"),
                               each = nsims), 3),
                 cost = do.call("c", c), qalys = do.call("c", e))

# net benefits by willingess to pay
krange <- seq(100000, 120000, 500)
kval1 <- sample(krange, 1)
kval2 <- sample(krange, 1)
ce[, nb1 := qalys * kval1 - cost]
ce[, nb2 := qalys * kval2 - cost]

# Functions to use for testing ------------------------------------------------
# most cost effective
mceR <- function(x, kval, grpname){
  x <- x[grp == grpname] 
  x[, nb := qalys * kval - cost]
  nb <- data.table::dcast(x, sim ~ arm, value.var = "nb")
  arms <- seq(1, ncol(nb) - 1)
  data.table::setnames(nb, colnames(nb), c("sim", paste0("nb", arms)))
  nb[, maxj := apply(nb[, -1, with = FALSE], 1, which.max)]
  nb[, maxj := factor(maxj, levels = arms)]
  return(as.numeric(prop.table(table(nb$maxj))))
}

# Test psa function -----------------------------------------------------------
## defaults
psa.dt <-  psa(ce, k = krange, sim = "sim", arm = "arm",
                grp = "grp", e = "qalys", c = "cost")

test_that("mce", {
  kval <- sample(krange, 1)
  
  # Group 1
  mce <- psa.dt$mce[grp == "Group 1" &  k == kval]
  mce.test <- mceR(ce, kval , "Group 1")
  expect_equal(mce$prob, mce.test)
  
  # Group 2
  mce <- psa.dt$mce[grp == "Group 2" &  k == kval]
  mce.test <- mceR(ce, kval , "Group 2")
  expect_equal(mce$prob, mce.test)
})