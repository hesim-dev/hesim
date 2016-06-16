# sim_msm and sim_los ----------------------------------------------------------
## Simulation without Age
# estimate multi-state model
library("flexsurv")
library("data.table")
dat <- bosms3
dat$treat <- rbinom(nrow(dat), 1, .5)
mod <- vector(3, mode = "list")
dist <- rep("weibull", 3)
for (i in 1:3){
  mod[[i]] <- flexsurvreg(Surv(years, status) ~ treat,
                             subset = (trans == i),  data = dat, dist = dist[i])
}

# flexsurv simulation
tmat <- rbind(c(NA, 1, 2), c(NA, NA, 3), c(NA, NA, NA))
maxt <- 30
N <- 1e5
sim.flexsurv <- sim.fmsm(x = mod, trans = tmat, newdata = data.frame(treat = 1),
                         t = maxt, M = N)

# cea simulation
param <- sim_param(mod)
loc_x <- as.matrix(data.frame(int = 1, treat = 1))
loc_x <- loc_x[rep(seq_len(nrow(loc_x)), each = N), ]
sim.cea <- sim_msm(param$loc_beta, loc_x, dist, tmat, param$anc1, maxt)

# compare simulations
mean(sim.flexsurv$t[, ncol(sim.flexsurv$t)])
sim_los(sim.cea)

## Simulation with Age
dat$age <- rnorm(nrow(dat), 65, 10)
mod2 <- vector(3, mode = "list")
for (i in 1:3){
  mod2[[i]] <- flexsurvreg(Surv(years, status) ~ treat + age,
                          subset = (trans == i),  data = dat, dist = dist[i])
}
loc_beta2 <-  array(NA, c(1, 3, nrow(tmat)))
for (i in 1:nrow(tmat)) loc_beta2[, , i] <- coef(mod2[[i]])[c("scale", "treat", "age")]
loc_x2 <- as.matrix(data.frame(int = 1, treat = 1, age = 50))
loc_x2 <- loc_x2[rep(seq_len(nrow(loc_x2)), each = N), ]
par2.2 <- unlist(lapply(mod2, function (x) coef(x)["shape"]))
sim.cea2 <- sim_msm(loc_beta2, loc_x2, dist, tmat, par2.2, maxt, agevar = "age")
sim_los(sim.cea2)

# sim_msm_pv -------------------------------------------------------------------
## Simulation without Age
# cost data
x.cost <- loc_x
beta <- c(500, 500)
beta <- t(replicate(3, beta))
poly.beta <- c(1000, 0,  -.5, 1, -1.5, -1000)
poly.beta <- t(replicate(3, poly.beta))
poly.deg <- c(0, 3, 0)
poly.deg <- t(replicate(3, poly.deg))
knots <- c(1, 4)
knots <- t(replicate(3, knots))

# utility
x.util <- loc_x[, "int", drop = FALSE]
qol <- rep(1, 3)

# calculations
pv.cost <- sim_pv(sim.cea, x = x.cost, beta = beta, poly.beta = poly.beta,
       poly.deg = poly.deg, knots = knots)
pv.qol <- sim_pv(sim.cea, x = x.util, beta = qol)

## Simulation with Age
x.cost2 <- cbind(x.cost, rnorm(nrow(x.cost), 65, 10))
colnames(x.cost2)[3] <- "age"
beta2 <- cbind(beta, rep(1, 3))
pv.cost2 <- sim_pv(sim.cea2, x = x.cost2, beta = beta2, poly.beta = poly.beta,
                  poly.deg = poly.deg, knots = knots)
