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
pv.qol <- sim_pv(sim.cea, tmat, beta = c(1, 1, 0), r = 0, name = "qalys")
sim.cea <- cbind(sim.cea, pv.qol)

# compare simulations
sim.cea.qalys <- sim_los(sim.cea, timevar = "qalys")
sum(sim.cea.qalys$los)
mean(sim.flexsurv$t[, ncol(sim.flexsurv$t)])

## Simulation with Age
dat$age <- rnorm(nrow(dat), 65, 10)
mod2 <- vector(3, mode = "list")
for (i in 1:3){
  mod2[[i]] <- flexsurvreg(Surv(years, status) ~ treat + age,
                          subset = (trans == i),  data = dat, dist = dist[i])
}
param2 <- sim_param(mod2)
loc_x2 <- cbind(loc_x, rnorm(nrow(loc_x), 50, 5))
colnames(loc_x2)[3] <- "age"
sim.cea2 <- sim_msm(param2$loc_beta, loc_x2, dist, tmat, param2$anc1, maxt,
                    agevar = "age")
pv.qol2 <- sim_pv(sim.cea2, tmat, beta = c(1, 1, 0), r = 0, name = "qalys")
sim.cea2 <- cbind(sim.cea2, pv.qol2)
sim_los(sim.cea2, timevar = "qalys")

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

# calculations
pv.cost <- sim_pv(sim.cea, tmat, x = x.cost, beta = beta, poly.beta = poly.beta,
       poly.deg = poly.deg, knots = knots, r = 0, name = "costs")
sim.cea <- cbind(sim.cea, pv.cost)

## Simulation with Age
x.cost2 <- cbind(x.cost, rnorm(nrow(x.cost), 65, 10))
colnames(x.cost2)[3] <- "age"
beta2 <- cbind(beta, rep(1, 3))
pv.cost2 <- sim_pv(sim.cea2, tmat, x = x.cost2, beta = beta2, poly.beta = poly.beta,
                  poly.deg = poly.deg, knots = knots, r = 0, name = "costs")
sim.cea2 <- cbind(sim.cea2, pv.cost2)
