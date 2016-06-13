# sim_los ----------------------------------------------------------------------
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
loc_beta <-  array(NA, c(1, 2, nrow(tmat)))
for (i in 1:nrow(tmat)) loc_beta[, , i] <- coef(mod[[i]])[c("scale", "treat")]
loc_x <- as.matrix(data.frame(int = 1, treat = 1))
loc_x <- loc_x[rep(seq_len(nrow(loc_x)), each = N), ]
par2 <- unlist(lapply(mod, function (x) coef(x)["shape"]))
sim.cea <- sim_msm(loc_beta, loc_x, dist, tmat, par2, maxt)

# compare simulations
mean(sim.flexsurv$t[, ncol(sim.flexsurv$t)])
sim_los(sim.cea)

# sim_msm_pv -------------------------------------------------------------------
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
