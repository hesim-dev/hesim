context("distributions.R unit tests")
library("flexsurv")

# Weibull distribution for NMA -------------------------------------------------
shape <-  1.2715
scale <- 6.1914
scalePH <- scale^{-shape}
a0 <- log(shape * scalePH)
a1 <- shape - 1

# pdf
expect_equal(stats::dweibull(.5, shape, scale),
             hesim::dweibullNMA(.5, a0 , a1))
expect_equal(flexsurv::dweibullPH(.5, shape, scalePH),
             hesim::dweibullNMA(.5, a0 , a1))

# cdf
expect_equal(stats::pweibull(4, shape, scale),
             hesim::pweibullNMA(4, a0 , a1))

# quantile
expect_equal(stats::qweibull(.8, shape, scale),
             hesim::qweibullNMA(.8, a0 , a1))

# random
set.seed(3)
r1 <- stats::rweibull(1, shape, scale)
set.seed(3)
r2 <- hesim::rweibullNMA(1, a0, a1)
expect_equal(r1, r2)

# hazard
expect_equal(flexsurv::hweibull(.8, shape, scale),
             hesim::hweibullNMA(.8, a0 , a1))
expect_equal(flexsurv::hweibull(.8, shape, scale, log = TRUE),
             hesim::hweibullNMA(.8, a0 , a1, log = TRUE))

# cumhazard 
expect_equal(flexsurv::Hweibull(.8, shape, scale),
             hesim::HweibullNMA(.8, a0 , a1))
expect_equal(flexsurv::Hweibull(.8, shape, scale, log = TRUE),
             hesim::HweibullNMA(.8, a0 , a1, log = TRUE))

# statistical modeling
## intercept only
s1 <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data = bc,
            dist = hesim.survdists$weibullNMA)
s2 <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data = bc,
            dist = "weibullPH")
s1.a1 <- s1$res.t["a1", "est"]; s1.a0 <- s1$res.t["a0", "est"]
s2.shape <- s2$res.t["shape", "est"]; s2.scale <- s2$res.t["scale", "est"]
expect_equal(exp(s1.a1), exp(s2.shape) - 1)
expect_equal(s1.a0, log(exp(s2.shape) * exp(s2.scale)))

# covariates on a0
s1 <- flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc,
            dist = hesim.survdists$weibullNMA)
s2 <- flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc,
            dist = "weibullPH")
s1.a0 <- s1$res.t["a0", "est"]; s1.a1 <- s1$res.t["a1", "est"]
s2.shape <- s2$res.t["shape", "est"]; s2.scale <- s2$res.t["scale", "est"]
expect_equal(exp(s1.a1), 
             exp(s2.shape) - 1)
expect_equal(log(exp(s2.scale) * exp(s2.shape)), 
             s1.a0)
expect_equal(s2.scale + s2.shape, # since shape is a scalar
             s1.a0)
expect_equal(s1$res.t["groupMedium"],
             s2$res.t["groupMedium"])

# covariates on a0 and a1
s1 <- suppressWarnings(flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc,
            anc = list(a1 = ~ group),
            dist = hesim.survdists$weibullNMA))
s2 <- suppressWarnings(flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc,
                  anc = list(shape = ~ group),
            dist = "weibullPH"))
s1.a0 <- s1$res.t["a0", "est"]; s1.a1 <- s1$res.t["a1", "est"]
s2.shape <- s2$res.t["shape", "est"]; s2.scale <- s2$res.t["scale", "est"]
expect_equal(exp(s1.a1),
             exp(s2.shape) - 1,
             scale = 1, tol = .001)
expect_equal(s2.scale + s2.shape, # since shape is a scalar
             s1.a0, 
             scale = 1, tol = .001)
expect_equal(exp(sum(s1$res.t[c("a1", "a1(groupMedium)"), "est"])),
             exp(sum(s2$res.t[c("shape", "shape(groupMedium)"), "est"])) - 1,
             scale = 1, tol = .001)


