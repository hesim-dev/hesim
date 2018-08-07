# Data for a 3-state continuous time state transition models
rm(list = ls())
ctstm3_exdata <- list()

# State transitions
library("mstate")
data(prothr)
ctstm3_exdata[["transitions"]] <- as.data.frame(prothr)
ctstm3_exdata$transitions$strategy_id <- ifelse(ctstm3_exdata$transitions$treat == "Placebo", 
                                                1, 2)
ctstm3_exdata$transitions$treat <- NULL

# Costs
## Drugs
ctstm3_exdata$costs$drug_costs <- data.frame(strategy_id = c(1, 2),
                                             costs = c(5000, 10000))

## Medical
medcosts <- data.frame(gender = rep(c("female", "male"), each = 2),
                       state_id = rep(c(1, 2), 2),
                       mean = c(800, 1200, 1000, 1600),
                       se = sqrt(c(100, 160, 120, 200)))
ctstm3_exdata$costs$medcosts <- medcosts

# Utility
ctstm3_exdata$utility <- data.frame(state_id = c(1, 2),
                                    mean = c(.65, .85),
                                    se = sqrt(c(.03, .04)))

# Save
save(ctstm3_exdata, file = "../data/ctstm3_exdata.rda", compress = "bzip2")