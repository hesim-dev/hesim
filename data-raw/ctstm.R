# Data for a 3-state continuous time state transition models
rm(list = ls())
ctstm3_exdata <- list()

# State transitions
library("mstate")
library("data.table")
data(prothr)
ctstm3_exdata[["transitions"]] <- as.data.table(prothr)
ctstm3_exdata$transitions[, Tstart := Tstart/365.25]
ctstm3_exdata$transitions[, Tstop := (Tstop + 1)/365.25]
ctstm3_exdata$transitions[, years := Tstop - Tstart]
ctstm3_exdata$transitions[, strategy_id := ifelse(treat == "Placebo", 1, 2)]
ctstm3_exdata$transitions[, treat := NULL]
setnames(ctstm3_exdata$transitions, "id", "patient_id")
setorderv(ctstm3_exdata$transitions,
          c("strategy_id", "patient_id", "from", "to", "trans",
            "Tstart", "Tstop", "years", "status"))
ctstm3_exdata$transitions <- data.frame(ctstm3_exdata$transitions)

# Costs
## Drugs
ctstm3_exdata$costs$drugs <- data.frame(strategy_id = c(1, 2),
                                             costs = c(5000, 10000))

## Medical
medcosts <- data.frame(state_id = c(1, 2),
                       mean = c(1000, 1600),
                       se = sqrt(c(120, 200)))
ctstm3_exdata$costs$medical <- medcosts

# Utility
ctstm3_exdata$utility <- data.frame(state_id = c(1, 2),
                                    mean = c(.65, .85),
                                    se = sqrt(c(.03, .04)))

# Save
save(ctstm3_exdata, file = "../data/ctstm3_exdata.rda", compress = "bzip2")