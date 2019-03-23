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

## Create age and gender variables
ctstm3_exdata$transitions[, n := 1:.N, by = "id"]
ctstm3_exdata$transitions[, max_years := max(years), by = "id"]
pat_data <- ctstm3_exdata$transitions[n == 1]
pat_data[, age_mean := 75 * (max_years < 1) + 
                       65 * (max_years >= 1 & max_years < 2.5) +
                       50 * (max_years >= 2.5 & max_years < 5.5) +
                       40 * (max_years >= 5.5)]
pat_data[, female_prob := .70 * (max_years < 1) + 
                          .6 * (max_years >= 1 & max_years < 2.5) +
                          .40 * (max_years >= 2.5 & max_years < 5.5) +
                          .30 * (max_years >= 5.5)]
pat_data[, age := rnorm(nrow(pat_data), age_mean, 5)]
pat_data[, female := rbinom(nrow(pat_data), 1, female_prob)]
ctstm3_exdata$transitions <- merge(ctstm3_exdata$transitions,
                                   pat_data[, .(id, age, female)], 
                                   by = "id")
ctstm3_exdata$transitions[, age := ifelse(age >= 100, 100, age)]
ctstm3_exdata$transitions[, age := ifelse(age <= 18, 18, age)]
ctstm3_exdata$transitions[,  c("n", "N", "max_years") := NULL]
pat_data <- NULL

## Nice names for variables
setnames(ctstm3_exdata$transitions, "id", "patient_id")
setorderv(ctstm3_exdata$transitions,
          c("strategy_id", "patient_id", "age", "female",
            "from", "to", "trans",
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