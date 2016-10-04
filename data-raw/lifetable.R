lifetable <- read.csv("data-raw/lifetable.csv")
lifetable_male <- lifetable[lifetable$male == 1, ]
lifetable_female <- lifetable[lifetable$male == 0, ]
lifetable_male$male <- NULL
lifetable_female$male <- NULL
save(lifetable_male, file = "data/lifetable-male.rda", compress = "bzip2")
save(lifetable_female, file = "data/lifetable-female.rda", compress = "bzip2")
