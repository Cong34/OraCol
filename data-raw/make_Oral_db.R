##############################################################################################
# make_Oral_db.R
#
# Builds R/sysdata.rda: the oral taxa lists that ship INSIDE the package.
# Run this script whenever Oral_db_v2.rds is updated, with the package folder as your
# working directory (this happens automatically when you open the project's .Rproj file):
#
#   source("data-raw/make_Oral_db.R")
#
# R only auto-loads a file with the exact name R/sysdata.rda. Objects saved in it (here
# Oral_db_v2) are available to the package's own functions by name, but hidden from users.
##############################################################################################

Oral_db <- readRDS("Oral_db_v3.rds")

# Convert each list from a tibble to a plain data.frame so the package never needs tibble.
Oral_db <- lapply(Oral_db, as.data.frame, stringsAsFactors = FALSE)

save(Oral_db, file = "R/sysdata.rda", compress = "xz")
