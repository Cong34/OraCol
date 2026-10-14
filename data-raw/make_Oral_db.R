##############################################################################################
# make_Oral_db.R
#
# Builds data/Oral_db.rda: the oral databases that ship with the package. OraCol() loads this
# file with data(), and users can look at it with  data("Oral_db", package = "OraCol").
# Run this script whenever Oral_db_v4.rds is updated, with the package folder as your
# working directory (this happens automatically when you open the project's .Rproj file):
#
#   source("data-raw/make_Oral_db.R")
#
# Then reinstall the package (Ctrl+Shift+B) so it uses the new file.
##############################################################################################

Oral_db <- readRDS("Oral_db_v4.rds")

# Convert each list from a tibble to a plain data.frame so the package never needs tibble.
Oral_db <- lapply(Oral_db, as.data.frame, stringsAsFactors = FALSE)

save(Oral_db, file = "data/Oral_db.rda", compress = "xz")
