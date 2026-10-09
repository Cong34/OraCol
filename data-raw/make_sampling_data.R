##############################################################################################
# make_sampling_data.R
#
# Builds inst/extdata/sampling_data.csv: a fake species-level count table (species in rows,
# 40 samples in columns) used to show the input format and to test OraCol().
# Run from the package root (open the project's .Rproj file first):
#
#   source("data-raw/make_sampling_data.R")
#
# Species names are real "Genus_species" names taken from TaxSEA's offline NCBI lookup, so
# OraCol() can match every name without querying NCBI online. Counts are random.
# Once installed, users can find the file with:
#   system.file("extdata", "sampling_data.csv", package = "OraCol")
##############################################################################################

set.seed(2026)   # same seed -> same file every time

# Pool of species names TaxSEA knows offline, keeping strict "Genus_species" names only
# (drops entries like "Eubacterium_sp.", "Blautia_" or "Firmicutes_bacterium").
data("NCBI_ids", package = "TaxSEA")
NCBI_ids <- NCBI_ids[lengths(NCBI_ids) > 0]   # some names are stored with no ID; drop them
species_pool <- grep("^[A-Z][a-z]+_[a-z]+$", names(NCBI_ids), value = TRUE)
species_pool <- setdiff(species_pool, "Firmicutes_bacterium")

# Split the pool into oral (in any of the three oral lists) and non-oral species.
load("R/sysdata.rda")   # loads Oral_db_v2
oral_names <- unique(unlist(lapply(Oral_db_v2, function(d) d$Taxa)))
oral_pool     <- intersect(species_pool, oral_names)
non_oral_pool <- setdiff(species_pool, oral_names)

# Pick 300-500 species in total, about a quarter of them oral.
n_species <- sample(300:500, 1)
n_oral    <- round(n_species / 4)
species   <- c(sample(oral_pool, n_oral), sample(non_oral_pool, n_species - n_oral))
species   <- sample(species)   # shuffle so oral species are not all at the top

# Random counts 0-1000, then set ~60% to 0 (real tables have many absent species).
n_samples <- 40
counts <- matrix(sample(0:1000, n_species * n_samples, replace = TRUE),
                 nrow = n_species,
                 dimnames = list(species, sprintf("S%02d", 1:n_samples)))
counts[runif(length(counts)) < 0.6] <- 0

# Write with species names as the first column, called "Species".
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(Species = rownames(counts), counts, check.names = FALSE),
          "inst/extdata/sampling_data.csv", row.names = FALSE)

message("Wrote ", n_species, " species (", n_oral, " oral) x ", n_samples, " samples.")
