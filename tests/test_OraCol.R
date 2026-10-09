##############################################################################################
# tests/test_OraCol.R
#
# Automatic checks that OraCol() still gives the right answers. Base R only - no packages.
#
# How to run (in RStudio, with the OraCol project open):
#   1. Ctrl+Shift+B                      installs the current code and restarts R cleanly
#   2. source("tests/test_OraCol.R")
#
# Every check prints PASS or FAIL. If anything fails, the file stops with an error at the end,
# so devtools::check() / R CMD check (which run every .R file in tests/ automatically) also
# report the failure.
#
# All species names used here are in TaxSEA's offline lookup, so no internet is needed.
##############################################################################################


###################### Helpers ###############################################################
n_checks <- 0
n_failed <- 0

# check("what should be true", expression): prints PASS if the expression is TRUE, otherwise
# FAIL (plus the error message if the expression itself crashed). `<<-` updates the two
# counters above, which live outside the function.
check <- function(description, expr) {
  result <- tryCatch(isTRUE(expr), error = function(e) conditionMessage(e))
  n_checks <<- n_checks + 1
  if (isTRUE(result)) {
    cat("  PASS  ", description, "\n", sep = "")
  } else {
    n_failed <<- n_failed + 1
    cat("  FAIL  ", description, "\n", sep = "")
    if (is.character(result)) cat("        error: ", result, "\n", sep = "")
  }
}

# TRUE if equal up to tiny rounding errors (computers store decimals approximately).
same <- function(a, b) isTRUE(all.equal(a, b))

# Always run the PACKAGE version of OraCol (OraCol::), never a copy loaded with source(),
# and hide its "Using oral definition" message to keep the output readable.
run <- function(...) suppressMessages(OraCol::OraCol(...))

# The text of the error / warning an expression produces ("no error" / "no warning" if none).
error_message   <- function(expr) tryCatch({ expr; "no error" },   error   = conditionMessage)
warning_message <- function(expr) tryCatch({ expr; "no warning" }, warning = conditionMessage)

cat("Testing OraCol", format(packageVersion("OraCol")), "from", find.package("OraCol"), "\n")


###################### 1. Internal oral database (R/sysdata.rda) ############################
# OraCol() reads Oral_db from inside the package. These checks catch a database that was not
# rebuilt after a change, e.g. a new oral definition added to OraCol() but not to sysdata.rda.
# If you rename the database object, change "Oral_db" here too.
cat("\n== 1. Internal oral database ==\n")

ns <- asNamespace("OraCol")   # the package's private storage, where sysdata.rda is loaded
# inherits = FALSE: look ONLY inside the package, ignoring any Oral_db in your own workspace.
db_found <- exists("Oral_db", envir = ns, inherits = FALSE)
check("R/sysdata.rda contains an object called Oral_db", db_found)

db      <- if (db_found) get("Oral_db", envir = ns) else list()
choices <- eval(formals(OraCol::OraCol)$oral_definition)   # the choices OraCol() offers
missing_lists <- setdiff(choices, names(db))
check(paste0("every oral_definition choice (", paste(choices, collapse = ", "), ") is in Oral_db",
             if (length(missing_lists)) paste0(" - missing: ", paste(missing_lists, collapse = ", "))),
      length(missing_lists) == 0)
check("every list in Oral_db has rows and an ncbi_id column",
      length(db) > 0 &&
        all(vapply(db, function(d) nrow(d) > 0 && "ncbi_id" %in% names(d), logical(1))))

if (n_failed > 0)
  stop("The internal oral database is missing or out of date, so the remaining tests were skipped.\n",
       "  Rebuild it with source(\"data-raw/make_Oral_db.R\"), reinstall (Ctrl+Shift+B), then rerun.")


###################### 2. Calculations checked against hand-worked answers ###################
cat("\n== 2. Calculations vs hand-worked answers ==\n")

# Toy table: 4 species x 3 samples. matrix() fills one column (= one sample) at a time.
# Species 1-2 are oral in every oral list; species 3-4 are not oral in HOMD or eHOMD.
toy <- matrix(c(10,  30,  60,   0,     # S1, total 100
                25,  25,  25,  25,     # S2, total 100
                 0, 100,   0, 100),    # S3, total 200
              nrow = 4,
              dimnames = list(c("Streptococcus_mitis", "Rothia_mucilaginosa",
                                "Bacteroides_uniformis", "Faecalibacterium_prausnitzii"),
                              c("S1", "S2", "S3")))
res <- run(toy)                     # default oral definition: HOMD
smp <- res$Summary_table_samples
spp <- res$Summary_table_species

# Relative abundance = each value / its sample's total:
#                                       S1     S2     S3
#   Streptococcus_mitis          (oral) 0.10   0.25   0
#   Rothia_mucilaginosa          (oral) 0.30   0.25   0.50
#   Bacteroides_uniformis               0.60   0.25   0
#   Faecalibacterium_prausnitzii        0      0.25   0.50
check("species_table_oral holds the oral rows of the relative abundance table",
      same(res$species_table_oral,
           rbind(Streptococcus_mitis = c(S1 = 0.10, S2 = 0.25, S3 = 0),
                 Rothia_mucilaginosa = c(S1 = 0.30, S2 = 0.25, S3 = 0.50))))

# Oral share per sample: S1 = 0.1 + 0.3, S2 = 0.25 + 0.25, S3 = 0 + 0.5
check("Oral_RA_sum = 0.4, 0.5, 0.5", same(smp$Oral_RA_sum, c(0.4, 0.5, 0.5)))
check("Oral_richness = 2, 2, 1 (oral species above 0)", same(smp$Oral_richness, c(2, 2, 1)))
check("Total_richness = 3, 4, 2 (all species above 0)", same(smp$Total_richness, c(3, 4, 2)))

# Shannon H = -sum(p * log(p)), where p is rescaled to add up to 1 over the species used.
#   Oral  S1: 0.1, 0.3 -> p = 0.25, 0.75   S2: two equal -> log(2)   S3: one species -> 0
#   Total S1: 0.1, 0.3, 0.6                S2: four equal -> log(4)  S3: two equal -> log(2)
check("Oral_shannon matches hand calculation",
      same(smp$Oral_shannon, c(-(0.25 * log(0.25) + 0.75 * log(0.75)), log(2), 0)))
check("Total_shannon matches hand calculation",
      same(smp$Total_shannon,
           c(-(0.1 * log(0.1) + 0.3 * log(0.3) + 0.6 * log(0.6)), log(4), log(2))))

# Species summary. Rothia is in 3 of 3 samples, S. mitis in 2 of 3, so Rothia comes first.
#   Rothia:   0.30, 0.25, 0.50 -> mean 1.05 / 3 = 0.35, median 0.30, sum 1.05
#   S. mitis: 0.10, 0.25, 0    -> mean 0.35 / 3,        median 0.10, sum 0.35
check("species summary is sorted by prevalence (Rothia first)",
      identical(spp$Species, c("Rothia_mucilaginosa", "Streptococcus_mitis")))
check("Prevalence = 1 and 2/3", same(spp$Prevalence, c(1, 2 / 3)))
check("N_samples_detected = 3 and 2", same(spp$N_samples_detected, c(3, 2)))
check("Mean_RA = 0.35 and 0.35/3", same(spp$Mean_RA, c(0.35, 0.35 / 3)))
check("Median_RA = 0.30 and 0.10", same(spp$Median_RA, c(0.30, 0.10)))
check("Sum_RA = 1.05 and 0.35", same(spp$Sum_RA, c(1.05, 0.35)))
check("SD_RA is the SD of each species' relative abundances",
      same(spp$SD_RA, c(sd(c(0.30, 0.25, 0.50)), sd(c(0.10, 0.25, 0)))))
check("NCBI_ID comes from TaxSEA (43675, 28037)", identical(spp$NCBI_ID, c("43675", "28037")))

# CLR. Pseudocount = half the smallest non-zero relative abundance in the table: 0.1 / 2 = 0.05.
# CLR = log(value + 0.05) minus the mean of those logs over ALL 4 species in the sample
# (the whole community, not only the oral species - see the note at the top of R/OraCol.R).
clr <- function(x) log(x + 0.05) - mean(log(x + 0.05))
clr_S1 <- clr(c(0.10, 0.30, 0.60, 0))
clr_S2 <- clr(c(0.25, 0.25, 0.25, 0.25))   # all equal, so every CLR value is exactly 0
clr_S3 <- clr(c(0, 0.50, 0, 0.50))
expected_clr <- rbind(Streptococcus_mitis = c(S1 = clr_S1[1], S2 = clr_S2[1], S3 = clr_S3[1]),
                      Rothia_mucilaginosa = c(S1 = clr_S1[2], S2 = clr_S2[2], S3 = clr_S3[2]))
check("CLR values match hand calculation (whole-community CLR)",
      same(res$species_table_oral_CLR, expected_clr))
check("Oral_CLR_sum and Oral_CLR_mean add up each sample's oral CLR values",
      same(smp$Oral_CLR_sum, unname(colSums(expected_clr))) &&
        same(smp$Oral_CLR_mean, unname(colMeans(expected_clr))))
check("Mean_CLR and Median_CLR per species",
      same(spp$Mean_CLR,   c(mean(expected_clr[2, ]),   mean(expected_clr[1, ]))) &&
        same(spp$Median_CLR, c(median(expected_clr[2, ]), median(expected_clr[1, ]))))


###################### 3. Input handling #####################################################
cat("\n== 3. Input handling ==\n")

# The same data in a different format must give exactly the same results.
proportions <- sweep(toy, 2, colSums(toy), "/")
check("a data.frame gives the same results as a matrix", same(run(as.data.frame(toy)), res))
check("proportions (0-1) give the same results as counts", same(run(proportions), res))
check("percentages (0-100) give the same results as counts", same(run(100 * proportions), res))
check("samples as rows with taxa_are_rows = FALSE gives the same results",
      same(run(t(toy), taxa_are_rows = FALSE), res))

# A sample where nothing was detected is removed with a warning; the rest is unchanged.
toy_empty <- cbind(toy, S4 = 0)
check("an all-zero sample gives a warning naming it", grepl("S4", warning_message(run(toy_empty))))
check("...and is removed, leaving the other results unchanged",
      same(suppressWarnings(run(toy_empty)), res))

# Bad input must stop with a clear error message.
toy_na  <- toy; toy_na[1, 1]  <- NA
toy_neg <- toy; toy_neg[1, 1] <- -1
check("no table -> 'required' error", grepl("required", error_message(run(NULL))))
check("not a table -> 'must be a matrix' error", grepl("must be a matrix", error_message(run(1:10))))
check("species names in a text column -> 'non-numeric' error",
      grepl("non-numeric", error_message(run(data.frame(Species = rownames(toy), toy)))))
check("no row/column names -> 'names' error", grepl("names", error_message(run(unname(toy)))))
check("NA values -> 'NA values' error", grepl("NA values", error_message(run(toy_na))))
check("negative values -> 'negative' error", grepl("negative", error_message(run(toy_neg))))
check("every sample empty -> 'All samples' error",
      grepl("All samples", error_message(suppressWarnings(run(toy * 0)))))
check("unknown oral_definition -> 'should be one of' error",
      grepl("should be one of", error_message(run(toy, oral_definition = "XYZ"))))
check("no oral species at all -> 'No taxa' error", grepl("No taxa", error_message(run(toy[3:4, ]))))


###################### 4. Output structure ###################################################
cat("\n== 4. Output structure ==\n")
check("result has the 5 expected parts",
      identical(names(res), c("species_table_oral", "species_table_oral_CLR",
                              "Summary_table_species", "Summary_table_samples", "matching")))
check("Summary_table_samples has the expected columns",
      identical(names(smp), c("Sample", "Oral_RA_sum", "Oral_CLR_sum", "Oral_CLR_mean",
                              "Oral_richness", "Oral_shannon", "Total_richness", "Total_shannon")))
check("Summary_table_species has the expected columns",
      identical(names(spp), c("Species", "NCBI_ID", "Mean_RA", "Median_RA", "SD_RA", "Sum_RA",
                              "Prevalence", "N_samples_detected", "Mean_CLR", "Median_CLR")))
check("one Summary_table_samples row per sample", identical(smp$Sample, c("S1", "S2", "S3")))
check("matching has one row per input species, in input order",
      identical(res$matching$input_name, rownames(toy)))
check("matching marks exactly the 2 oral species",
      identical(res$matching$is_oral, c(TRUE, TRUE, FALSE, FALSE)))
check("matching records every NCBI ID",
      identical(res$matching$NCBI_ID, c("28037", "43675", "820", "853")))


###################### 5. Example data (inst/extdata/sampling_data.csv) #####################
cat("\n== 5. Example data ==\n")
# as.matrix(): data.frame input is already tested in section 3, so a problem there does not
# hide the results of these checks.
ex <- as.matrix(read.csv(system.file("extdata", "sampling_data.csv", package = "OraCol"),
                         row.names = 1, check.names = FALSE))
res_ex <- run(ex)
prev   <- res_ex$Summary_table_species$Prevalence
ra_sum <- res_ex$Summary_table_samples$Oral_RA_sum

check("every species gets an NCBI ID (no online lookups needed)", !anyNA(res_ex$matching$NCBI_ID))
check("one matching row per species and one summary row per sample",
      nrow(res_ex$matching) == nrow(ex) && nrow(res_ex$Summary_table_samples) == ncol(ex))
check("one species summary row per oral species",
      nrow(res_ex$Summary_table_species) == sum(res_ex$matching$is_oral))
check("Oral_RA_sum is between 0 and 1", all(ra_sum >= 0 & ra_sum <= 1))
# diff() = difference between neighbours; all <= 0 means each value is no larger than the last.
check("Prevalence is between 0 and 1 and sorted high to low",
      all(prev >= 0 & prev <= 1) && all(diff(prev) <= 0))
for (definition in choices)
  check(paste0("runs with oral_definition = \"", definition, "\" and finds oral species"),
        any(run(ex, oral_definition = definition)$matching$is_oral))


###################### Summary ###############################################################
cat("\n", n_checks - n_failed, " of ", n_checks, " checks passed.\n", sep = "")
if (n_failed > 0) stop(n_failed, " check(s) failed - see the FAIL lines above.")
