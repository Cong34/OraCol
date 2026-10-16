##############################################################################################
# tests/test_OraCol.R
#
# Checks that OraCol() still gives the right answers. Base R only, no extra packages.
#
# How to run (in RStudio, with the OraCol project open):
#   1. Ctrl+Shift+B                      install the current code
#   2. source("tests/test_OraCol.R")
#
# Each check prints PASS or FAIL. If any check fails, the file ends with an error, so
# devtools::check() and R CMD check (which run every file in tests/) report it too.
##############################################################################################


###################### Helpers ###############################################################
failed <- 0

# check("what should be true", test): prints PASS or FAIL. A test that crashes counts as a
# FAIL and its error is printed. `<<-` adds to the `failed` counter outside the function.
check <- function(description, test) {
  ok <- tryCatch(isTRUE(test), error = function(e) { cat("        error:", conditionMessage(e), "\n"); FALSE })
  cat(if (ok) "  PASS  " else "  FAIL  ", description, "\n", sep = "")
  if (!ok) failed <<- failed + 1
}
same     <- function(a, b) isTRUE(all.equal(a, b))                        # equal, allowing tiny rounding
run      <- function(...) suppressMessages(OraCol::OraCol(...))           # the package's OraCol, quietly
error_of <- function(code) tryCatch({ code; "" }, error = conditionMessage)   # its error message


###################### 1. Databases ##########################################################
cat("\n== 1. Databases ==\n")
data("Oral_db", package = "OraCol")
data("NCBI_lookup", package = "OraCol")

check("Oral_db has the four oral definitions",
      all(c("HOMD", "eHOMD", "OES", "mBodyMap") %in% names(Oral_db)))
check("every oral database has the columns OraCol() uses: ncbi_id, name, new_name, old_names",
      all(sapply(Oral_db, function(d) all(c("ncbi_id", "name", "new_name", "old_names") %in% names(d)))))
check("NCBI_lookup has the columns OraCol() uses: name, ncbi_id",
      all(c("name", "ncbi_id") %in% names(NCBI_lookup)))

if (failed > 0)
  stop("The databases in data/ are missing or incomplete, so the other tests were skipped.")


###################### 2. Calculations checked against hand-worked answers ###################
cat("\n== 2. Calculations ==\n")
# Toy table: 4 species x 3 samples. matrix() fills one column (= one sample) at a time.
# Species 1-2 are oral in every database; species 3-4 are not oral in HOMD (the default).
toy <- matrix(c(10,  30,  60,   0,     # S1, total 100
                25,  25,  25,  25,     # S2, total 100
                 0, 100,   0, 100),    # S3, total 200
              nrow = 4,
              dimnames = list(c("Streptococcus_mitis", "Rothia_mucilaginosa",
                                "Bacteroides_uniformis", "Faecalibacterium_prausnitzii"),
                              c("S1", "S2", "S3")))
res <- run(toy)
smp <- res$Summary_table_samples
spp <- res$Summary_table_species

# Relative abundance = each value / its sample's total:
#                                       S1     S2     S3
#   Streptococcus_mitis          (oral) 0.10   0.25   0
#   Rothia_mucilaginosa          (oral) 0.30   0.25   0.50
#   Bacteroides_uniformis               0.60   0.25   0
#   Faecalibacterium_prausnitzii        0      0.25   0.50
check("oral share per sample = 0.1 + 0.3, 0.25 + 0.25, 0 + 0.5",
      same(smp$Oral_RA_sum, c(0.4, 0.5, 0.5)))
check("oral richness = 2, 2, 1 and total richness = 3, 4, 2",
      same(smp$Oral_richness, c(2, 2, 1)) && same(smp$Total_richness, c(3, 4, 2)))

# Shannon = -sum(p * log(p)), with p rescaled to add up to 1.
#   Oral S1: 0.1 and 0.3 -> 0.25 and 0.75.  S2: two equal -> log(2).  S3: one species -> 0.
check("oral and total Shannon diversity",
      same(smp$Oral_shannon, c(-(0.25 * log(0.25) + 0.75 * log(0.75)), log(2), 0)) &&
        same(smp$Total_shannon, c(-(0.1 * log(0.1) + 0.3 * log(0.3) + 0.6 * log(0.6)), log(4), log(2))))

# Rothia: 0.30, 0.25, 0.50 (3 of 3 samples). S. mitis: 0.10, 0.25, 0 (2 of 3), so Rothia is first.
check("species summary: order, prevalence and samples detected",
      identical(spp$Species, c("Rothia_mucilaginosa", "Streptococcus_mitis")) &&
        same(spp$Prevalence, c(1, 2 / 3)) && same(spp$N_samples_detected, c(3, 2)))
check("species summary: mean, median, sum and SD of relative abundance",
      same(spp$Mean_RA, c(0.35, 0.35 / 3)) && same(spp$Median_RA, c(0.30, 0.10)) &&
        same(spp$Sum_RA, c(1.05, 0.35)) &&
        same(spp$SD_RA, c(sd(c(0.30, 0.25, 0.50)), sd(c(0.10, 0.25, 0)))))

# CLR: pseudocount = half the smallest non-zero relative abundance = 0.05.
# CLR = log(value + 0.05) minus the mean of those logs over ALL 4 species in the sample.
clr <- function(x) log(x + 0.05) - mean(log(x + 0.05))
expected_clr <- cbind(S1 = clr(c(0.10, 0.30, 0.60, 0))[1:2],
                      S2 = clr(c(0.25, 0.25, 0.25, 0.25))[1:2],    # all equal, so 0
                      S3 = clr(c(0, 0.50, 0, 0.50))[1:2])
rownames(expected_clr) <- c("Streptococcus_mitis", "Rothia_mucilaginosa")
check("CLR values (whole-community CLR, then the oral rows)",
      same(res$species_table_oral_CLR, expected_clr))
check("oral CLR sum and mean per sample, and mean CLR per species",
      same(smp$Oral_CLR_sum, unname(colSums(expected_clr))) &&
        same(smp$Oral_CLR_mean, unname(colMeans(expected_clr))) &&
        same(spp$Mean_CLR, unname(rowMeans(expected_clr))[2:1]))


###################### 3. Input formats and errors ###########################################
cat("\n== 3. Input formats and errors ==\n")
proportions <- sweep(toy, 2, colSums(toy), "/")
check("a data.frame gives the same results as a matrix", same(run(as.data.frame(toy)), res))
check("proportions and percentages give the same results as counts",
      same(run(proportions), res) && same(run(100 * proportions), res))
check("samples in rows with taxa_are_rows = FALSE gives the same results",
      same(run(t(toy), taxa_are_rows = FALSE), res))

toy_empty <- cbind(toy, S4 = 0)
check("an all-zero sample is removed, with a warning naming it",
      grepl("S4", tryCatch(run(toy_empty), warning = conditionMessage)) &&
        same(suppressWarnings(run(toy_empty)), res))

toy_na  <- toy; toy_na[1, 1]  <- NA
toy_neg <- toy; toy_neg[1, 1] <- -1
check("no table -> error",                       grepl("required", error_of(run(NULL))))
check("not a table -> error",                    grepl("must be a matrix", error_of(run(1:10))))
check("species names in a text column -> error",
      grepl("non-numeric", error_of(run(data.frame(Species = rownames(toy), toy)))))
check("no row or column names -> error",         grepl("names", error_of(run(unname(toy)))))
check("NA values -> error",                      grepl("NA values", error_of(run(toy_na))))
check("negative values -> error",                grepl("negative", error_of(run(toy_neg))))
check("every sample empty -> error",             grepl("All samples", error_of(suppressWarnings(run(toy * 0)))))
check("unknown oral_definition -> error",        grepl("should be one of", error_of(run(toy, oral_definition = "XYZ"))))
check("no oral species at all -> error",         grepl("No taxa", error_of(run(toy[3:4, ]))))


###################### 4. Detection threshold ################################################
cat("\n== 4. Detection threshold ==\n")
# S. mitis in S2 is 5 of 1,000,000 = 0.0005%, below the default threshold of 0.001%.
low <- matrix(c(100000, 300000, 600000,      0,     # S1, total 1,000,000
                     5, 499995, 250000, 250000),    # S2, total 1,000,000
              nrow = 4, dimnames = list(rownames(toy), c("S1", "S2")))
kept    <- run(low, detection_threshold = 0)$Summary_table_samples
dropped <- run(low)$Summary_table_samples

check("threshold 0 keeps the low value (S2: oral richness 2, oral share 0.5)",
      same(kept$Oral_richness, c(2, 2)) && same(kept$Oral_RA_sum, c(0.4, 0.5)))
# Default: the 5 becomes 0 and S2 is rescaled to the 999,995 that remain.
check("the default threshold drops it (S2: oral richness 1, oral share 499995/999995)",
      same(dropped$Oral_richness, c(2, 1)) && same(dropped$Oral_RA_sum, c(0.4, 499995 / 999995)))

filtered <- low
filtered["Streptococcus_mitis", "S2"] <- 0
check("a raw table and an already-filtered table give the same results", same(run(filtered), run(low)))
check("a negative threshold -> error", grepl("detection_threshold", error_of(run(toy, detection_threshold = -1))))
check("a threshold that removes everything -> error", grepl("nothing left", error_of(run(toy, detection_threshold = 99))))


###################### 5. Name matching ######################################################
cat("\n== 5. Name matching ==\n")
# NCBI IDs of the 4 toy species. 28037 and 43675 come from the oral database; 820 and 853 are
# not oral in HOMD, so they come from OraCol's NCBI name list (NCBI_lookup).
toy_ids <- c("28037", "43675", "820", "853")
ids_for <- function(names, definition = "HOMD") {
  m <- toy
  rownames(m) <- names
  run(m, oral_definition = definition)$matching
}
check("names from the oral database and from NCBI_lookup get their NCBI IDs",
      identical(res$matching$NCBI_ID, toy_ids))
check("spaces, s__ prefixes, full lineages and capital letters do not matter",
      identical(ids_for(c("Streptococcus mitis", "s__Rothia_mucilaginosa",
                          "k__Bacteria|g__Bacteroides|s__Bacteroides_uniformis",
                          "FAECALIBACTERIUM PRAUSNITZII"))$NCBI_ID, toy_ids))
toy_as_ids <- toy
rownames(toy_as_ids) <- toy_ids
check("NCBI IDs as row names give the same results as species names",
      same(run(toy_as_ids)$Summary_table_samples, smp))

# MetaPhlAn writes "Veillonella_sp_T11011_6"; the OES database writes "Veillonella sp. T11011-6".
check("punctuation does not matter (Veillonella_sp_T11011_6 is oral in OES)",
      ids_for(c("Veillonella_sp_T11011_6", rownames(toy)[2:4]), "OES")$is_oral[1])
# Actinomyces odontolyticus is the old name of Schaalia odontolytica (1660), which is in HOMD.
old_name <- ids_for(c("Actinomyces_odontolyticus", rownames(toy)[2:4]))
check("an old species name is matched (Actinomyces odontolyticus -> 1660, oral)",
      old_name$NCBI_ID[1] == "1660" && old_name$is_oral[1])
unknown <- ids_for(c(rownames(toy)[1:3], "GGB1234_SGB5678"))
check("a name in no list gets NCBI_ID NA and counts as not oral",
      is.na(unknown$NCBI_ID[4]) && !unknown$is_oral[4])


###################### 6. Output structure ###################################################
cat("\n== 6. Output structure ==\n")
check("the result has the 5 expected parts",
      identical(names(res), c("species_table_oral", "species_table_oral_CLR",
                              "Summary_table_species", "Summary_table_samples", "matching")))
check("Summary_table_samples has the expected columns, one row per sample",
      identical(names(smp), c("Sample", "Oral_RA_sum", "Oral_CLR_sum", "Oral_CLR_mean", "Oral_richness",
                              "Oral_shannon", "Total_richness", "Total_shannon")) &&
        identical(smp$Sample, colnames(toy)))
check("Summary_table_species has the expected columns",
      identical(names(spp), c("Species", "NCBI_ID", "Mean_RA", "Median_RA", "SD_RA", "Sum_RA",
                              "Prevalence", "N_samples_detected", "Mean_CLR", "Median_CLR")))
check("matching has one row per input species, in input order",
      identical(res$matching$input_name, rownames(toy)) &&
        identical(res$matching$is_oral, c(TRUE, TRUE, FALSE, FALSE)))


###################### 7. Example data (inst/extdata/sampling_data.csv) #####################
cat("\n== 7. Example data ==\n")
example <- read.csv(system.file("extdata", "sampling_data.csv", package = "OraCol"),
                    row.names = 1, check.names = FALSE)
res_ex <- run(example)
prev   <- res_ex$Summary_table_species$Prevalence
ra_sum <- res_ex$Summary_table_samples$Oral_RA_sum
# NCBI no longer recognises the name Ruminococcus bicirculans (not oral), so it has no ID.
check("every species gets an NCBI ID, except Ruminococcus bicirculans",
      all(res_ex$matching$input_name[is.na(res_ex$matching$NCBI_ID)] %in% "Ruminococcus_bicirculans"))
check("oral share is between 0 and 1, and prevalence is sorted from high to low",
      all(ra_sum >= 0 & ra_sum <= 1) && all(diff(prev) <= 0))
for (definition in c("HOMD", "eHOMD", "OES", "mBodyMap"))
  check(paste0("runs with oral_definition = \"", definition, "\" and finds oral species"),
        any(run(example, oral_definition = definition)$matching$is_oral))


###################### Summary ###############################################################
if (failed > 0) stop(failed, " check(s) failed - see the FAIL lines above.")
cat("\nAll checks passed.\n")
