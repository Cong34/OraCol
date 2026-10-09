##############################################################################################
# OraCol.R
#
# A standalone R script that estimates how much oral bacteria is present in any
# gut-metagenomic sample set. One function, one taxa table in, one list of results out.
#
#   Output (a named list):
#     $species_table_oral      Relative abundance table (proportions, 0-1) of oral species only.
#                              Abundances are relative to the WHOLE community, so column sums
#                              = fraction of each gut sample that is oral.
#     $species_table_oral_CLR  CLR values of oral species. CLR is computed on the FULL community
#                              first and then subset (see note below).
#     $Summary_table_species   One row per oral species: mean / median / SD / sum of relative
#                              abundance, prevalence, mean and median CLR.
#     $Summary_table_samples   One row per sample: oral relative abundance sum, oral CLR sum and
#                              mean, oral richness + Shannon, whole-community richness + Shannon.

# Dependencies: TaxSEA (NCBI taxon ID lookup).
#
# Note on CLR: CLR is taken over all taxa in the input table and only then subset to oral
# taxa. If you CLR-transform the oral-only subset, each value becomes a log-ratio against the
# geometric mean of the *oral* taxa, which says nothing about how much oral bacteria sits in
# the gut. Taking it over the whole community means an oral taxon's CLR describes its
# abundance relative to the typical taxon in that gut sample.
##############################################################################################


###################### Main function #########################################################
#' @param species_table   data.frame or matrix of species-level or NCBI_ID abundances (counts,
#'                        proportions or percentages).
#' @param oral_definition Which oral database to use: "HOMD" (default), "eHOMD", "OES" or "mBodyMap".
#' @param taxa_are_rows   TRUE if species are rows and samples are columns (MetaPhlAn style).
#'
OraCol <- function(species_table = NULL,
                   oral_definition = c("HOMD", "eHOMD", "OES", "mBodyMap"),
                   taxa_are_rows = TRUE) {
  # First Define the Oral definition:
  oral_definition <- match.arg(oral_definition)
  Oral_taxa <- Oral_db[[oral_definition]]

  # Some lists contain the same NCBI ID more than once, so count unique IDs only.
  oral_taxa_id <- unique(trimws(as.character(Oral_taxa$ncbi_id)))
  oral_taxa_id <- oral_taxa_id[!is.na(oral_taxa_id) & oral_taxa_id != ""]
  message("Using oral definition: ", oral_definition, " (", length(oral_taxa_id), " unique taxa)")

  # Check basic conditions:
  if (is.null(species_table)) stop("`species_table` is required.")
  if (!is.matrix(species_table) && !is.data.frame(species_table))
    stop("`species_table` must be a matrix")

  # Turn the input into matrix:
  species_table <- as.matrix(species_table)
  if (!is.numeric(species_table))
    stop("`species_table` has non-numeric columns. Move species names into row names first.")

  # Put species in rows and samples in columns from here on.
  if (!taxa_are_rows) species_table <- t(species_table)

  if (is.null(rownames(species_table)) || is.null(colnames(species_table)))
    stop("`species_table` needs species names and sample names (row and column names). ",
         "If species are columns, set taxa_are_rows = FALSE.")
  if (anyNA(species_table))
    stop("`species_table` contains NA values.")
  if (any(species_table < 0))
    stop("`species_table` contains negative values.")

  # Remove 0 columns, Samples with no bacteria detected.
  n_samples_input <- ncol(species_table)
  empty <- colSums(species_table) == 0
  if (any(empty)) {
    warning(sum(empty), " sample(s) with zero total abundance removed: ",
            paste(colnames(species_table)[empty], collapse = ", "))
    species_table <- species_table[, !empty, drop = FALSE]
  }
  if (ncol(species_table) == 0) stop("All samples have zero total abundance.")

  # Relative abundance:
  rel_ab <- sweep(species_table, 2, colSums(species_table), "/")

  ############## Substitute the oral taxa into the species table ###############
  species_table_taxa  <- rownames(species_table)

  # Check to see if the rownames are NCBI ID
  is_id <- grepl("^[0-9]+$", species_table_taxa)
  # Sort the ID into one list (kept out of the search)
  species_ids <- ifelse(is_id, species_table_taxa, NA_character_)
  # Search the rest:
  if (any(!is_id)) {
    names_to_look_up <- species_table_taxa[!is_id]
    ncbi_lookup <- TaxSEA::get_ncbi_taxon_ids(names_to_look_up)
    species_ids[!is_id] <- vapply(ncbi_lookup[names_to_look_up],
                                  function(v) if (is.null(v)) NA_character_ else as.character(v[1]),
                                  character(1), USE.NAMES = FALSE)
  }

  is_oral <- !is.na(species_ids) & species_ids %in% oral_taxa_id

  if (!any(is_oral)) {
    stop("No taxa in `species_table` matched the oral definition. ")}

  # Creating a matching table for downstream rechecking:
  # NCBI_ID is NA when TaxSEA could not find the name.
  matching <- data.frame(input_name = species_table_taxa,
                         NCBI_ID    = species_ids,
                         is_oral    = is_oral,
                         stringsAsFactors = FALSE)

  # Build the Oral only table:
  species_table_oral <- rel_ab[is_oral, , drop = FALSE]

  ############## CLR transform the species table:  #############################
  pseudocount <- min(rel_ab[rel_ab > 0]) / 2

  log_ab <- log(rel_ab + pseudocount)
  clr_full <- sweep(log_ab, 2, colMeans(log_ab), "-")

  # Build the Oral centered log transformed table
  species_table_oral_CLR <- clr_full[is_oral, , drop = FALSE]


  # Making a summary table for species:
  #         + species       : species name
  #         + Mean_RA       : Average relative abundance across samples
  #         + Median_RA     : Middle value.
  #         + SD_RA         : How much abundance varies between samples
  #         + Sum_RA        : Relative abundance added up across all samples
  #         + Prevalence    : Fraction of samples where it was detected
  #         + N_samples_detected: Number of samples where it was detected
  #         + Mean_CLR      : Mean CLR value
  #         + Median_CLR    : Middle CLR value
  present <- species_table_oral > 0 # species detected (abundance above 0)

  Summary_table_species <- data.frame(
    Species        = rownames(species_table_oral),
    NCBI_ID        = species_ids[is_oral],
    Mean_RA        = rowMeans(species_table_oral),
    Median_RA      = apply(species_table_oral, 1, median),
    SD_RA          = apply(species_table_oral, 1, sd),
    Sum_RA         = rowSums(species_table_oral),
    Prevalence     = rowMeans(present),
    N_samples_detected = rowSums(present),
    Mean_CLR       = rowMeans(species_table_oral_CLR),
    Median_CLR     = apply(species_table_oral_CLR, 1, median),
    row.names = NULL, stringsAsFactors = FALSE
  )
  Summary_table_species <- Summary_table_species[order(-Summary_table_species$Prevalence), ]
  rownames(Summary_table_species) <- NULL

  # Making a summary table for samples:
  #         + Samples       : Samples name
  #         + Oral_RA_sum   : total relative abundance of Oral bacteria (proportion, 0-1)
  #         + Oral_CLR_sum  : Oral CLR value sums
  #         + Oral_CLR_mean : Average oral CLR value
  #         + Oral_richness : Number of oral species detected in the sample
  #         + Oral_shannon  : Shannon diversity of the oral species only
  #         + Total_richness: Number of species detected in the whole community
  #         + Total_shannon : Shannon diversity of the whole community

  # Ai generated Shannon diversity function: benchmark against VEGAN with great precision.
  shannon <- function(x) {
    x <- x[x > 0]                          # drop absent species
    if (length(x) == 0) return(NA_real_)   # nothing detected -> NA
    p <- x / sum(x)                        # rescale so they add to 1
    -sum(p * log(p))                       # Shannon formula
  }
  Oral_shannon = apply(species_table_oral, 2, shannon)

  Summary_table_samples <- data.frame(
    Sample            = colnames(species_table_oral),
    Oral_RA_sum       = colSums(species_table_oral),
    Oral_CLR_sum      = colSums(species_table_oral_CLR),
    Oral_CLR_mean     = colMeans(species_table_oral_CLR),
    Oral_richness     = colSums(species_table_oral > 0),
    Oral_shannon      = Oral_shannon,
    Total_richness    = colSums(rel_ab > 0),
    Total_shannon     = apply(rel_ab, 2, shannon),
    row.names = NULL, stringsAsFactors = FALSE
  )


  results_object = list("species_table_oral" = species_table_oral,
                        "species_table_oral_CLR" = species_table_oral_CLR,
                        "Summary_table_species" = Summary_table_species,
                        "Summary_table_samples" = Summary_table_samples,
                        "matching" = matching)
  return(results_object)
}

