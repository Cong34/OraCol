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
#                              Samples with no oral species are left out by default.
#     $matching                Every input species, its NCBI ID and whether it counted as oral.
#     $no_oral_samples         Samples with no oral species detected (same columns as above).
# Dependencies: none beyond R itself (stats and utils ship with R).
##############################################################################################
###################### Main function #########################################################
#' Summarise oral bacteria in a metagenomic species table
#'
#' Finds which species in a species table are oral bacteria, using one of four oral
#' reference databases, and returns per-species and per-sample summaries: relative
#' abundance, prevalence, centred log-ratio (CLR) values, richness and Shannon diversity.
#'
#' @param species_table A matrix or data.frame of species-level abundances (counts,
#'   proportions or percentages). Row names are species names or NCBI taxonomy IDs;
#'   column names are sample names.
#' @param oral_definition The oral database to use: `"HOMD"` (default), `"eHOMD"`,
#'   `"OES"` or `"mBodyMap"`. See the section **Oral definitions**.
#' @param taxa_are_rows `TRUE` (default) if species are rows and samples are columns.
#'   Use `FALSE` if it is the other way round.
#' @param detection_threshold Lowest relative abundance at which a species is considered
#'   detected in the sample. The default `0.001` means 0.001%.
#' @param remove_no_oral `TRUE` (default) removes samples with no oral species detected
#'   from `Summary_table_samples`. They are listed in `no_oral_samples` instead. Use
#'   `FALSE` to keep them in `Summary_table_samples`.
#'
#' @details
#' OraCol first take a Species table and convert it to relative abundance to apply its threshold.
#' After removing all the species that are not detected, it will re-scale the table to add up to 1.
#'
#' The Detection are divided into 3 steps asking basic questions:
#'
#' 1. Is the species detected in NCBI ID format
#' 2. Does its name matched with the Oral databases (old or new)
#' 3. Is it in OraCol's NCBI name list ([NCBI_lookup]), which was built once from NCBI so
#'    no internet connection is needed
#'
#' With these steps, the aim is to convert the bacteria species name into ncbi id for easier comparison.
#' A species counts as oral when its ID is in the chosen oral database.
#'
#' CLR values are calculated over the whole community and then subset to the oral
#' species. Zeros are handled by adding half the smallest non-zero relative abundance.
#' In a sample with no oral species detected, every oral species gets the same small
#' negative CLR value. That value comes only from the non-oral species, so it carries
#' no oral signal. For this reason such samples are left out of `Summary_table_samples`
#' by default (see `remove_no_oral`). They stay in `species_table_oral`,
#' `species_table_oral_CLR` and the prevalence values in `Summary_table_species`.
#'
#' @section Oral definitions:
#' * `"HOMD"`: Human Oral Microbiome Database (Chen et al. 2010).
#' * `"eHOMD"`: expanded HOMD, which adds nasal and upper-airway taxa (Escapa et al. 2018).
#' * `"OES"`: Oral Enrichment Score at 1% (Manghi 2025).
#' * `"mBodyMap"`: mBodyMap, oral sites of healthy people only (Jin et al. 2022).
#'
#' @return A list with six elements:
#' * `species_table_oral`: relative abundance (0 to 1) of the oral species in each sample.
#' * `species_table_oral_CLR`: CLR values of the oral species.
#' * `Summary_table_species`: one row per oral species, sorted by prevalence.
#' * `Summary_table_samples`: one row per sample, with the oral relative abundance,
#'   richness and Shannon diversity. By default, samples with no oral species detected
#'   are left out.
#' * `matching`: every input species, its NCBI ID (`NA` if none was found) and whether it
#'   counted as oral.
#' * `no_oral_samples`: the samples with no oral species detected, with the same columns
#'   as `Summary_table_samples`. Empty if every sample has oral species.
#'
#' @examples
#' # Example table shipped with the package: 420 species x 40 samples
#' f <- system.file("extdata", "sampling_data.csv", package = "OraCol")
#' species_table <- read.csv(f, row.names = 1, check.names = FALSE)
#'
#' res <- OraCol(species_table, oral_definition = "HOMD")
#' head(res$Summary_table_samples)
#' head(res$Summary_table_species)
#'
#' @importFrom stats median sd
#' @export
OraCol <- function(species_table = NULL,
                   oral_definition = c("HOMD", "eHOMD", "OES", "mBodyMap"),
                   taxa_are_rows = TRUE,
                   detection_threshold = 0.001,
                   remove_no_oral = TRUE) {
  # First Define the Oral definition:
  oral_definition <- match.arg(oral_definition)
  # Load the oral databases and the NCBI name list that ship in the package's data/ folder.
  utils::data("Oral_db", package = "OraCol", envir = environment())
  utils::data("NCBI_lookup", package = "OraCol", envir = environment())
  Oral_taxa <- Oral_db[[oral_definition]]

  # Some lists contain the same NCBI ID more than once, so count unique IDs only.
  oral_taxa_id <- unique(trimws(as.character(Oral_taxa$ncbi_id)))
  oral_taxa_id <- oral_taxa_id[!is.na(oral_taxa_id) & oral_taxa_id != ""]
  message("Using oral definition: ", oral_definition)

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
    stop("`species_table` needs species names and sample names (row and column names). ")
  if (anyNA(species_table))
    stop("`species_table` contains NA values.")
  if (any(species_table < 0))
    stop("`species_table` contains negative values.")
  if (!is.numeric(detection_threshold) || length(detection_threshold) != 1 ||
      is.na(detection_threshold) || detection_threshold < 0 || detection_threshold >= 100)
    stop("`detection_threshold` must be between 0% and 100%.")
  if (!is.logical(remove_no_oral) || length(remove_no_oral) != 1 || is.na(remove_no_oral))
    stop("`remove_no_oral` must be TRUE or FALSE.")

  # Remove 0 columns, Samples with no bacteria detected.
  n_samples_input <- ncol(species_table)
  empty <- colSums(species_table) == 0
  if (any(empty)) {
    warning(sum(empty), " sample(s) with zero total abundance removed: ",
            paste(colnames(species_table)[empty], collapse = ", "))
    species_table <- species_table[, !empty, drop = FALSE]
  }
  if (ncol(species_table) == 0) stop("All samples have zero total abundance.")

  # Relative abundance (RA):
  rel_ab <- sweep(species_table, 2, colSums(species_table), "/")

  # Setting up Detection threshold AFTER RA is calculated.
  # It is then re-scaled up to 1 again.
  below_threshold <- rel_ab > 0 & rel_ab < detection_threshold / 100
  if (any(below_threshold)) {
    rel_ab[below_threshold] <- 0
    if (any(colSums(rel_ab) == 0))
      stop("`detection_threshold` is so high that some samples have nothing left.")
    rel_ab <- sweep(rel_ab, 2, colSums(rel_ab), "/")
  }

  ###############################################################################
  ############## Substitute the oral taxa into the species table ################
  ###############################################################################
  species_table_taxa  <- rownames(species_table)

  ### Clean up rownames:
  cleaned_taxa_name <- sub(".*\\|", "", species_table_taxa)
  cleaned_taxa_name <- sub("^[a-z]__", "", cleaned_taxa_name)
  cleaned_taxa_name <- trimws(gsub("[_ ]+", " ", cleaned_taxa_name))

  # Names are compared ignoring capitals and punctuation, so "Veillonella_sp_T11011_6"
  # (MetaPhlAn) matches "Veillonella sp. T11011-6" (NCBI).
  name_key <- function(x) tolower(trimws(gsub("[^A-Za-z0-9]+", " ", x)))

  species_ids         <- rep(NA_character_, length(cleaned_taxa_name))

  # Detection process:
  #          - First check if it is NCBI ID already (save time)
  #          - Second check if species name match with reference database.
  #          - Third check OraCol's NCBI name list (NCBI_lookup). No internet is needed.

  # Step 1: row names that are only digits are already NCBI IDs.
  is_id <- grepl("^[0-9]+$", species_table_taxa)
  species_ids[is_id] <- species_table_taxa[is_id]

  # Step 2: names in the chosen oral database (current, new or old names): take the ID from it.
  old_names   <- strsplit(Oral_taxa$old_names, ";")
  known_names <- c(Oral_taxa$name, Oral_taxa$new_name, unlist(old_names))
  known_ids   <- c(Oral_taxa$ncbi_id, Oral_taxa$ncbi_id, rep(Oral_taxa$ncbi_id, lengths(old_names)))
  todo <- is.na(species_ids)
  species_ids[todo] <- known_ids[match(name_key(cleaned_taxa_name[todo]), name_key(known_names))]

  # Step 3: everything else is looked up in OraCol's NCBI name list (no internet needed).
  todo <- is.na(species_ids)
  species_ids[todo] <- NCBI_lookup$ncbi_id[match(name_key(cleaned_taxa_name[todo]),
                                                 name_key(NCBI_lookup$name))]
  if (anyNA(species_ids))
    message(sum(is.na(species_ids)), " of ", length(species_ids), " row name(s) are not in OraCol's ",
            "NCBI name list and count as not oral (NCBI_ID is NA in $matching).")

  is_oral <- !is.na(species_ids) & species_ids %in% oral_taxa_id

  if (!any(is_oral)) {
    stop("No taxa in `species_table` matched the oral definition. ")}

  # Matching dataframe for rechecking
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

  # Samples with no oral species detected: their oral CLR values only reflect the rest of
  # the community (every oral zero gets the same small negative CLR), so by default they are
  # moved out of the sample summary and listed in $no_oral_samples instead.
  no_oral <- Summary_table_samples$Oral_richness == 0
  no_oral_samples <- Summary_table_samples[no_oral, , drop = FALSE]
  rownames(no_oral_samples) <- NULL
  if (any(no_oral)) {
    if (remove_no_oral) {
      message(sum(no_oral), " sample(s) with no oral species detected were removed from ",
              "$Summary_table_samples and listed in $no_oral_samples.")
      Summary_table_samples <- Summary_table_samples[!no_oral, , drop = FALSE]
      rownames(Summary_table_samples) <- NULL
      if (nrow(Summary_table_samples) == 0)
        warning("No sample has any oral species detected, so $Summary_table_samples is empty.")
    } else {
      message(sum(no_oral), " sample(s) have no oral species detected (listed in $no_oral_samples).")
    }
  }

  results_object = list("species_table_oral" = species_table_oral,
                        "species_table_oral_CLR" = species_table_oral_CLR,
                        "Summary_table_species" = Summary_table_species,
                        "Summary_table_samples" = Summary_table_samples,
                        "matching" = matching,
                        "no_oral_samples" = no_oral_samples)
  return(results_object)
}

