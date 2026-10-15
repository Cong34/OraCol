#' Oral reference databases used by OraCol
#'
#' The four oral definitions that [OraCol()] can use. Each one lists the taxa that the
#' database considers oral.
#'
#' @format A named list of four data frames (`OES`, `mBodyMap`, `HOMD`, `eHOMD`), each with
#' one row per taxon and these columns:
#' * `ncbi_id`: NCBI taxonomy ID.
#' * `name`: scientific name.
#' * `rank`: taxonomic rank, for example `species` or `genus`.
#' * `Taxa`: the name with underscores instead of spaces.
#' * `new_name`: the current NCBI scientific name of that ID.
#' * `old_names`: earlier names of the same taxon, separated by `;`.
#'
#' @details
#' **HOMD**: Human Oral Microbiome Database (Chen et al. 2010)
#'
#' * DB size: 340 taxa
#' * Samples size: 33k 16S RNA clone
#' * Sequencing method: Full length 16S
#' * Species definition: 98.5% 16S identity
#' * Status: Healthy and diseases
#'
#' **eHOMD**: expanded HOMD (Add Nasal) (Escapa et al. 2018)
#'
#' * DB size: 404 taxa
#' * Samples size: 44k 16S RNA clone
#' * Sequencing method: Full length 16S
#' * Species definition: 98.5% 16S identity
#' * Status: Healthy and diseases
#'
#' **OES**: Oral Enrichment Score at 1% (Manghi, 2025)
#'
#' * DB size: 305 taxa
#' * Samples size: 857 metagenomic sample
#' * Sequencing method: Shotgun
#' * Species definition: 1% threshold prevalence
#' * Status: Healthy and diseases
#'
#' **mBodyMap**: mBodyMap_Oral Healthy (Jin et al. 2022)
#'
#' * DB size: 351 taxa
#' * Samples size: 8077 samples
#' * Sequencing method: 16S + Shotgun
#' * Species definition: NCBI taxonomy (>=2 samples + 0.01% RA)
#' * Status: Healthy
#'
#' @source The references for each database are given under Details. Names were taken from
#' the NCBI Taxonomy database.
#'
#' @examples
#' data("Oral_db", package = "OraCol")
#' names(Oral_db)
#' head(Oral_db$HOMD)
"Oral_db"


#' OraCol's list of taxon names and their NCBI taxonomy IDs
#'
#' [OraCol()] uses this list to find the NCBI ID of any name that is not in the chosen oral
#' database, so it never needs an internet connection. The list was built once by asking NCBI
#' about every name in it, and can be extended with the script
#' `data-raw/update_NCBI_lookup.R` in the package's GitHub repository.
#'
#' @format A data frame with one row per name and these columns:
#' * `name`: the taxon name, cleaned the same way as `OraCol()` cleans row names
#'   (lineage and `s__` prefix removed, underscores replaced by spaces).
#' * `ncbi_id`: its NCBI taxonomy ID, or `NA` if NCBI did not recognise the name.
#' * `ncbi_name`: the current NCBI scientific name for that ID. It differs from `name` when
#'   the taxon has been renamed.
#' * `rank`: the NCBI rank, for example `species` or `genus`.
#' * `source`: where the name came from: TaxSEA's name list (only the names were used; every
#'   ID comes from NCBI) or the oral samples of curatedMetagenomicData.
#'
#' @source NCBI Taxonomy, queried with the NCBI Datasets service.
#'
#' @examples
#' data("NCBI_lookup", package = "OraCol")
#' nrow(NCBI_lookup)
#' head(NCBI_lookup[!is.na(NCBI_lookup$ncbi_id), ])
"NCBI_lookup"

# Oral_db and NCBI_lookup are loaded with data() inside OraCol(). This line only tells R's
# package checker that the names are expected, so it does not report them as undefined.
utils::globalVariables(c("Oral_db", "NCBI_lookup"))
