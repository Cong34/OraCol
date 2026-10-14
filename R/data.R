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
#' The four databases:
#'
#' | Database | Oral taxa list | Reference |
#' |---|---|---|
#' | `HOMD` | Human Oral Microbiome Database | Chen et al. 2010 |
#' | `eHOMD` | expanded HOMD (Add Nasal) | Escapa et al. 2018 |
#' | `OES` | Oral Enrichment Score at 1% | Manghi, 2025 |
#' | `mBodyMap` | mBodyMap_Oral Healthy | Jin et al. 2022 |
#'
#' More details on how each one was built:
#'
#' | Database | DB size | Samples size | Sequencing method | Species definition | Status |
#' |---|---|---|---|---|---|
#' | `HOMD` | 340 taxa | 33k 16S RNA clone | Full length 16S | 98.5% 16S identity | Healthy and diseases |
#' | `eHOMD` | 404 taxa | 44k 16S RNA clone | Full length 16S | 98.5% 16S identity | Healthy and diseases |
#' | `OES` | 305 taxa | 857 metagenomic sample | Shotgun | 1% threshold prevalence | Healthy and diseases |
#' | `mBodyMap` | 351 taxa | 8077 samples | 16S + Shotgun | NCBI taxonomy (>=2 samples + 0.01% RA) | Healthy |
#'
#' @source HOMD: Chen et al. 2010. eHOMD: Escapa et al. 2018. OES: Manghi 2025.
#' mBodyMap: Jin et al. 2022. Names were taken from the NCBI Taxonomy database.
#'
#' @examples
#' data("Oral_db", package = "OraCol")
#' names(Oral_db)
#' head(Oral_db$HOMD)
"Oral_db"

# Oral_db is loaded with data() inside OraCol(). This line only tells R's package checker
# that the name is expected, so it does not report it as an undefined variable.
utils::globalVariables("Oral_db")
