#' Oral reference databases used by OraCol
#'
#' The seven oral definitions that [OraCol()] can use: four published oral databases
#' (`HOMD`, `eHOMD`, `OES`, `mBodyMap`) and three that combine eHOMD and OES (`eHOMD_OES`,
#' `eHOMD_OES_oral`, `eHOMD_OES_oropharyngeal`). Each one lists the taxa that the database
#' considers oral. `HOMD` is the default.
#'
#' @format A named list of seven data frames (`OES`, `mBodyMap`, `HOMD`, `eHOMD`,
#' `eHOMD_OES`, `eHOMD_OES_oropharyngeal`, `eHOMD_OES_oral`), each with one row per taxon
#' and these columns:
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
#' **Combined eHOMD + OES databases**
#'
#' In addition, for more precise detection, eHOMD and OES are combined into one database
#' (`eHOMD_OES`) with two curated versions: `eHOMD_OES_oral` for mainly oral microbes, and
#' `eHOMD_OES_oropharyngeal` for microbes of the mouth and nasal passages. The curated
#' versions remove "generalists": species that are common at several body sites but whose
#' main niche is somewhere else (see **How generalists were found** below).
#'
#' **eHOMD_OES**: eHOMD and OES combined
#'
#' * DB size: 525 taxa
#' * Built from: merging OES and eHOMD.
#' * Status: Healthy and diseases
#'
#' **eHOMD_OES_oral**: mainly oral microbes
#'
#' * DB size: 504 taxa
#' * Built from: `eHOMD_OES` minus the 21 generalists whose main niche is not the mouth
#' * Status: Healthy and diseases
#'
#' **eHOMD_OES_oropharyngeal**: microbes of the mouth and nasal passages
#'
#' * DB size: 517 taxa
#' * Built from: `eHOMD_OES` minus the 8 generalists whose main niche is the gut.
#' * Status: Healthy and diseases
#'
#' **How generalists were found**
#'
#' Every eHOMD and OES species was checked in curatedMetagenomicData (CMD),
#' using all its samples from five body sites: oral cavity (857), gut (21,030), nasal cavity
#' (93), skin (504) and vagina (96). Milk was left out because it has only 8 samples.
#'
#' * Detection threshold: relative abundance above 0.1% in a sample.
#' * Prevalence of generalists: prevalence of at least 1% at 3 or more of the 5 body sites.
#' * Niche: the body site where its prevalence is highest.
#'
#' @source The references for each database are given under Details. Names were taken from
#' the NCBI Taxonomy database. Generalists were found with curatedMetagenomicData.
#'
#' @examples
#' data("Oral_db", package = "OraCol")
#' names(Oral_db)
#' sapply(Oral_db, nrow)   # number of taxa in each database
#' head(Oral_db$HOMD)
#'
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
