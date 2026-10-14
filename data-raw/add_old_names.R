##############################################################################################
# add_old_names.R
#
# Adds two columns to every oral list in Oral_db_v3.rds and saves the result as Oral_db_v4.rds:
#
#     new_name    the current NCBI scientific name of that NCBI ID
#     old_names   earlier names of the same taxon (NCBI "synonyms"), separated by ";"
#                 e.g. Schaalia odontolytica  ->  "Actinomyces odontolyticus"
#
# With these, a species table that still uses an old name (older MetaPhlAn / Kraken databases)
# is matched by name, with no internet lookup.
#
# Needs the internet (4 requests to NCBI) and the xml2 package. xml2 is only used by this
# script; it is NOT a dependency of OraCol. Run from the package root, then rebuild sysdata:
#
#   source("data-raw/add_old_names.R")
#   source("data-raw/make_Oral_db.R")      # after pointing it at Oral_db_v4.rds
##############################################################################################

Oral_db <- lapply(readRDS("Oral_db_v3.rds"), as.data.frame, stringsAsFactors = FALSE)
all_ids <- unique(unlist(lapply(Oral_db, function(d) d$ncbi_id)))

###################### 1. Fetch every ID's names from NCBI ##################################
new_name  <- character(0)   # NCBI ID -> current scientific name
synonyms  <- list()         # NCBI ID -> names NCBI lists as synonyms (raw text)
id_update <- character(0)   # outdated NCBI ID -> the ID NCBI merged it into

for (batch in split(all_ids, ceiling(seq_along(all_ids) / 200))) {
  url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&retmode=xml&id=",
                paste(batch, collapse = ","))
  for (taxon in xml2::xml_find_all(xml2::read_xml(url), "/TaxaSet/Taxon")) {
    id <- xml2::xml_text(xml2::xml_find_first(taxon, "./TaxId"))
    new_name[id]   <- xml2::xml_text(xml2::xml_find_first(taxon, "./ScientificName"))
    synonyms[[id]] <- xml2::xml_text(xml2::xml_find_all(taxon, "./OtherNames/Synonym"))
    for (old_id in xml2::xml_text(xml2::xml_find_all(taxon, "./AkaTaxIds/TaxId"))) id_update[old_id] <- id
  }
  Sys.sleep(1)   # be polite to NCBI
}

###################### 2. Turn NCBI's synonym text into clean old names #####################
# Keep a synonym only if it is a proper "Genus species" (or "Genus species subsp. x") name,
# alone or followed by author text, e.g. "Prevotella oralis (Loesche et al. 1964) Shah and
# Collins 1990" -> "Prevotella oralis". Anything else (strain labels, "sp.", "genomosp. C1",
# "group II" ...) gives NA and is dropped.
bare_name <- function(x) {
  pattern <- "^((Candidatus )?\\[?[A-Z][a-z]+\\]? [a-z-]+( subsp\\. [a-z-]+)?)( corrig\\.)?( (\\(|[A-Z]).*)?$"
  ifelse(grepl(pattern, x), sub(pattern, "\\1", x), NA_character_)
}
old_names <- lapply(names(synonyms), function(id) {
  old <- unique(bare_name(synonyms[[id]]))
  old[!is.na(old) & tolower(old) != tolower(new_name[id])]
})
names(old_names) <- names(synonyms)

# Drop ambiguous old names: used for two different IDs, or equal to another ID's current name.
pairs <- data.frame(id = rep(names(old_names), lengths(old_names)), old = tolower(unlist(old_names)))
n_ids <- tapply(pairs$id, pairs$old, function(v) length(unique(v)))
ambiguous <- union(names(n_ids)[n_ids > 1], intersect(pairs$old, tolower(new_name)))
old_names <- lapply(old_names, function(old) old[!tolower(old) %in% ambiguous])

###################### 3. Add the two columns to every list and save #######################
Oral_db <- lapply(Oral_db, function(d) {
  current_id  <- ifelse(d$ncbi_id %in% names(new_name), d$ncbi_id, id_update[d$ncbi_id])
  d$new_name  <- unname(new_name[current_id])
  d$old_names <- vapply(current_id, function(id) paste(old_names[[id]], collapse = ";"), character(1))
  d
})
saveRDS(Oral_db, "Oral_db_v4.rds")

for (n in names(Oral_db))
  message(sprintf("%-9s %3d rows, %3d with at least one old name", n, nrow(Oral_db[[n]]),
                  sum(nzchar(Oral_db[[n]]$old_names))))
message("Ambiguous old names dropped: ", paste(ambiguous, collapse = ", "))
outdated <- intersect(all_ids, names(id_update))
if (length(outdated))
  message("Outdated NCBI IDs in the lists (NCBI merged them): ",
          paste(sprintf("%s -> %s", outdated, id_update[outdated]), collapse = ", "))
