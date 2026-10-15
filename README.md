# OraCol

Summarise, detect and categorise Oral Colonisation of bacteria in the gut. 

OraCol finds which species in a metagenomic species table are oral bacteria, using one of four
curated oral taxa reference databases, and returns per-species and per-sample summaries: relative abundance,
prevalence, centred log-ratio (CLR) values, richness and Shannon diversity.

## Installation

R version 4.5.0 or newer. 

Dependency: base R
```r
install.packages("remotes")
remotes::install_github("Cong34/OraCol")
```

## Quick start

```r
library(OraCol)

# Example table: 
example_data <- system.file("extdata", "sampling_data.csv", package = "OraCol")
species_table <- read.csv(example_data, row.names = 1, check.names = FALSE)

OraCol_res <- OraCol(species_table, oral_definition = "HOMD")

# The main results of OraCol: 
OraCol_res$Summary_table_samples   # one row per sample
OraCol_res$Summary_table_species   # one row per oral species

# It also currently output supplementary data, check OraCol_res other output. 
```

## Input

- A matrix or data.frame of **species-level** abundances. Counts, proportions and percentages
  all give the same results. (Also accept NCBI ID as input)
- Relative abundances below `detection_threshold` (default 0.001%) count as not detected and are
  set to 0, so a raw table and an already-filtered table give the same results. Use
  `detection_threshold = 0` to keep every value.

Species names are matched to NCBI taxonomy IDs with TaxSEA. Names missing from TaxSEA's built-in
list are looked up online at NCBI, which is slower and needs an internet connection.

## Oral definitions

| `oral_definition` | Oral taxa list | Reference |
|---|---|---|
| `"HOMD"` | Human Oral Microbiome Database | Chen et al. 2010 |
| `"eHOMD"` | expanded HOMD (Add Nasal) | Escapa et al. 2018 |
| `"OES"` | Oral Enrichment Score at 1% | Manghi, 2025 |
| `"mBodyMap"` | mBodyMap_Oral Healthy | Jin et al. 2022 |

## More Details: 
| Oral_db | DB size | Samples size | Sequencing method | Species definition | Status |
|---|---|---|---|---|---|
| HOMD | 340 taxa | 33k 16S RNA clone | Full length 16S | 98.5% 16S identity | Healthy and diseases |
| eHOMD | 404 taxa | 44k 16S RNA clone | Full length 16S | 98.5% 16S identity | Healthy and diseases |
| OES | 305 taxa | 857 metagenomic sample | Shotgun | 1% threshold prevalence | Healthy and diseases |
| mBodyMap | 351 taxa | 8077 samples | 16S + Shotgun | NCBI taxonomy (>=2 samples + 0.01% RA) | Healthy |

The databases ship with the package, so you can look at exactly which taxa each one contains:

```r
data("Oral_db", package = "OraCol")
names(Oral_db)       # "OES" "mBodyMap" "HOMD" "eHOMD"
?Oral_db             # description of the columns
```

## Output

`OraCol()` returns a list:

| Element | Contents |
|---|---|
| `species_table_oral` | Relative abundance of the oral species in each sample |
| `species_table_oral_CLR` | Center Log ratio values of the oral species |
| `Summary_table_species` | Summary of species: relative abundance, prevalence, number of samples detected, and CLR |
| `Summary_table_samples` | Summary of sample: oral relative abundance, oral CLR sum, richness and Shannon diversity |
| `matching` | Every input species, its NCBI ID (`NA` if not found) and whether it counted as oral |



## Example: oral samples from the Human Microbiome Project

This example runs OraCol on the HMP_2012 dataset from
[curatedMetagenomicData](https://bioconductor.org/packages/curatedMetagenomicData)

```r
library(OraCol)
library(curatedMetagenomicData)

# Species relative abundance of HMP_2012 (downloaded the first time, then stored on your computer)
cmd <- curatedMetagenomicData(pattern = "2021-03-31.HMP_2012.relative_abundance",
                                                 counts = TRUE, dryrun = FALSE)
# Count table: 
count_table <- assay(cmd$`2021-03-31.HMP_2012.relative_abundance`)
species_table <- species_table[rowSums(species_table) > 0, ]  

# Optional: shorten the labels (OraCol already do this)
# rownames(species_table) <- sub(".*s__", "", rownames(species_table))

HMP_res <- OraCol(species_table, oral_definition = "HOMD")

head(HMP_res$Summary_table_samples)   # per sample: oral share, richness and diversity
head(HMP_res$Summary_table_species)   # per species: sort by prevalence
```






## License

MIT © Cong Pham
