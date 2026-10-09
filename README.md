## OraCol

Summarise, detect and categorise Oral Colonisation of bacteria in the gut. 

OraCol finds which species in a metagenomic species table are oral bacteria, using one of four
curated oral taxa reference databases, and returns per-species and per-sample summaries: relative abundance,
prevalence, centred log-ratio (CLR) values, richness and Shannon diversity.

## Installation

R version 4.5.0 or newer. 

Dependency: [TaxSEA](https://bioconductor.org/packages/TaxSEA) (from Bioconductor).

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
- Species names as row names, written `Genus_species` (e.g. `Streptococcus_mitis`), and sample
  names as column names. If species are columns instead, use `taxa_are_rows = FALSE`.
- Otherwise, provide NCBI ID for row also works wonderfully. 
- No missing (`NA`) or negative values.

Species names are matched to NCBI taxonomy IDs with TaxSEA. Names missing from TaxSEA's built-in
list are looked up online at NCBI, which is slower and needs an internet connection.

## Oral definitions

| `oral_definition` | Oral taxa list | Reference |
|---|---|---|
| `"HOMD"` (default) | Human Oral Microbiome Database | Chen et al. 2010 |
| `"eHOMD"` | expanded HOMD (adds nasal and upper-airway taxa) | Escapa et al. 2018 |
| `"OES"` | Oral Enrichment Score at 1% | Manghi, 2025 |
| `"mBodyMap"` | mBodyMap_Oral Healthy only| Jin et al. 2022 |

## Output

`OraCol()` returns a list:

| Element | Contents |
|---|---|
| `species_table_oral` | Relative abundance (0-1) of the oral species in each sample |
| `species_table_oral_CLR` | CLR values of the oral species |
| `Summary_table_species` | One row per oral species: mean, median, SD and sum of relative abundance, prevalence, number of samples detected, mean and median CLR (sorted by prevalence) |
| `Summary_table_samples` | One row per sample: oral relative abundance, oral CLR sum and mean, oral and whole-community richness and Shannon diversity |
| `matching` | Every input species, its NCBI ID (`NA` if not found) and whether it counted as oral |


## License

MIT © Cong Pham
