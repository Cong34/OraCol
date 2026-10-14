# Function for: Cleaning the species name: remove the "_" and other metaphlan format
clean_name <- function(x) {
  x <- sub(".*\\|", "", x)
  x <- sub("^[a-z]__", "", x)
  trimws(gsub("[_ ]+", " ", x))
}
