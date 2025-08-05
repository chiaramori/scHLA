# Load required packages
library(readr)
library(dplyr)

# Load the raw TSV file (assuming it's tab-separated)
file_path <- "databases/epitope_full_v3.tsv"  # replace with your actual file path
raw_data <- read_tsv(file_path, col_names = FALSE)

# Extract the first two rows as headers
header1 <- as.character(raw_data[1, ])
header2 <- as.character(raw_data[2, ])

# Process headers: replace spaces with underscores instead of dashes
new_headers <- sapply(seq_along(header1), function(i) {
  h1 <- gsub(" ", "_", header1[i])  # Replace spaces with underscores
  h2 <- gsub(" ", "_", header2[i])  # Replace spaces with underscores
  
  if (i == 1) {
    # First column is "EpitopeID" (remove any underscores from h1 to keep it clean)
    h1 <- gsub("_", "", h1)         
    paste0(h1, "_", h2)
  } else {
    # Other columns: rename as ColumnX_Header
    h1 <- paste0("Column", i - 1)
    paste0(h1, "_", h2)
  }
})

# Drop the first two rows (header rows) and assign new column names
clean_data <- raw_data[-c(1, 2), ]
colnames(clean_data) <- new_headers

# Optionally, write to a new cleaned file
write_tsv(clean_data, "databases/epitope_full_v3_cleaned.tsv")

# Reload cleaned data for further processing
df <- read_tsv("databases/epitope_full_v3_cleaned.tsv")

#---- SIMPLIFIED PROCESSING ----

# Store the rows where Column9_Source_Molecule is NA
na_mask <- is.na(df$Column9_Source_Molecule)

# Concatenate safely with underscores if both are present, else just use what's present
concatenated_value <- ifelse(
  na_mask,
  paste0(
    coalesce(df$Column17_Epitope_Relation, ""),
    ifelse(!is.na(df$Column17_Epitope_Relation) & !is.na(df$Column24_Source_Molecule), "_", ""),
    coalesce(df$Column24_Source_Molecule, "")
  ),
  df$Column9_Source_Molecule
)

# Apply the replacements
df <- df %>%
  mutate(
    Column9_Source_Molecule = concatenated_value,
    Column10_Source_Molecule_IRI = if_else(
      na_mask,
      Column25_Source_Molecule_IRI,
      Column10_Source_Molecule_IRI
    ),
    Column13_Source_Organism = if_else(
      na_mask,
      Column28_Source_Organism,
      Column13_Source_Organism
    )
  )

# Save simplified output
write_tsv(df, "databases/epitope_full_v3_simplified.tsv")

