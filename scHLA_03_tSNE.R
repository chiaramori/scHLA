#############################################################
# Section 1: Setup & Dependencies
#############################################################

library(data.table)
library(dplyr)
library(readr)
library(stringr)

# Paths
input_dir <- "MINRUN025_scHLA_inserts/trimmed/trimmed_selected_4tSNE"  # peptide input files
output_dir <- "MINRUN025_scHLA_inserts/trimmed/trimmed_selected_4tSNE_blasted"  # for saving annotated versions
ref_db_path <- "databases/epitope_full_v3_simplified.tsv"  # cleaned reference

if (!dir.exists(output_dir)) dir.create(output_dir)

#############################################################
# Section 2: Load and Prepare Reference Epitope DB
#############################################################
# The database was pre-cleaned:
#   - Column2_Name = peptide AA sequence
#   - Column9_Source_Molecule = protein
#   - Column13_Source_Organism = species

ref_db <- read_tsv(ref_db_path, show_col_types = FALSE)

cols_to_select <- c("Column2_Name", "Column9_Source_Molecule", "Column13_Source_Organism")

ref_db_simple <- dplyr::select(ref_db, dplyr::all_of(cols_to_select)) %>%
  dplyr::distinct() %>%
  dplyr::rename(
    Peptide = Column2_Name,
    SourceMolecule = Column9_Source_Molecule,
    Species = Column13_Source_Organism
  )


#############################################################
# Section 3: Loop Through Peptide Files and Annotate
#############################################################

input_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(input_files) == 0) stop("❌ No peptide CSVs found in input directory!")

for (f in input_files) {
  cat("\n🔍 Processing:", basename(f), "\n")
  df <- fread(f)
  
  if (!("amino_acid" %in% colnames(df))) {
    cat("⚠️ Skipping: no 'amino_acid' column found.\n")
    next
  }
  
  # Left join on peptide sequence to bring in species & molecule info
  annotated_df <- df %>%
    left_join(ref_db_simple, by = c("amino_acid" = "Peptide"))
  
  # Save _blasted file (peptide + species + source molecule)
  out_file <- file.path(output_dir, str_replace(basename(f), "\\.csv$", "_blasted.csv"))
  fwrite(annotated_df, out_file)
  cat("✅ Saved:", out_file, "\n")
  
}

cat("\n🎉 All done! Annotated and species files are in:", output_dir, "\n")

#############################################################
# Section 4: Modular and Independent tSNE & UMAP Execution
#############################################################
# --- Load libraries ---
library(dplyr)
library(ggplot2)
library(viridis)
library(umap)
library(Rtsne)
library(readr)

# --- Define directories ---
input_dir <- "MINRUN025_scHLA_inserts/trimmed/trimmed_selected_4tSNE_blasted"
dist_dir <- "MINRUN025_scHLA_inserts/distances"
tsne_dir <- "MINRUN025_scHLA_inserts/tsne"
umap_dir <- "MINRUN025_scHLA_inserts/umap"
plot_dir <- "MINRUN025_scHLA_inserts/plots"

# Create output directories if they don't exist
dir.create(dist_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tsne_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(umap_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

max_peptides <- 1000

# --- Helper: Compute and save distance matrix for one file ---
compute_distance_for_file <- function(file, dist_dir) {
  data <- read_csv(file)
  peptides <- data$amino_acid
  n <- length(peptides)
  dist_mat <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n) {
    for(j in i:n) {
      dist <- sum(strsplit(peptides[i], "")[[1]] != strsplit(peptides[j], "")[[1]])
      dist_mat[i,j] <- dist
      dist_mat[j,i] <- dist
    }
  }
  rownames(dist_mat) <- peptides
  colnames(dist_mat) <- peptides
  dist_file <- file.path(dist_dir, paste0(tools::file_path_sans_ext(basename(file)), "_dist.rds"))
  saveRDS(dist_mat, dist_file)
}

# --- Run distance computation for all files ---
run_compute_distances_all <- function(input_dir, dist_dir) {
  blasted_files <- list.files(input_dir, pattern = "_blasted\\.csv$", full.names = TRUE)
  if(length(blasted_files) == 0) stop("❌ No blasted CSV files found in input directory!")
  for(file in blasted_files) {
    message("▶ Computing distances for: ", basename(file))
    tryCatch({
      compute_distance_for_file(file, dist_dir)
      message("✅ Done computing distances for ", basename(file))
    }, error = function(e) {
      message("❌ Error computing distances for ", basename(file), ": ", e$message)
    })
  }
}

# --- Helper: Run tSNE for one file ---
run_tsne_for_file <- function(blasted_file) {
  base_name <- tools::file_path_sans_ext(basename(blasted_file))
  dist_file <- file.path(dist_dir, paste0(base_name, "_dist.rds"))
  if(!file.exists(dist_file)) stop("Distance file not found for: ", base_name)
  dist_mat <- readRDS(dist_file)
  
  blast_df <- read_csv(blasted_file) %>%
    slice_max(order_by = count, n = max_peptides) %>%
    arrange(count)
  
  peptides <- blast_df$amino_acid
  dist_mat <- dist_mat[peptides, peptides]
  
  set.seed(42)
  n <- nrow(dist_mat)
  perplexity <- min(30, floor((n - 1) / 3))
  if (perplexity < 1) stop("Too few samples for tSNE: ", n)
  
  tsne_out <- Rtsne(as.dist(dist_mat), is_distance = TRUE, dims = 2, perplexity = perplexity, verbose = FALSE)
  
  
  emb_df <- data.frame(
    amino_acid = peptides,
    TSNE1 = tsne_out$Y[,1],
    TSNE2 = tsne_out$Y[,2]
  )
  
  merged <- emb_df %>%
    left_join(blast_df, by = "amino_acid") %>%
    mutate(
      TotalCount = count,
      HumanNonHuman = case_when(
        Species == "Homo sapiens" ~ "Human",
        !is.na(Species) ~ "Non-human",
        TRUE ~ NA_character_
      )
    ) %>%
    arrange(TotalCount)
  
  tsne_file <- file.path(tsne_dir, paste0(base_name, "_tsne.csv"))
  write_csv(merged, tsne_file)
  message("✅ tSNE saved for ", base_name)
  return(merged)
}


# --- Run tSNE for all files ---
run_tsne_all <- function() {
  blasted_files <- list.files(input_dir, pattern = "_blasted\\.csv$", full.names = TRUE)
  if(length(blasted_files) == 0) stop("❌ No blasted CSV files found in input directory!")
  for(file in blasted_files) {
    message("▶ Running tSNE for: ", basename(file))
    tryCatch({
      merged <- run_tsne_for_file(file)
      plot_embedding(merged, "tsne", tools::file_path_sans_ext(basename(file)), plot_dir)
    }, error = function(e) {
      message("❌ Error in tSNE for file: ", basename(file))
      message(e$message)
    })
  }
}

# --- Helper: Run UMAP for one file ---
run_umap_for_file <- function(blasted_file) {
  base_name <- tools::file_path_sans_ext(basename(blasted_file))
  dist_file <- file.path(dist_dir, paste0(base_name, "_dist.rds"))
  if(!file.exists(dist_file)) stop("Distance file not found for: ", base_name)
  dist_mat <- readRDS(dist_file)
  
  blast_df <- read_csv(blasted_file) %>%
    slice_max(order_by = count, n = max_peptides) %>%
    arrange(count)
  
  peptides <- blast_df$amino_acid
  dist_mat <- dist_mat[peptides, peptides]
  
  set.seed(42)
  n <- nrow(dist_mat)
  if (n < 5) {
    warning("Skipping file ", base_name, ": too few peptides (", n, ") for UMAP.")
    return(NULL)
  }
  n_neighbors <- min(15, n - 1)
  if (n_neighbors < 2) {
    warning("Skipping file ", base_name, ": n_neighbors too small (", n_neighbors, ")")
    return(NULL)
  }
  
  umap_config <- umap::umap.defaults
  umap_config$n_neighbors <- n_neighbors
  
  umap_out <- umap::umap(as.matrix(dist_mat), config = umap_config)
  
  
  emb_df <- data.frame(
    amino_acid = peptides,
    UMAP1 = umap_out$layout[,1],
    UMAP2 = umap_out$layout[,2]
  )
  
  merged <- emb_df %>%
    left_join(blast_df, by = "amino_acid") %>%
    mutate(
      TotalCount = count,
      HumanNonHuman = case_when(
        Species == "Homo sapiens" ~ "Human",
        !is.na(Species) ~ "Non-human",
        TRUE ~ NA_character_
      )
    ) %>%
    arrange(TotalCount)
  
  umap_file <- file.path(umap_dir, paste0(base_name, "_umap.csv"))
  write_csv(merged, umap_file)
  message("✅ UMAP saved for ", base_name)
  return(merged)
}


# --- Run UMAP for all files ---
run_umap_all <- function() {
  blasted_files <- list.files(input_dir, pattern = "_blasted\\.csv$", full.names = TRUE)
  if(length(blasted_files) == 0) stop("❌ No blasted CSV files found in input directory!")
  for(file in blasted_files) {
    message("▶ Running UMAP for: ", basename(file))
    tryCatch({
      merged <- run_umap_for_file(file)
      plot_embedding(merged, "umap", tools::file_path_sans_ext(basename(file)), plot_dir)
    }, error = function(e) {
      message("❌ Error in UMAP for file: ", basename(file))
      message(e$message)
    })
  }
}

# --- Helper: Plot embedding colored by count, species, and human/non-human ---
plot_embedding <- function(merged_df, method, base_name, out_dir) {
  x_col <- paste0(toupper(substr(method, 1, 4)), "1")
  y_col <- paste0(toupper(substr(method, 1, 4)), "2")
  
  p_count <- ggplot(merged_df, aes(.data[[x_col]], .data[[y_col]], color = TotalCount)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_viridis_c(na.value = "grey") +
    theme_bw() +
    ggtitle(paste0(toupper(method), " - Colored by Count"))
  
  p_species <- ggplot(merged_df, aes(.data[[x_col]], .data[[y_col]], color = Species)) +
    geom_point(size = 2, alpha = 0.8) +
    theme_bw() +
    ggtitle(paste0(toupper(method), " - Colored by Species"))
  
  p_human <- ggplot(merged_df, aes(.data[[x_col]], .data[[y_col]], color = HumanNonHuman)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_manual(values = c("Human" = "blue", "Non-human" = "red"), na.value = "grey") +
    theme_bw() +
    ggtitle(paste0(toupper(method), " - Colored by Human/Non-human"))
  
  ggsave(filename = file.path(out_dir, paste0(base_name, "_", method, "_count.png")), plot = p_count, width = 6, height = 5)
  ggsave(filename = file.path(out_dir, paste0(base_name, "_", method, "_species.png")), plot = p_species, width = 6, height = 5)
  ggsave(filename = file.path(out_dir, paste0(base_name, "_", method, "_human.png")), plot = p_human, width = 6, height = 5)
  
  message("✅ Plots saved for ", base_name, " (", toupper(method), ")")
}

# 1. Compute all distance matrices
run_compute_distances_all(input_dir, dist_dir)

# 2. Run tSNE
run_tsne_all()

# 3. Run UMAP
run_umap_all()
