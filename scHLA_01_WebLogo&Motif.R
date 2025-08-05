# =============================
# Peptide Motif Analysis Script (Flat Execution)
# Description:
# - Processes 11aa insert sequences from barcode folders
# - Trims and saves to `trimmed/`
# - Plots sequence logos
# - Classifies by AA class motif
# - Computes motif frequencies and plots barplots
# - Computes and plots dendrograms of motifs
# =============================

# === Load Libraries ===
library(ggseqlogo)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(stringdist)
library(purrr)
library(dendextend)

# === Prepare Output Folders ===
dir.create("trimmed", showWarnings = FALSE)
dir.create("trimmed/MINRUN025_motif", showWarnings = FALSE)
dir.create("cluster_motif", showWarnings = FALSE)
dir.create("cluster_motif/dendrogram", showWarnings = FALSE)

# === STEP 1: Trim inserts and save trimmed CSVs ===
barcodes <- list.dirs(".", recursive = FALSE, full.names = FALSE)
file_paths <- file.path(barcodes, "11aa_inserts.csv")
names(file_paths) <- barcodes

valid_barcodes <- names(file_paths)[file.exists(file_paths)]
sample_dfs <- lapply(valid_barcodes, function(barcode) {
  df <- read.csv(file_paths[barcode], stringsAsFactors = FALSE)
  df <- df %>% select(amino_acid, count)
  df <- df[!grepl("\\*", df$amino_acid), ]  # Remove sequences with asterisks
  df$amino_acid <- substr(df$amino_acid, 2, nchar(df$amino_acid) - 1)  # Trim to 9-mer
  df$barcode <- barcode
  df
})

for (i in seq_along(valid_barcodes)) {
  out_path <- file.path("trimmed", paste0(valid_barcodes[i], ".csv"))
  write.csv(sample_dfs[[i]], out_path, row.names = FALSE)
}

# === STEP 2: Sequence Logo Plotting ===
files <- list.files("trimmed", pattern = "\\.csv$", full.names = TRUE)

for (f in files) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  expanded_seqs <- rep(df$amino_acid, df$count)
  if (length(expanded_seqs) == 0) next
  
  p <- ggseqlogo(expanded_seqs, method = 'prob') +
    theme_bw() +
    theme(
      axis.line = element_line(colour = 'black'),
      axis.ticks = element_line(colour = 'black'),
      axis.title = element_text(colour = 'black'),
      axis.text = element_text(colour = 'black')
    )
  
  out_png <- sub("\\.csv$", ".png", f)
  ggsave(out_png, p, width = 6, height = 3, dpi = 300)
}

# === STEP 3: Amino Acid Class Motif Barplots ===

aa_categories <- list(
  acidic      = c('D','E'),
  basic       = c('K','R','H'),
  hydrophobic = c('A','V','I','L','M','F','W','Y'),
  neutral     = c('G','P'),
  polar       = c('S','T','N','Q','C')
)

category_colors <- c(
  acidic      = "red",
  basic       = "blue",
  hydrophobic = "black",
  neutral     = "violet",
  polar       = "green"
)

# Create a lookup from amino acid to category
aa_category_map <- setNames(rep(names(aa_categories), sapply(aa_categories, length)), unlist(aa_categories))

for (f in files) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  expanded_seqs <- rep(df$amino_acid, df$count)
  if (length(expanded_seqs) == 0) next
  
  n_positions <- nchar(expanded_seqs[1])
  mat <- do.call(rbind, strsplit(expanded_seqs, split = ""))
  seq_df <- as.data.frame(mat)
  colnames(seq_df) <- paste0("Pos", 1:n_positions)
  seq_df$SequenceID <- seq_along(expanded_seqs)
  
  long_df <- seq_df %>%
    pivot_longer(cols = starts_with("Pos"), names_to = "Position", values_to = "AA") %>%
    mutate(Category = aa_category_map[AA]) %>%
    na.omit()
  
  prop_df <- long_df %>%
    group_by(Position, Category) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(Position) %>%
    mutate(freq = count / sum(count)) %>%
    ungroup()
  
  sample_name <- tools::file_path_sans_ext(basename(f))
  out_csv <- file.path("trimmed/MINRUN025_motif", paste0(sample_name, "_motif.csv"))
  write.csv(prop_df, out_csv, row.names = FALSE)
  
  p <- ggplot(prop_df, aes(x = Position, y = freq, fill = Category)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = category_colors) +
    scale_y_continuous(expand = c(0,0), name = "Probability (class)") +
    xlab("Amino Acid Position") +
    ggtitle(paste0(sample_name, " AA Category by Position")) +
    theme_bw() +
    theme(
      axis.line = element_line(colour = 'black'),
      axis.ticks = element_line(colour = 'black'),
      axis.title = element_text(colour = 'black'),
      axis.text = element_text(colour = 'black'),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  out_png <- file.path("trimmed/MINRUN025_motif", paste0(sample_name, "_motif.png"))
  ggsave(out_png, p, width = 6, height = 3, dpi = 300)
}


# === STEP 4: Dendrograms using motif clustering with dendextend ===
library(dendextend)
library(RColorBrewer)
library(dplyr)
library(stringdist)

# Load query peptides mapping
query_map <- read.csv("query_map.csv", stringsAsFactors = FALSE)
query_lookup <- setNames(query_map$Query_Peptide, query_map$Sample)

# List motif files
motif_files <- list.files("trimmed", pattern = "\\.csv$", full.names = TRUE)

# Create output dirs if missing
dir.create("cluster_motif/dendrogram", showWarnings = FALSE, recursive = TRUE)
dir.create("cluster_motif/legend", showWarnings = FALSE, recursive = TRUE)

# AA category code mapping
aa_category_code <- list(
  A = c("D", "E"),        # Acidic
  B = c("K", "R", "H"),  # Basic
  H = c("A", "V", "I", "L", "M", "F", "W", "Y"),  # Hydrophobic
  N = c("G", "P"),        # Neutral
  P = c("S", "T", "N", "Q", "C")  # Polar
)
code_map <- setNames(rep(names(aa_category_code), sapply(aa_category_code, length)), unlist(aa_category_code))

for (f in motif_files) {
  sample_name <- tools::file_path_sans_ext(basename(f))
  if (!sample_name %in% names(query_lookup)) {
    warning("Sample not found in query_map: ", sample_name)
    next
  }
  
  df <- read.csv(f, stringsAsFactors = FALSE)
  if (!"amino_acid" %in% colnames(df)) next
  
  # Convert AA sequences to motif codes
  df$motif <- sapply(df$amino_acid, function(seq) {
    paste0(sapply(strsplit(seq, "")[[1]], function(aa) {
      code <- code_map[[aa]]
      if (is.null(code)) return("X") else return(code)
    }), collapse = "")
  })
  
  # Summarize motif counts
  motif_freq <- df %>%
    group_by(motif) %>%
    summarise(motif_count = sum(count), .groups = "drop") %>%
    arrange(desc(motif_count))
  
  motifs <- motif_freq$motif
  if (length(motifs) < 2) {
    message("Skipping sample ", sample_name, ": fewer than 2 motifs")
    next
  }
  
  # Check motif lengths uniform
  motif_lengths <- nchar(motifs)
  if (length(unique(motif_lengths)) > 1) {
    warning("Skipping sample ", sample_name, ": motifs not all same length")
    next
  }
  
  # Compute Hamming distance and hierarchical clustering
  dist_matrix <- stringdistmatrix(motifs, motifs, method = "hamming")
  hc <- hclust(as.dist(dist_matrix), method = "average")
  dend <- as.dendrogram(hc)
  
  # Cut tree into k clusters
  k <- 4
  clusters <- cutree(hc, k = k)
  
  # Define distinct colors for clusters
  cluster_colors <- brewer.pal(max(k, 3), "Dark2")[1:k]
  
  # Color dendrogram branches by cluster
  dend <- color_branches(dend, k = k, col = cluster_colors)
  
  # Prepare query motif
  query_seq <- query_lookup[[sample_name]]
  query_motif <- paste0(sapply(strsplit(query_seq, "")[[1]], function(aa) {
    code <- code_map[[aa]]
    if (is.null(code)) return("X") else return(code)
  }), collapse = "")
  
  # Scale label sizes by motif frequency
  freq_vec <- setNames(motif_freq$motif_count, motif_freq$motif)
  label_cex <- pmax(1.0, freq_vec[labels(dend)] / max(freq_vec) * 1.5)
  dend <- set(dend, "labels_cex", label_cex)
  
  # Label colors: red for query motif, cluster colors for others
  label_cols <- cluster_colors[clusters[labels(dend)]]
  names(label_cols) <- labels(dend)
  if (query_motif %in% names(label_cols)) {
    label_cols[query_motif] <- "red"
  }
  dend <- set(dend, "labels_col", label_cols)
  
  # Relabel tips by cluster number (e.g. "C1", "C2") to avoid clutter
  motif_to_cluster <- clusters[labels(dend)]
  labels(dend) <- paste0("C", motif_to_cluster[labels(dend)])
  
  # Save dendrogram plot with legend
  png(file.path("cluster_motif/dendrogram", paste0(sample_name, "_motif_dendextend.png")), width = 1000, height = 800)
  par(mar = c(5, 8, 4, 10))
  plot(dend, horiz = TRUE, main = paste0(sample_name, ": Motif Dendrogram"))
  
  # Prepare legend text with top motif per cluster
  cluster_summary <- data.frame(
    Cluster = 1:k,
    Color = cluster_colors,
    Top_Motif = sapply(1:k, function(cl) {
      top <- motif_freq$motif[clusters[motif_freq$motif] == cl][1]
      if (is.na(top)) return("") else return(top)
    })
  )
  
  legend("topright", legend = paste0("C", cluster_summary$Cluster, ": ", cluster_summary$Top_Motif),
         fill = cluster_summary$Color, border = NA, cex = 0.8, box.col = NA)
  dev.off()
  
  # Create export table with cluster info per peptide
  df_export <- df %>%
    left_join(motif_freq, by = "motif") %>%
    mutate(
      Cluster = clusters[motif],
      Cluster_Size = sapply(Cluster, function(cl) sum(clusters == cl)),
      Freq_in_Cluster = motif_count / Cluster_Size
    )
  
  # Save cluster info CSV
  write.csv(df_export, file.path("cluster_motif/legend", paste0(sample_name, "_cluster_table.csv")), row.names = FALSE)
  
  message("Processed sample: ", sample_name)
}

  
  
