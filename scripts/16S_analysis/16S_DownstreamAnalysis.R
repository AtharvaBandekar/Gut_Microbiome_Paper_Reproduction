# --- Project: Gut Microbiome Fiber Study Reproduction ---
# --- Analysis: 16S Amplicon Data Downstream Analysis (using authors' processed files) ---

# --- 1. Load Necessary Libraries ---
library(tidyverse) # Data manipulation (dplyr, tibble) and visualization (ggplot2)
library(phyloseq)  # Microbiome data objects and analysis
library(vegan)     # Ecological diversity metrics
library(ape)       # For reading phylogenetic tree files

# --- 2. Load Data Files ---
# ASV table (samples as rows, features/ASVs as columns)
asv_counts_df <- read.csv("data/metadata/asv_table_full.csv", row.names = 1, header = TRUE)

# Fix: Remove 'X' prefix from ASV IDs (column names) if present.
# R prepends 'X' to column names that start with numbers (like hashed ASV IDs).
if (any(startsWith(colnames(asv_counts_df), "X"))) {
  colnames(asv_counts_df) <- sub("^X", "", colnames(asv_counts_df))
}

# Taxonomy table (features/ASVs as rows, taxonomic ranks as columns)
# Loaded as raw for common feature intersection; not fully integrated into basic phyloseq object.
taxonomy_df_raw <- read.csv("data/metadata/taxonomy_table_full.csv", row.names = 1, header = TRUE)

# Fix: Remove 'X' prefix from taxonomy_df row names (if present).
if (any(startsWith(rownames(taxonomy_df_raw), "X"))) {
  rownames(taxonomy_df_raw) <- sub("^X", "", rownames(taxonomy_df_raw))
}
taxonomy_matrix_raw <- as.matrix(taxonomy_df_raw)

# Metadata table (samples as rows, metadata variables as columns)
metadata_df <- read.table("data/metadata/meta_full.tsv", sep = "\t", header = TRUE, row.names = 1, quote = "")

# --- 3. Prepare ASV Table and Metadata for Basic Phyloseq Object ---
# Convert ASV counts to matrix and transpose for phyloseq: features (ASVs) as rows, samples as columns.
asv_counts_matrix_phylo <- as.matrix(t(asv_counts_df))

# Critical check: Ensure sample IDs match between ASV table and Metadata.
if (!all(colnames(asv_counts_matrix_phylo) %in% rownames(metadata_df))) {
  stop("Sample IDs in ASV table are missing in metadata. Please check matching IDs!")
}
if (!all(rownames(metadata_df) %in% colnames(asv_counts_matrix_phylo))) {
  stop("Sample IDs in metadata are missing in ASV table. Please check matching IDs!")
}

# Create basic phyloseq object with only OTU table (ASVs) and Sample Data (metadata).
# Taxonomy is omitted here due to identified feature ID mismatch between author's files for full integration.
OTU <- otu_table(asv_counts_matrix_phylo, taxa_are_rows = TRUE)
SAM <- sample_data(metadata_df)
phyloseq_obj_basic <- phyloseq(OTU, SAM)

message("Phyloseq basic object (ASVs + Sample Data) created.")
print(phyloseq_obj_basic)

# --- 4. Data Normalization (Rarefaction) ---
message("Normalizing data via rarefaction...")

# Calculate sample sums before rarefaction.
sample_sums_before_rarefaction <- sample_sums(phyloseq_obj_basic)
message(paste("Minimum sequence count per sample before rarefaction: ", min(sample_sums_before_rarefaction)))

# Define rarefaction depth based on minimum sample depth (80% of min).
rarefaction_depth <- floor(min(sample_sums_before_rarefaction) * 0.8)
if (rarefaction_depth == 0) {
  rarefaction_depth <- 1 # Ensure minimum is at least 1.
}
message(paste("Rarefying to depth: ", rarefaction_depth))

# Prune samples below rarefaction depth (should be none).
phyloseq_rarefied_basic <- prune_samples(sample_sums(phyloseq_obj_basic) >= rarefaction_depth, phyloseq_obj_basic)
# Perform rarefaction (random subsampling) for normalization.
phyloseq_rarefied_basic <- rarefy_even_depth(phyloseq_rarefied_basic, sample.size = rarefaction_depth, rngseed = 123, replace = FALSE)

message("Phyloseq object after rarefaction:")
print(phyloseq_rarefied_basic)

# --- 5. Alpha Diversity Analysis (Shannon) ---
message("Calculating Shannon alpha diversity...")

# Estimate Shannon diversity index.
alpha_diversity_shannon <- estimate_richness(phyloseq_rarefied_basic, measures = c("Shannon"))
message("Alpha Diversity (Shannon) results:")
print(head(alpha_diversity_shannon))

# Merge alpha diversity results with metadata for plotting.
alpha_diversity_df <- data.frame(sample_data(phyloseq_rarefied_basic), alpha_diversity_shannon)

# Plot Shannon diversity by Diet.
p_alpha_shannon <- ggplot(alpha_diversity_df, aes(x = Diet, y = Shannon, fill = Diet)) +
  geom_boxplot() + # Box plots show distribution.
  geom_point(position = position_jitter(width = 0.2), alpha = 0.6) + # Jittered points for individual data.
  theme_bw() + # Clean theme.
  labs(title = "Alpha Diversity (Shannon) by Diet", y = "Shannon Diversity Index") + # Labels.
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis text.
print(p_alpha_shannon)

message("Alpha diversity plot (Shannon) generated.")

# --- 6. Beta Diversity Analysis (Jaccard) ---
message("Calculating Jaccard beta diversity...")

# Calculate Jaccard distance matrix (non-phylogenetic).
dist_jaccard <- distance(phyloseq_rarefied_basic, method = "jaccard", binary = FALSE)
message("Jaccard Distance Matrix created.")

# Perform Principal Coordinates Analysis (PCoA) for dimensionality reduction.
pcoa_jaccard <- ordinate(phyloseq_rarefied_basic, method = "PCoA", distance = dist_jaccard)

# Plot PCoA colored by Diet.
p_beta_jaccard <- plot_ordination(phyloseq_rarefied_basic, pcoa_jaccard, color = "Diet") +
  geom_point(size = 3, alpha = 0.7) +
  theme_bw() +
  labs(title = "Beta Diversity (Jaccard) PCoA by Diet") +
  stat_ellipse(type = "norm", aes(color = Diet)) # Add ellipses for visual grouping.
print(p_beta_jaccard)

message("Beta diversity PCoA plot (Jaccard) generated.")

# --- 7. Basic Taxonomic Composition (using common features for plotting) ---
message("Attempting basic taxonomic composition visualization for common features...")

# Identify features (ASVs) common to both the ASV table (columns) and Taxonomy table (rows).
# This is crucial for linking counts to taxonomy for visualization.
features_in_asv_cols_and_tax_rows <- intersect(colnames(asv_counts_df), rownames(taxonomy_df_raw))

if (length(features_in_asv_cols_and_tax_rows) > 0) {
  message("Common features found. Generating Phylum level bar plot.")
  
  # Prune ASV table to include only the common features.
  asv_counts_common_features <- asv_counts_df[, features_in_asv_cols_and_tax_rows, drop = FALSE]
  
  # Prune taxonomy table to include only the common features.
  taxonomy_common_features <- taxonomy_df_raw[features_in_asv_cols_and_tax_rows, , drop = FALSE]
  
  # Create a temporary phyloseq object with common features for barplot generation.
  temp_otu <- otu_table(as.matrix(t(asv_counts_common_features)), taxa_are_rows = TRUE)
  temp_tax <- tax_table(as.matrix(taxonomy_common_features))
  temp_sam <- sample_data(metadata_df)
  
  phyloseq_obj_common_features <- phyloseq(temp_otu, temp_tax, temp_sam)
  
  # Aggregate counts to Phylum level.
  phyloseq_phylum <- tax_glom(phyloseq_obj_common_features, taxrank = "Phylum")
  
  # Transform counts to relative abundance for comparison across samples.
  phyloseq_phylum_relabund <- transform_sample_counts(phyloseq_phylum, function(x) x / sum(x))
  
  # Plot Phylum level bar plot.
  p_phylum_barplot <- plot_bar(phyloseq_phylum_relabund, fill = "Phylum", x = "Diet") +
    geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
    labs(x = "", y = "Relative Abundance") +
    facet_wrap(~ Diet, scales = "free_x") + # Create separate panels for each Diet.
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p_phylum_barplot)
  message("Phylum level bar plot generated for common features.")
  
} else {
  message("No common features found between ASV table and Taxonomy table for bar plot.")
  message("Skipping taxonomic composition visualization.")
}

message("\n16S Amplicon analysis complete (with data limitations noted). Ready for saving results!\n")

# --- 8. Save Results (Optional) ---
# Ensure 'results/16S_results/figures/' directory exists before running to save plots.
ggsave("results/16S_results/figures/alpha_diversity_shannon.png", p_alpha_shannon)
ggsave("results/16S_results/figures/beta_diversity_jaccard.png", p_beta_jaccard)
ggsave("results/16S_results/figures/phylum_barplot.png", p_phylum_barplot)
saveRDS(phyloseq_obj_basic, "data/processed/16S_processed/phyloseq_object_basic.rds")
saveRDS(phyloseq_rarefied_basic, "data/processed/16S_processed/phyloseq_object_rarefied_basic.rds")
saveRDS(phyloseq_obj_common_features, "data/processed/16S_processed/phyloseq_object_common_features_for_taxa.rds")