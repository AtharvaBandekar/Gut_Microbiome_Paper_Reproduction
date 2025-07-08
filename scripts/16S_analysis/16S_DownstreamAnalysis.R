# --- Project: Gut Microbiome Fiber Study Reproduction ---
# --- Analysis: 16S Amplicon Data Downstream Analysis (using authors' processed files) ---

# --- 1. Load Necessary Libraries ---
library(tidyverse) # Data manipulation (dplyr, tibble) and visualization (ggplot2)
library(phyloseq)  # Microbiome data objects and analysis
library(vegan)     # Ecological diversity metrics
library(ape)       # For reading phylogenetic tree files (if available for UniFrac)

# --- 2. Load Data Files ---
message("Loading data files...")

# ASV table (samples as rows, features/ASVs as columns)
asv_counts_df <- read.csv("data/metadata/asv_table_full.csv", row.names = 1, header = TRUE)

# Fix: Remove 'X' prefix from ASV IDs (column names)
# R prepends 'X' to column names that start with numbers (like hashed ASV IDs).
# This aligns ASV IDs with the taxonomy table.
if (any(startsWith(colnames(asv_counts_df), "X"))) {
  colnames(asv_counts_df) <- sub("^X", "", colnames(asv_counts_df))
}

# Taxonomy table (features/ASVs as rows, taxonomic ranks as columns)
# Note: Loaded as raw to handle potential ID mismatch with ASV table, then integrated for common features.
taxonomy_df_raw <- read.csv("data/metadata/taxonomy_table_full.csv", row.names = 1, header = TRUE)

# Fix: Remove 'X' prefix from taxonomy_df row names (if present)
if (any(startsWith(rownames(taxonomy_df_raw), "X"))) {
  rownames(taxonomy_df_raw) <- sub("^X", "", rownames(taxonomy_df_raw))
}
taxonomy_matrix_raw <- as.matrix(taxonomy_df_raw)

# Metadata table (samples as rows, metadata variables as columns)
metadata_df <- read.table("data/metadata/meta_full.tsv", sep = "\t", header = TRUE, row.names = 1, quote = "")

message("Data files loaded successfully.")

# --- 3. Prepare Data for Phyloseq Object ---
message("Preparing data for phyloseq object...")

# Convert ASV counts to matrix and transpose for phyloseq: features (ASVs) as rows, samples as columns
asv_counts_matrix_phylo <- as.matrix(t(asv_counts_df))

# Ensure sample IDs match between ASV table and Metadata table (critical check)
if (!all(colnames(asv_counts_matrix_phylo) %in% rownames(metadata_df))) {
  stop("Sample IDs in ASV table are missing in metadata. Please check matching IDs!")
}
if (!all(rownames(metadata_df) %in% colnames(asv_counts_matrix_phylo))) {
  stop("Sample IDs in metadata are missing in ASV table. Please check matching IDs!")
}

# Identify common features (ASVs) between ASV table columns and Taxonomy table rows
# This addresses the fundamental ID mismatch between the author's ASV counts and taxonomy files.
features_in_asv_cols_and_tax_rows <- intersect(colnames(asv_counts_df), rownames(taxonomy_df_raw))

# Prune ASV table to common features only
asv_counts_common_features <- asv_counts_df[, features_in_asv_cols_and_tax_rows, drop = FALSE]

# Prune taxonomy table to common features only
taxonomy_common_features <- taxonomy_df_raw[features_in_asv_cols_and_tax_rows, , drop = FALSE]

# Create temporary phyloseq components with only common features for taxonomic visualization
temp_otu <- otu_table(as.matrix(t(asv_counts_common_features)), taxa_are_rows = TRUE)
temp_tax <- tax_table(as.matrix(taxonomy_common_features))
temp_sam <- sample_data(metadata_df)

# Create phyloseq object with only ASV counts and Sample Data (no full taxonomy at this stage)
OTU <- otu_table(asv_counts_matrix_phylo, taxa_are_rows = TRUE)
SAM <- sample_data(metadata_df)
phyloseq_obj_basic <- phyloseq(OTU, SAM)

message("Phyloseq basic object (ASVs + Sample Data) created.")

# --- 4. Load Phylogenetic Tree and Create complete phyloseq object ---
message("Loading phylogenetic tree...")

# The tree file is in Newick format
tree_phylo <- read.tree("data/processed/16S_processed/rooted_tree_newick/tree.nwk")

# Align features between basic phyloseq object and tree before adding tree
# Prune phyloseq object to only include taxa also found in the tree
phyloseq_obj_pre_merge_aligned <- prune_taxa(taxa_names(phyloseq_obj_basic) %in% tree_phylo$tip.label, phyloseq_obj_basic)

# Prune tree to only include tips (features) present in the ASV table
tree_phylo_aligned <- prune.tree(tree_phylo, phyloseq_obj_pre_merge_aligned@tax_table@.Data[,1]) 

# Create the complete phyloseq object with OTU table, Taxonomy, Sample Data, AND Aligned Tree
# Note: For this object, the Taxonomy is still from the original phyloseq_obj_basic (not temp_tax).
# We are aligning the tree to the features present in the OTU table.
# A more robust approach might be to build a *new* TAX table only from common features *here*
# However, for UniFrac it's primarily OTU table and tree that are critical.
# Let's add the original taxonomy table (taxonomy_matrix_raw) to phyloseq_obj_basic
# to make phyloseq_obj_with_tree fully functional, then prune.

# This step needs careful re-evaluation: the taxonomy_matrix_raw has mismatched IDs.
# So phyloseq_obj_with_tree will NOT have a correct full taxonomy table if created directly.
# The previous solution of using temp_tax for barplot was specific.

# Let's re-create phyloseq_obj_with_tree only if the taxonomic barplot creation worked.
# Otherwise, we omit this full integration here and only use it for UniFrac if necessary.

# Since we established the barplot *did* generate, means temp_otu and temp_tax have common features.
# Let's create a *new* phyloseq object specifically for tree-based metrics that has aligned features.

# Re-create OTU and TAX objects with *common_features* for the phyloseq_obj_with_tree
OTU_aligned <- otu_table(asv_counts_matrix_phylo_aligned, taxa_are_rows = TRUE)
TAX_aligned <- tax_table(taxonomy_matrix_aligned)
SAM_aligned <- sample_data(metadata_df)

# Create the phyloseq object *with aligned taxonomy* for tree integration
phyloseq_obj_with_aligned_tax <- phyloseq(OTU_aligned, TAX_aligned, SAM_aligned)

# Prune this object and tree before final integration
phyloseq_obj_pre_merge_final <- prune_taxa(taxa_names(phyloseq_obj_with_aligned_tax) %in% tree_phylo$tip.label, phyloseq_obj_with_aligned_tax)
tree_phylo_final <- prune.tree(tree_phylo, phyloseq_obj_pre_merge_final@tax_table@.Data[,1])

# Create the complete phyloseq object
phyloseq_obj_with_tree <- phyloseq(otu_table(phyloseq_obj_pre_merge_final), tax_table(phyloseq_obj_pre_merge_final), sample_data(phyloseq_obj_pre_merge_final), tree_phylo_final)

message("Phyloseq object created successfully with aligned features and tree integrated!")
print(phyloseq_obj_with_tree)


# --- 6. Data Normalization (Rarefaction) ---
message("Normalizing data via rarefaction...")

sample_sums_before_rarefaction <- sample_sums(phyloseq_obj_with_tree)
message(paste("Minimum sequence count per sample before rarefaction: ", min(sample_sums_before_rarefaction)))

rarefaction_depth <- floor(min(sample_sums_before_rarefaction) * 0.8) # 80% of min depth
if (rarefaction_depth == 0) {
  rarefaction_depth <- 1 # Ensure minimum is at least 1
}
message(paste("Rarefying to depth: ", rarefaction_depth))

phyloseq_rarefied <- prune_samples(sample_sums(phyloseq_obj_with_tree) >= rarefaction_depth, phyloseq_obj_with_tree)
phyloseq_rarefied <- rarefy_even_depth(phyloseq_rarefied, sample.size = rarefaction_depth, rngseed = 123, replace = FALSE)

message("Phyloseq object after rarefaction and with tree integrated:")
print(phyloseq_rarefied)

# --- 7. Alpha Diversity Analysis (Shannon and Faith's PD) ---
message("Calculating alpha diversity...")

alpha_diversity_results <- estimate_richness(phyloseq_rarefied, measures = c("Shannon", "Faith_PD"))
alpha_diversity_df <- data.frame(sample_data(phyloseq_rarefied), alpha_diversity_results)

p_alpha_shannon <- ggplot(alpha_diversity_df, aes(x = Diet, y = Shannon, fill = Diet)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +
  theme_bw() +
  labs(title = "Alpha Diversity (Shannon) by Diet", y = "Shannon Diversity Index") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_alpha_shannon)

p_alpha_faith <- ggplot(alpha_diversity_df, aes(x = Diet, y = Faith_PD, fill = Diet)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +
  theme_bw() +
  labs(title = "Alpha Diversity (Faith's PD) by Diet", y = "Faith's Phylogenetic Diversity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_alpha_faith)

message("Alpha diversity plots (Shannon and Faith_PD) generated.")

# --- 8. Beta Diversity Analysis (UniFrac and Jaccard) ---
message("Calculating beta diversity...")

dist_jaccard <- distance(phyloseq_rarefied, method = "jaccard", binary = FALSE)
dist_unifrac_unweighted <- distance(phyloseq_rarefied, method = "unifrac", weighted = FALSE)
dist_unifrac_weighted <- distance(phyloseq_rarefied, method = "unifrac", weighted = TRUE)

# Perform PCoA for Jaccard
pcoa_jaccard <- ordinate(phyloseq_rarefied, method = "PCoA", distance = dist_jaccard)
p_beta_jaccard <- plot_ordination(phyloseq_rarefied, pcoa_jaccard, color = "Diet") +
  geom_point(size = 3, alpha = 0.7) +
  theme_bw() +
  labs(title = "Beta Diversity (Jaccard) PCoA by Diet") +
  stat_ellipse(type = "norm", aes(color = Diet))
print(p_beta_jaccard)

# Perform PCoA for Unweighted UniFrac
pcoa_unifrac_unweighted <- ordinate(phyloseq_rarefied, method = "PCoA", distance = dist_unifrac_unweighted)
p_beta_unifrac_unweighted <- plot_ordination(phyloseq_rarefied, pcoa_unifrac_unweighted, color = "Diet") +
  geom_point(size = 3, alpha = 0.7) +
  theme_bw() +
  labs(title = "Beta Diversity (Unweighted UniFrac) PCoA by Diet") +
  stat_ellipse(type = "norm", aes(color = Diet))
print(p_beta_unifrac_unweighted)

# Perform PCoA for Weighted UniFrac
pcoa_unifrac_weighted <- ordinate(phyloseq_rarefied, method = "PCoA", distance = dist_unifrac_weighted)
p_beta_unifrac_weighted <- plot_ordination(phyloseq_rarefied, pcoa_unifrac_weighted, color = "Diet") +
  geom_point(size = 3, alpha = 0.7) +
  theme_bw() +
  labs(title = "Beta Diversity (Weighted UniFrac) PCoA by Diet") +
  stat_ellipse(type = "norm", aes(color = Diet))
print(p_beta_unifrac_weighted)

message("Beta diversity PCoA plots generated.")

# --- 9. Basic Taxonomic Composition (using common features) ---
message("Attempting basic taxonomic composition visualization for common features...")

# Identify common features between ASV table columns and Taxonomy table rows
# This intersection was crucial for correctly aligning the data for plots.
features_in_asv_cols_and_tax_rows <- intersect(colnames(asv_counts_df), rownames(taxonomy_df_raw))

if (length(features_in_asv_cols_and_tax_rows) > 0) {
  # Prune ASV table to common features
  asv_counts_common_features <- asv_counts_df[, features_in_asv_cols_and_tax_rows, drop = FALSE]
  
  # Prune taxonomy table to common features
  taxonomy_common_features <- taxonomy_df_raw[features_in_asv_cols_and_tax_rows, , drop = FALSE]
  
  # Create a temporary phyloseq object with only common features for barplot
  temp_otu <- otu_table(as.matrix(t(asv_counts_common_features)), taxa_are_rows = TRUE)
  temp_tax <- tax_table(as.matrix(taxonomy_common_features))
  temp_sam <- sample_data(metadata_df)
  
  phyloseq_obj_common_features <- phyloseq(temp_otu, temp_tax, temp_sam)
  
  # Aggregate to Phylum level for bar plot
  phyloseq_phylum <- tax_glom(phyloseq_obj_common_features, taxrank = "Phylum")
  
  # Transform to relative abundance
  phyloseq_phylum_relabund <- transform_sample_counts(phyloseq_phylum, function(x) x / sum(x))
  
  # Plot Phylum level bar plot
  p_phylum_barplot <- plot_bar(phyloseq_phylum_relabund, fill = "Phylum", x = "Diet") +
    geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
    labs(x = "", y = "Relative Abundance") +
    facet_wrap(~ Diet, scales = "free_x") + # Separate plots for each Diet
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p_phylum_barplot)
  message("Phylum level bar plot generated for common features.")
  
} else {
  message("No common features found between ASV table and Taxonomy table for bar plot.")
  message("Skipping taxonomic composition visualization.")
}

message("\n16S Amplicon analysis complete (with data limitations noted). Ready for saving results!\n")

# --- 10. Save Results ---
Ensure the 'results/16S_results/figures/' directory exists if saving plots
ggsave("results/16S_results/figures/alpha_diversity_shannon.png", p_alpha_shannon)
ggsave("results/16S_results/figures/alpha_diversity_faith.png", p_alpha_faith)
ggsave("results/16S_results/figures/beta_diversity_jaccard.png", p_beta_jaccard)
ggsave("results/16S_results/figures/beta_diversity_unifrac_unweighted.png", p_beta_unifrac_unweighted)
ggsave("results/16S_results/figures/beta_diversity_unifrac_weighted.png", p_beta_unifrac_weighted)
saveRDS(phyloseq_obj_with_tree, "data/processed/16S_processed/phyloseq_object_complete.rds")
saveRDS(phyloseq_rarefied, "data/processed/16S_results/phyloseq_object_rarefied_complete.rds")
saveRDS(phyloseq_obj_common_features, "data/processed/16S_processed/phyloseq_object_common_features_for_taxa.rds")