# --- Project: Gut Microbiome Fiber Study Reproduction ---
# --- Analysis: Metabolomics Data Analysis ---

# --- 1. Load Necessary Libraries ---
library(tidyverse) # For data manipulation (dplyr, tibble) and visualization (ggplot2)

# --- 2. Load Metabolomics Data Files ---
message("Loading metabolomics data files...")

# Load main metadata (prepared from meta_full.tsv)
metadata_full_df <- read.table("data/metadata/meta_full.tsv", sep = "\t", header = TRUE, row.names = 1, quote = "")

# Load SCFA data (from meta_scfa.csv)
scfa_df <- read.csv("data/metadata/meta_scfa.csv", header = TRUE, row.names = 1)

# Load subset metadata (from meta_subset.csv)
meta_subset_df <- read.csv("data/metadata/meta_subset.csv", header = TRUE, row.names = 1)

message("Metabolomics data files loaded successfully.")

# --- 3. Initial Data Inspection ---
# Print dimensions and head of loaded dataframes for verification.
message("\nDimensions and head of metadata_full_df:")
print(dim(metadata_full_df))
print(head(metadata_full_df))

message("\nDimensions, head, and column names of scfa_df:")
print(dim(scfa_df))
print(head(scfa_df))
print(colnames(scfa_df))

message("\nDimensions, head, and column names of meta_subset_df:")
print(dim(meta_subset_df))
print(head(meta_subset_df))
print(colnames(meta_subset_df))

# --- 4. Prepare and Plot SCFA Data (Figure S9 Reproduction) ---
message("Preparing and plotting SCFA data (Figure S9 reproduction)...")

# Ensure 'Diet' column is a factor for consistent plotting order.
meta_scfa_df$Diet <- factor(meta_scfa_df$Diet, levels = c("Cellulose", "Inulin", "Pectin", "Assorted Fiber"))

# Reshape SCFA data to long format for ggplot2.
scfa_long_df <- meta_scfa_df %>%
  select(Diet, Acetate, Propionate, Butyrate, Total.SCFA) %>%
  pivot_longer(
    cols = c(Acetate, Propionate, Butyrate, Total.SCFA),
    names_to = "SCFA_Type",
    values_to = "Concentration"
  )

# Generate boxplots for each SCFA type by Diet.
p_scfa <- ggplot(scfa_long_df, aes(x = Diet, y = Concentration, fill = Diet)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.7) +
  facet_wrap(~ SCFA_Type, scales = "free_y") + # Separate plots by SCFA type with independent y-axes.
  theme_bw() +
  labs(title = "Cecal Short-Chain Fatty Acids by Diet (Figure S9 Reproduction)", y = "Concentration (µmoles/g wet weight)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_scfa)
message("SCFA plot (Figure S9 reproduction) generated.")

# --- 5. Prepare and Plot Host Phenotypes (Figure 2 Reproduction) ---
message("Preparing and plotting Host Phenotypes (Figure 2 reproduction)...")

# Select and reshape host phenotype data (Epididymal Fat Pad, Liver Triglycerides, Glucose).
host_phenotypes_df <- metadata_full_df %>%
  select(Diet, GDATdivBW, LiverTriglycerides, Glucose) %>%
  pivot_longer(
    cols = c(GDATdivBW, LiverTriglycerides, Glucose),
    names_to = "Phenotype_Type",
    values_to = "Value"
  )

# Rename phenotype types for clearer plot labels.
host_phenotypes_df$Phenotype_Type <- factor(host_phenotypes_df$Phenotype_Type,
                                            levels = c("GDATdivBW", "LiverTriglycerides", "Glucose"),
                                            labels = c("Epididymal Fat Pad (% of BW)", "Liver Triglycerides (nmol/g)", "Glucose (arbitrary units)"))

# Generate boxplots for each host phenotype by Diet.
p_phenotypes <- ggplot(host_phenotypes_df, aes(x = Diet, y = Value, fill = Diet)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.7) +
  facet_wrap(~ Phenotype_Type, scales = "free_y") + # Separate plots by phenotype type with independent y-axes.
  theme_bw() +
  labs(title = "Host Metabolic Phenotypes by Diet (Figure 2 Reproduction)", y = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_phenotypes)
message("Host phenotypes plot (Figure 2 reproduction) generated.")

# --- 6. Save Results (Optional) ---
# Ensure the 'results/metabolomics_results/figures/' directory exists before saving plots.
ggsave("results/metabolomics_results/figures/scfa_boxplot.png", p_scfa, width = 12, height = 8, units = "in", dpi = 300)
ggsave("results/metabolomics_results/figures/host_phenotypes_boxplot.png", p_phenotypes, width = 12, height = 8, units = "in", dpi = 300)

message("\nMetabolomics analysis complete. Results saved and ready for documentation!")