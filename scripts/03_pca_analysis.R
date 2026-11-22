# =============================================================================
# Script 03: PCA - Principal Component Analysis
# Author: Niharika Joshi
# Description: Runs PCA using PLINK output and produces publication-quality
#              plots coloured by population and superpopulation.
# Requirements: ggplot2, dplyr, RColorBrewer, cowplot
# =============================================================================

library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(cowplot)

# --- Configuration ---
QC_DIR    <- "data/qc"
DATA_DIR  <- "data"
FIGS_DIR  <- "results/figures"
CHR       <- 22

dir.create(FIGS_DIR, recursive = TRUE, showWarnings = FALSE)

cat("============================================\n")
cat("  PCA Analysis\n")
cat("============================================\n")

# =============================================================================
# Step 1: Run PCA using PLINK (called from R for reproducibility)
# =============================================================================
cat("[1/4] Running PCA with PLINK...\n")

plink_cmd <- sprintf(
  "plink --bfile %s/chr%s_pruned --pca 10 --out %s/chr%s_pca --allow-extra-chr",
  QC_DIR, CHR, QC_DIR, CHR
)
system(plink_cmd)
cat("  PCA complete\n")

# =============================================================================
# Step 2: Load PCA results and population labels
# =============================================================================
cat("[2/4] Loading PCA results...\n")

# Load eigenvectors (PC scores per individual)
eigenvec <- read.table(
  file.path(QC_DIR, paste0("chr", CHR, "_pca.eigenvec")),
  header = FALSE,
  col.names = c("FID", "IID", paste0("PC", 1:10))
)

# Load eigenvalues (variance explained by each PC)
eigenval <- read.table(
  file.path(QC_DIR, paste0("chr", CHR, "_pca.eigenval")),
  header = FALSE, col.names = "eigenvalue"
)

# Calculate % variance explained
pct_var <- round((eigenval$eigenvalue / sum(eigenval$eigenvalue)) * 100, 2)
cat(sprintf("  PC1: %.2f%% | PC2: %.2f%% variance explained\n", pct_var[1], pct_var[2]))

# Load population panel
panel <- read.table(
  file.path(DATA_DIR, "integrated_call_samples_v3.20130502.ALL.panel"),
  header = TRUE
)

# Merge PCA scores with population labels
pca_data <- eigenvec %>%
  left_join(panel %>% select(sample, pop, super_pop),
            by = c("IID" = "sample"))

# =============================================================================
# Step 3: Define colours for superpopulations
# =============================================================================
superpop_colours <- c(
  "AFR" = "#E41A1C",   # Red
  "AMR" = "#FF7F00",   # Orange
  "EAS" = "#4DAF4A",   # Green
  "EUR" = "#377EB8",   # Blue
  "SAS" = "#984EA3"    # Purple
)

# =============================================================================
# Step 4: Plot PCA coloured by superpopulation
# =============================================================================
cat("[3/4] Generating PCA plots...\n")

p1 <- ggplot(pca_data, aes(x = PC1, y = PC2, colour = super_pop)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_colour_manual(
    values = superpop_colours,
    name = "Superpopulation"
  ) +
  labs(
    title = "Global Human Population Structure",
    subtitle = "1000 Genomes Phase 3, chromosome 22 only: 2,504 individuals, 26 populations",
    x = sprintf("PC1 (%.2f%% of top-10 PC variance)", pct_var[1]),
    y = sprintf("PC2 (%.2f%% of top-10 PC variance)", pct_var[2]),
    caption = "Analysis: PLINK 1.9 | Visualisation: R/ggplot2 | Author: Niharika Joshi"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(colour = "grey50", size = 10),
    legend.position = "right",
    legend.title  = element_text(face = "bold"),
    plot.caption  = element_text(colour = "grey60", size = 8)
  )

# PC3 vs PC4 (additional variation)
p2 <- ggplot(pca_data, aes(x = PC3, y = PC4, colour = super_pop)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_colour_manual(values = superpop_colours, name = "Superpopulation") +
  labs(
    title = "PC3 vs PC4",
    x = sprintf("PC3 (%.2f%% of top-10 PC variance)", pct_var[3]),
    y = sprintf("PC4 (%.2f%% of top-10 PC variance)", pct_var[4])
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

# Variance explained scree plot
scree_data <- data.frame(PC = 1:10, pct_var = pct_var)
p3 <- ggplot(scree_data, aes(x = PC, y = pct_var)) +
  geom_col(fill = "#377EB8", alpha = 0.8) +
  geom_line(colour = "grey40") +
  geom_point(size = 2, colour = "grey40") +
  scale_x_continuous(breaks = 1:10) +
  labs(title = "Variance Among Top 10 PCs", x = "Principal Component", y = "% of top-10 PC variance") +
  theme_classic(base_size = 11)

# Combined figure
combined <- plot_grid(
  p1,
  plot_grid(p2, p3, ncol = 2, labels = c("B", "C")),
  nrow = 2, labels = c("A", ""), rel_heights = c(1.6, 1)
)

# Save
ggsave(
  filename = file.path(FIGS_DIR, "01_pca_analysis.png"),
  plot = combined,
  width = 12, height = 10, dpi = 300, bg = "white"
)

cat("  PCA plot saved: results/figures/01_pca_analysis.png\n")

# =============================================================================
# Step 5: Save PCA data for downstream use
# =============================================================================
cat("[4/4] Saving PCA data table...\n")
write.csv(pca_data, file.path(QC_DIR, "pca_results_with_populations.csv"), row.names = FALSE)

cat("\n============================================\n")
cat("  PCA complete. Output:\n")
cat("  results/figures/01_pca_analysis.png\n")
cat("  data/qc/pca_results_with_populations.csv\n")
cat("  Next: bash scripts/04_admixture.sh\n")
cat("============================================\n")
