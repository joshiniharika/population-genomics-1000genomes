# =============================================================================
# Script 05: ROH â Runs of Homozygosity Analysis
# Author: Niharika Joshi
# Description: Detects runs of homozygosity (ROH) using PLINK and analyses
#              their distribution across superpopulations. ROH are genomic
#              regions inherited identically from both parents, reflecting
#              recent shared ancestry and inbreeding.
#
#              This analysis is directly connected to methods used in
#              ancient DNA research (hapROH, kinship inference) â where
#              ROH patterns are used to infer relatedness and pedigree
#              structure in prehistoric populations.
# Requirements: ggplot2, dplyr, tidyr, RColorBrewer, cowplot
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(cowplot)

# --- Configuration ---
QC_DIR   <- "data/qc"
DATA_DIR <- "data"
FIGS_DIR <- "results/figures"
CHR      <- 22

dir.create(FIGS_DIR, recursive = TRUE, showWarnings = FALSE)

cat("============================================\n")
cat("  ROH Analysis\n")
cat("============================================\n")

# =============================================================================
# Step 1: Run PLINK --homozyg to detect ROH
# =============================================================================
cat("[1/4] Detecting ROH with PLINK...\n")

# PLINK homozyg parameters:
# --homozyg-window-snp 50   : sliding window of 50 SNPs
# --homozyg-snp 50          : minimum 50 SNPs in ROH
# --homozyg-kb 500          : minimum ROH length 500 kb
# --homozyg-gap 1000        : maximum gap within ROH 1000 kb
# --homozyg-density 50      : at least 1 SNP per 50 kb
# --homozyg-het 1           : allow max 1 heterozygous call per window

plink_cmd <- sprintf(
  paste(
    "plink --bfile %s/chr%s_qc",
    "--homozyg",
    "--homozyg-window-snp 50",
    "--homozyg-snp 50",
    "--homozyg-kb 500",
    "--homozyg-gap 1000",
    "--homozyg-density 50",
    "--homozyg-het 1",
    "--out %s/chr%s_roh",
    "--allow-extra-chr"
  ),
  QC_DIR, CHR, QC_DIR, CHR
)

system(plink_cmd)
cat("  â ROH detection complete\n")

# =============================================================================
# Step 2: Load ROH results and population panel
# =============================================================================
cat("[2/4] Loading ROH results...\n")

# Load ROH segments per individual
roh <- read.table(
  file.path(QC_DIR, paste0("chr", CHR, "_roh.hom")),
  header = TRUE
)

# Load ROH summary per individual (total ROH length, count)
roh_indiv <- read.table(
  file.path(QC_DIR, paste0("chr", CHR, "_roh.hom.indiv")),
  header = TRUE
)

# Load population labels
panel <- read.table(
  file.path(DATA_DIR, "integrated_call_samples_v3.20130502.ALL.panel"),
  header = TRUE
)

# Merge with population info
roh_indiv <- roh_indiv %>%
  left_join(panel %>% select(sample, pop, super_pop),
            by = c("IID" = "sample"))

# Calculate FROH (fraction of genome in ROH â proxy for inbreeding)
# For chr22 only: chr22 is ~51 Mb
CHR22_LENGTH_KB <- 51000
roh_indiv <- roh_indiv %>%
  mutate(FROH = KB / CHR22_LENGTH_KB)

cat(sprintf("  Loaded %d individuals with ROH data\n", nrow(roh_indiv)))

# =============================================================================
# Step 3: Define colour scheme
# =============================================================================
superpop_colours <- c(
  "AFR" = "#E41A1C",
  "AMR" = "#FF7F00",
  "EAS" = "#4DAF4A",
  "EUR" = "#377EB8",
  "SAS" = "#984EA3"
)

# =============================================================================
# Step 4: Plot ROH distributions
# =============================================================================
cat("[3/4] Generating ROH plots...\n")

# Plot A: Total ROH length per individual by superpopulation (violin + boxplot)
p_violin <- ggplot(roh_indiv, aes(x = super_pop, y = KB / 1000, fill = super_pop)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.12, fill = "white", outlier.size = 0.5, alpha = 0.8) +
  scale_fill_manual(values = superpop_colours, guide = "none") +
  labs(
    title = "Total ROH Length per Individual by Superpopulation",
    subtitle = "1000 Genomes Phase 3 â Chr22",
    x = "Superpopulation",
    y = "Total ROH Length (Mb)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey50", size = 10)
  )

# Plot B: Number of ROH segments per individual
p_count <- ggplot(roh_indiv, aes(x = super_pop, y = NSEG, fill = super_pop)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.12, fill = "white", outlier.size = 0.5, alpha = 0.8) +
  scale_fill_manual(values = superpop_colours, guide = "none") +
  labs(
    title = "Number of ROH Segments per Individual",
    x = "Superpopulation",
    y = "Number of ROH segments"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# Plot C: FROH distribution (inbreeding proxy)
p_froh <- ggplot(roh_indiv, aes(x = FROH, fill = super_pop, colour = super_pop)) +
  geom_density(alpha = 0.4, size = 0.8) +
  scale_fill_manual(values = superpop_colours, name = "Superpopulation") +
  scale_colour_manual(values = superpop_colours, guide = "none") +
  labs(
    title = expression(F[ROH]~"Distribution by Superpopulation"),
    subtitle = "Proportion of chromosome 22 in ROH (proxy for recent inbreeding)",
    x = expression(F[ROH]),
    y = "Density"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey50", size = 9),
    legend.position = "right"
  )

# Plot D: ROH segment length distribution (all segments)
roh_with_pop <- roh %>%
  left_join(panel %>% select(sample, super_pop), by = c("IID" = "sample"))

p_length <- ggplot(roh_with_pop, aes(x = KB / 1000, fill = super_pop, colour = super_pop)) +
  geom_density(alpha = 0.4, size = 0.7) +
  scale_x_log10(labels = scales::comma) +
  scale_fill_manual(values = superpop_colours, name = "Superpopulation") +
  scale_colour_manual(values = superpop_colours, guide = "none") +
  labs(
    title = "ROH Segment Length Distribution",
    subtitle = "Log scale â shorter ROH = older shared ancestry",
    x = "ROH Length (Mb, log scale)",
    y = "Density"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey50", size = 9)
  )

# Combined figure
combined <- plot_grid(
  p_violin, p_count, p_froh, p_length,
  ncol = 2, labels = c("A", "B", "C", "D")
)

ggsave(
  filename = file.path(FIGS_DIR, "02_roh_analysis.png"),
  plot = combined,
  width = 14, height = 11, dpi = 300, bg = "white"
)

cat("  â ROH plot saved: results/figures/02_roh_analysis.png\n")

# =============================================================================
# Step 5: Summary statistics table
# =============================================================================
cat("[4/4] Computing summary statistics...\n")

summary_table <- roh_indiv %>%
  group_by(super_pop) %>%
  summarise(
    n_individuals    = n(),
    mean_roh_mb      = round(mean(KB) / 1000, 3),
    median_roh_mb    = round(median(KB) / 1000, 3),
    mean_n_segments  = round(mean(NSEG), 1),
    mean_froh        = round(mean(FROH), 5),
    .groups = "drop"
  )

print(summary_table)
write.csv(summary_table,
          file.path("results", "roh_summary_by_superpopulation.csv"),
          row.names = FALSE)

cat("\n============================================\n")
cat("  ROH analysis complete. Output:\n")
cat("  â results/figures/02_roh_analysis.png\n")
cat("  â results/roh_summary_by_superpopulation.csv\n")
cat("\n  Connection to ancient DNA:\n")
cat("  ROH detection is the modern-genome equivalent of\n")
cat("  hapROH â which adapts this method for low-coverage\n")
cat("  ancient DNA to infer pedigree and kinship.\n")
cat("============================================\n")
