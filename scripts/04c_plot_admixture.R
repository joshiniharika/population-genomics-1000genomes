# =============================================================================
# Script 04c: Plot fastmixture ancestry proportions and CV errors
# Author: Niharika Joshi
# =============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

QC_DIR   <- "data/qc"
ADM_DIR  <- "results/admixture"
FIGS_DIR <- "results/figures"
CHR      <- 22
SEED     <- 1
K_VALUES <- 2:7
dir.create(FIGS_DIR, recursive = TRUE, showWarnings = FALSE)

# Sample order = order in the .fam file (same as the Q files)
fam <- read.table(file.path(QC_DIR, paste0("chr", CHR, "_pruned.fam")),
                  col.names = c("FID", "IID", "PAT", "MAT", "SEX", "PHENO"))
panel <- read.table("data/integrated_call_samples_v3.20130502.ALL.panel",
                    header = TRUE, sep = "\t", fill = TRUE)
meta <- fam %>%
  select(IID) %>%
  left_join(panel %>% select(sample, pop, super_pop), by = c("IID" = "sample")) %>%
  mutate(idx = row_number())

# Order individuals by superpopulation, then population
meta <- meta %>% arrange(super_pop, pop, idx) %>% mutate(order = row_number())

# Read Q matrices into one long table
read_q <- function(K) {
  f <- file.path(ADM_DIR, sprintf("chr%d.K%d.s%d.Q", CHR, K, SEED))
  q <- read.table(f)
  colnames(q) <- paste0("C", seq_len(K))
  q$idx <- seq_len(nrow(q))
  q %>%
    pivot_longer(starts_with("C"), names_to = "component", values_to = "prop") %>%
    mutate(K = paste("K =", K))
}
long <- bind_rows(lapply(K_VALUES, read_q)) %>%
  left_join(meta, by = "idx")

cols <- brewer.pal(8, "Set1")[1:7]

p <- ggplot(long, aes(x = order, y = prop, fill = component)) +
  geom_col(width = 1) +
  facet_grid(K ~ ., switch = "y") +
  scale_fill_manual(values = cols, guide = "none") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Ancestry proportions (fastmixture)",
       subtitle = "1000 Genomes Phase 3, chromosome 22 only; individuals ordered by superpopulation",
       x = "Individuals (AFR, AMR, EAS, EUR, SAS)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        strip.text.y.left = element_text(angle = 0),
        panel.spacing = unit(0.2, "lines"))

ggsave(file.path(FIGS_DIR, "03_admixture_barplot.png"),
       p, width = 12, height = 8, dpi = 300, bg = "white")

# Cross-validation error plot
cv <- read.csv(file.path(ADM_DIR, "cv_errors.csv"))
p_cv <- ggplot(cv, aes(x = K, y = cv_error)) +
  geom_point(alpha = 0.6) +
  stat_summary(fun = mean, geom = "line") +
  scale_x_continuous(breaks = unique(cv$K)) +
  labs(title = "Cross-validation error by K", y = "CV error", x = "K") +
  theme_classic(base_size = 12)
ggsave(file.path(FIGS_DIR, "04_admixture_cv.png"),
       p_cv, width = 6, height = 4, dpi = 300, bg = "white")

cat("Saved figures to", FIGS_DIR, "\n")