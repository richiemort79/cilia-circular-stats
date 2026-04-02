# ============================================================
#  Primary cilia angle analysis — circular statistics
#  Two conditions (Control, ES) x four timepoints (T0, T4, T8, T12)
#  Two datasets: 100 mV/mm and 25 mV/mm electrical stimulation
#
#  Run from the project root:
#    source("R/cilia_circular_analysis.R")
#
#  Output:
#    - Console: sample sizes, Rayleigh tests, U2 results, summary stats
#    - figures/rose_plots_100mVmm.pdf
#    - figures/rose_plots_25mVmm.pdf
#    - figures/rose_plots_combined.pdf
# ============================================================

library(circular)

# ------------------------------------------------------------
# 1. Helper functions
# ------------------------------------------------------------

to_circ <- function(x) {
  x     <- na.omit(x)
  x_360 <- ifelse(x < 0, x + 360, x)
  circular(x_360, units = "degrees", modulo = "2pi")
}

# Permutation-based Watson U2 test.
# Shuffles group labels n_perm times to build a null distribution,
# then computes an exact p-value as the proportion of permuted
# statistics >= the observed statistic.
# The +1 correction avoids p = 0 and keeps the estimate unbiased.
circular_perm_test <- function(x, y, n_perm = 9999) {
  observed_U2 <- watson.two.test(x, y)$statistic
  combined    <- c(as.numeric(x), as.numeric(y))
  nx          <- length(x)
  n_total     <- length(combined)

  perm_U2 <- replicate(n_perm, {
    idx    <- sample(n_total, nx)
    x_perm <- circular(combined[idx],  units = "degrees", modulo = "2pi")
    y_perm <- circular(combined[-idx], units = "degrees", modulo = "2pi")
    watson.two.test(x_perm, y_perm)$statistic
  })

  p_val <- (sum(perm_U2 >= observed_U2) + 1) / (n_perm + 1)
  list(statistic = observed_U2, p.value = p_val, n_perm = n_perm)
}

# Load and parse a dataset CSV into a named list of circular vectors
load_dataset <- function(path) {
  raw           <- read.csv(path, header = TRUE)
  colnames(raw) <- c("Control_T0", "Control_T4", "Control_T8", "Control_T12",
                     "ES_T0",      "ES_T4",      "ES_T8",      "ES_T12")
  raw           <- raw[-1, ]
  raw[]         <- lapply(raw, function(x) as.numeric(as.character(x)))

  timepoints <- c("T0", "T4", "T8", "T12")
  conditions <- c("Control", "ES")

  angles <- list()
  for (cond in conditions) {
    angles[[cond]] <- list()
    for (tp in timepoints) {
      angles[[cond]][[tp]] <- to_circ(raw[[paste0(cond, "_", tp)]])
    }
  }
  angles
}

# Run full analysis for one dataset, print results, return summary_df
run_analysis <- function(angles, dataset_label, n_perm = 9999) {
  timepoints <- c("T0", "T4", "T8", "T12")
  conditions <- c("Control", "ES")

  cat("\n", strrep("=", 60), "\n", sep = "")
  cat(" Dataset:", dataset_label, "\n")
  cat(strrep("=", 60), "\n", sep = "")

  # Sample sizes
  cat("\nSample sizes:\n")
  for (cond in conditions)
    for (tp in timepoints)
      cat(sprintf("  %s %s : n = %d\n", cond, tp,
                  length(angles[[cond]][[tp]])))

  # Rayleigh tests
  cat("\n--- Rayleigh tests (non-uniformity) ---\n")
  cat("Control groups are expected to be uniform (no preferred direction).\n")
  cat("Watson-Williams is not appropriate; permutation Watson U2 is used.\n\n")
  for (cond in conditions) {
    for (tp in timepoints) {
      p <- rayleigh.test(angles[[cond]][[tp]])$p.value
      cat(sprintf("  %s %s: p = %.4f %s\n", cond, tp, p,
                  ifelse(p < 0.05, "(concentrated)", "(uniform)")))
    }
  }

  # Permutation Watson U2
  cat("\n--- Permutation Watson U2 tests (Control vs ES) ---\n")
  cat(sprintf("Running %d permutations per timepoint...\n\n", n_perm))

  perm_results <- data.frame(timepoint = character(),
                              U2        = numeric(),
                              p_value   = numeric(),
                              stringsAsFactors = FALSE)

  for (tp in timepoints) {
    result <- circular_perm_test(angles[["Control"]][[tp]],
                                  angles[["ES"]][[tp]],
                                  n_perm = n_perm)
    cat(sprintf("  %s: U2 = %.4f, p = %.4f (%d permutations)\n",
                tp, result$statistic, result$p.value, result$n_perm))
    perm_results <- rbind(perm_results,
                           data.frame(timepoint = tp,
                                      U2        = round(result$statistic, 4),
                                      p_value   = result$p.value,
                                      stringsAsFactors = FALSE))
  }

  perm_results$p_bonferroni <- p.adjust(perm_results$p_value, method = "bonferroni")
  perm_results$p_fdr        <- p.adjust(perm_results$p_value, method = "BH")

  cat("\n--- Results with multiple comparison correction ---\n")
  print(perm_results, row.names = FALSE)
  cat(sprintf("\nNote: p = %.4f is the minimum resolvable value with %d permutations.\n",
              1 / (n_perm + 1), n_perm))

  # Summary statistics
  cat("\n--- Summary statistics ---\n")
  summary_df <- data.frame()
  for (cond in conditions) {
    for (tp in timepoints) {
      a <- angles[[cond]][[tp]]
      summary_df <- rbind(summary_df,
                           data.frame(condition = cond,
                                      timepoint = tp,
                                      n         = length(a),
                                      mean_dir  = round(as.numeric(mean(a)), 2),
                                      rho       = round(rho.circular(a), 3),
                                      circ_sd   = round(as.numeric(sd.circular(a)), 2)))
    }
  }
  print(summary_df, row.names = FALSE)
  cat("\nrho (resultant length): 0 = uniform, 1 = all angles identical.\n")

  invisible(list(perm_results = perm_results, summary_df = summary_df))
}

# Rose plot function for one dataset
make_rose_plots <- function(angles, dataset_label, outfile) {
  timepoints <- c("T0", "T4", "T8", "T12")

  pdf(outfile, width = 14, height = 7)
  par(mfrow = c(2, 4),
      mar   = c(1, 1, 2.5, 1),
      oma   = c(0, 3, 3, 0))

  for (tp in timepoints) {
    plot(angles[["Control"]][[tp]],
         stack = TRUE, bins = 24,
         col = "steelblue", border = "white",
         main = tp, shrink = 1.5)
    arrows.circular(mean(angles[["Control"]][[tp]]),
                    col = "steelblue", lwd = 2)
    mtext(bquote(rho == .(round(rho.circular(angles[["Control"]][[tp]]), 3))),
          side = 1, line = -1, cex = 0.8)
  }

  for (tp in timepoints) {
    plot(angles[["ES"]][[tp]],
         stack = TRUE, bins = 24,
         col = "coral", border = "white",
         main = tp, shrink = 1.5)
    arrows.circular(mean(angles[["ES"]][[tp]]),
                    col = "coral", lwd = 2)
    mtext(bquote(rho == .(round(rho.circular(angles[["ES"]][[tp]]), 3))),
          side = 1, line = -1, cex = 0.8)
  }

  mtext("Control", side = 2, outer = TRUE, line = 1.5,
        at = 0.75, cex = 1.1, font = 2, col = "steelblue")
  mtext("ES",      side = 2, outer = TRUE, line = 1.5,
        at = 0.25, cex = 1.1, font = 2, col = "coral")
  mtext(paste("Primary cilia orientation —", dataset_label),
        side = 3, outer = TRUE, line = 1, cex = 1.2, font = 2)

  dev.off()
  cat("Saved:", outfile, "\n")
}

# Combined rose plot — all datasets as rows
make_combined_plots <- function(all_angles, all_labels, outfile) {
  timepoints  <- c("T0", "T4", "T8", "T12")
  n_datasets  <- length(all_angles)
  n_rows      <- n_datasets * 2   # Control + ES per dataset
  cond_cols   <- c("Control" = "steelblue", "ES" = "coral")

  pdf(outfile, width = 14, height = 3.5 * n_datasets)
  par(mfrow = c(n_rows, 4),
      mar   = c(1, 1, 2, 1),
      oma   = c(0, 5, 3, 0))

  for (i in seq_along(all_angles)) {
    for (cond in c("Control", "ES")) {
      for (tp in timepoints) {
        plot(all_angles[[i]][[cond]][[tp]],
             stack = TRUE, bins = 24,
             col = cond_cols[[cond]], border = "white",
             main = if (cond == "Control") tp else "",
             shrink = 1.5)
        arrows.circular(mean(all_angles[[i]][[cond]][[tp]]),
                        col = cond_cols[[cond]], lwd = 2)
        mtext(bquote(rho == .(round(rho.circular(all_angles[[i]][[cond]][[tp]]), 3))),
              side = 1, line = -1, cex = 0.75)
      }
      # Row label: "100 mV/mm — Control" etc.
      mtext(paste0(all_labels[[i]], "\n", cond),
            side = 2, outer = FALSE, line = -1,
            at = 0.5, cex = 0.8, font = 2,
            col = cond_cols[[cond]],
            las = 3)
    }
  }

  mtext("Primary cilia orientation — Control vs electrical stimulation",
        side = 3, outer = TRUE, line = 1, cex = 1.2, font = 2)

  dev.off()
  cat("Saved:", outfile, "\n")
}

# ------------------------------------------------------------
# 2. Load datasets
# ------------------------------------------------------------

datasets <- list(
  list(
    path  = "data/100mvmm_Cilia_angles_25_9_24.csv",
    label = "100 mV/mm",
    file  = "figures/rose_plots_100mVmm.pdf"
  ),
  list(
    path  = "data/25mvmm_Cilia_angles_30_9_24.csv",
    label = "25 mV/mm",
    file  = "figures/rose_plots_25mVmm.pdf"
  )
)

dir.create("figures", showWarnings = FALSE)

# ------------------------------------------------------------
# 3. Run analysis and generate plots for each dataset
# ------------------------------------------------------------

all_angles <- list()
all_labels <- list()

for (ds in datasets) {
  angles <- load_dataset(ds$path)
  run_analysis(angles, ds$label)
  make_rose_plots(angles, ds$label, ds$file)
  all_angles <- c(all_angles, list(angles))
  all_labels <- c(all_labels, list(ds$label))
}

# ------------------------------------------------------------
# 4. Combined figure
# ------------------------------------------------------------

make_combined_plots(all_angles, all_labels,
                    "figures/rose_plots_combined.pdf")

cat("\nAll done. Figures written to figures/\n")
