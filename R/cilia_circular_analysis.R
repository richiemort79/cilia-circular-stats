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

# ------------------------------------------------------------
# SET THIS FLAG:
#   TRUE  = regenerate plots only (fast, no statistics)
#   FALSE = run full analysis with permutation tests (slow)
# ------------------------------------------------------------
PLOT_ONLY <- FALSE

library(circular)
library(ggplot2)
library(gridExtra)

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

# Run full analysis for one dataset, print results, return results list
run_analysis <- function(angles, dataset_label, n_perm = 9999) {
  timepoints <- c("T0", "T4", "T8", "T12")
  conditions <- c("Control", "ES")

  cat("\n", strrep("=", 60), "\n", sep = "")
  cat(" Dataset:", dataset_label, "\n")
  cat(strrep("=", 60), "\n", sep = "")

  cat("\nSample sizes:\n")
  for (cond in conditions)
    for (tp in timepoints)
      cat(sprintf("  %s %s : n = %d\n", cond, tp,
                  length(angles[[cond]][[tp]])))

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

# ------------------------------------------------------------
# 2. ggplot rose plot for a single group
# ------------------------------------------------------------

make_gg_rose <- function(angles_vec, title = "", subtitle = "", caption = "",
                          col_fill = "#ddddff", col_border = "blue") {

  # Histogram on 0-360 in 20-degree bins
  breaks <- seq(0, 360, by = 20)
  raw    <- as.numeric(angles_vec)
  h      <- hist(raw, breaks = breaks, plot = FALSE)

  d <- data.frame(Angle     = h$mids,
                  Frequency = as.numeric(h$counts))

  # Mean direction arrow — length scaled by rho relative to tallest bar
  circ_vec  <- circular(raw, units = "degrees", modulo = "2pi")
  mean_dir  <- as.numeric(mean(circ_vec))
  rho_val   <- rho.circular(circ_vec)
  arrow_len <- max(d$Frequency) * rho_val

  ggplot(d, aes(x = Angle, y = Frequency)) +
    ggtitle(label = title, subtitle = subtitle) +
    labs(caption = caption) +
    coord_polar(theta = "x", start = pi / 2) +
    geom_bar(stat      = "identity",
             fill      = col_fill,
             color     = col_border,
             linewidth = 0.25) +
    geom_segment(aes(x = mean_dir, xend = mean_dir,
                     y = 0,        yend = arrow_len),
                 color       = col_border,
                 linewidth   = 1,
                 arrow       = arrow(length = unit(0.2, "cm"), type = "closed"),
                 inherit.aes = FALSE) +
    scale_x_continuous(breaks = c(0, 45, 90, 135, 180, 225, 270, 315),
                       expand = c(0.002, 0),
                       limits = c(0, 360)) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title       = element_text(size = 10, face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(size = 8,  hjust = 0.5, color = "grey40"),
      axis.title       = element_blank(),
      axis.text.y      = element_blank(),
      axis.text.x      = element_text(size = 7, color = "grey40"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.margin      = margin(4, 4, 4, 4),
      plot.caption     = element_text(size = 7, hjust = 0.5, color = "grey40")
    )
}

# Build a 2 x 4 page (Control row + ES row) and save to PDF
make_rose_plots <- function(angles, perm_results, dataset_label, outfile) {
  timepoints <- c("T0", "T4", "T8", "T12")
  plot_list  <- list()

  # Control row
  for (tp in timepoints) {
    n_cilia   <- length(angles[["Control"]][[tp]])
    rho_val   <- round(rho.circular(angles[["Control"]][[tp]]), 3)
    rayleigh_p <- round(rayleigh.test(angles[["Control"]][[tp]])$p.value, 3)
    plot_list <- c(plot_list, list(
      make_gg_rose(angles[["Control"]][[tp]],
                   title    = paste("Control -", tp),
                   subtitle = paste0("n = ", n_cilia, ",  rho = ", rho_val,
                                     ",  Rayleigh p = ", format.pval(rayleigh_p, digits = 2)),
                   col_fill   = "#ddeeff",
                   col_border = "steelblue")
    ))
  }

  # ES row
  for (tp in timepoints) {
    n_cilia    <- length(angles[["ES"]][[tp]])
    rho_val    <- round(rho.circular(angles[["ES"]][[tp]]), 3)
    rayleigh_p <- round(rayleigh.test(angles[["ES"]][[tp]])$p.value, 3)
    pr         <- perm_results[perm_results$timepoint == tp, ]
    sub <- paste0("n = ", n_cilia, ",  rho = ", rho_val,
                  ",  Rayleigh p = ", format.pval(rayleigh_p, digits = 2))
    cap <- if (!is.na(pr$U2)) paste0("Watson U2 p (Bonferroni) = ", format.pval(pr$p_bonferroni, digits = 2)) else ""
    plot_list <- c(plot_list, list(
      make_gg_rose(angles[["ES"]][[tp]],
                   title    = paste("ES -", tp),
                   subtitle = sub,
                   caption  = cap,
                   col_fill   = "#ffeeee",
                   col_border = "coral")
    ))
  }

  pdf(outfile, width = 14, height = 7)
  grid.arrange(grobs = plot_list, nrow = 2,
               top = paste("Primary cilia orientation -", dataset_label))
  dev.off()
  cat("Saved:", outfile, "\n")
}

# Combined PDF — one pair of rows per dataset
make_combined_plots <- function(all_angles, all_perm, all_labels, outfile) {
  timepoints <- c("T0", "T4", "T8", "T12")
  all_plots  <- list()

  for (i in seq_along(all_angles)) {
    angles       <- all_angles[[i]]
    perm_results <- all_perm[[i]]
    label        <- all_labels[[i]]

    for (tp in timepoints) {
      n_cilia    <- length(angles[["Control"]][[tp]])
      rho_val    <- round(rho.circular(angles[["Control"]][[tp]]), 3)
      rayleigh_p <- round(rayleigh.test(angles[["Control"]][[tp]])$p.value, 3)
      all_plots <- c(all_plots, list(
        make_gg_rose(angles[["Control"]][[tp]],
                     title    = paste0("Control - ", tp, " [", label, "]"),
                     subtitle = paste0("n = ", n_cilia, ",  rho = ", rho_val,
                                       ",  Rayleigh p = ", format.pval(rayleigh_p, digits = 2)),
                     col_fill   = "#ddeeff",
                     col_border = "steelblue")
      ))
    }

    for (tp in timepoints) {
      n_cilia    <- length(angles[["ES"]][[tp]])
      rho_val    <- round(rho.circular(angles[["ES"]][[tp]]), 3)
      rayleigh_p <- round(rayleigh.test(angles[["ES"]][[tp]])$p.value, 3)
      pr         <- perm_results[perm_results$timepoint == tp, ]
      sub <- paste0("n = ", n_cilia, ",  rho = ", rho_val,
                    ",  Rayleigh p = ", format.pval(rayleigh_p, digits = 2))
      cap <- if (!is.na(pr$U2)) paste0("Watson U2 p (Bonferroni) = ", format.pval(pr$p_bonferroni, digits = 2)) else ""
      all_plots <- c(all_plots, list(
        make_gg_rose(angles[["ES"]][[tp]],
                     title    = paste0("ES - ", tp, " [", label, "]"),
                     subtitle = sub,
                     caption  = cap,
                     col_fill   = "#ffeeee",
                     col_border = "coral")
      ))
    }
  }

  n_rows <- length(all_angles) * 2
  pdf(outfile, width = 14, height = 3.5 * length(all_angles))
  grid.arrange(grobs = all_plots, nrow = n_rows,
               top = "Primary cilia orientation - Control vs electrical stimulation")
  dev.off()
  cat("Saved:", outfile, "\n")
}

# ------------------------------------------------------------
# 3. Plot-only function — skips permutations, just redraws PDFs
#    Use after sourcing the script to tweak plots without waiting
# ------------------------------------------------------------

plot_only <- function() {
  dir.create("figures", showWarnings = FALSE)
  dummy_perm <- data.frame(
    timepoint    = c("T0", "T4", "T8", "T12"),
    U2           = NA, p_value      = NA,
    p_bonferroni = NA, p_fdr        = NA
  )
  for (i in seq_along(datasets)) {
    ang <- load_dataset(datasets[[i]]$path)
    make_rose_plots(ang, dummy_perm, datasets[[i]]$label, datasets[[i]]$file)
  }
  cat("Figures regenerated (no p-values shown — run full script for statistics).\n")
}

# ------------------------------------------------------------
# 4. Define datasets
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
# 5. Run analysis and generate plots for each dataset
# ------------------------------------------------------------

all_angles <- list()
all_perm   <- list()
all_labels <- list()

dummy_perm <- data.frame(
  timepoint    = c("T0", "T4", "T8", "T12"),
  U2           = NA, p_value      = NA,
  p_bonferroni = NA, p_fdr        = NA
)

for (ds in datasets) {
  angles <- load_dataset(ds$path)
  all_angles <- c(all_angles, list(angles))
  all_labels <- c(all_labels, list(ds$label))

  if (PLOT_ONLY) {
    all_perm <- c(all_perm, list(dummy_perm))
    make_rose_plots(angles, dummy_perm, ds$label, ds$file)
  } else {
    results <- run_analysis(angles, ds$label)
    all_perm <- c(all_perm, list(results$perm_results))
    make_rose_plots(angles, results$perm_results, ds$label, ds$file)
  }
}

# ------------------------------------------------------------
# 6. Combined figure
# ------------------------------------------------------------

make_combined_plots(all_angles, all_perm, all_labels,
                    "figures/rose_plots_combined.pdf")

cat("\nAll done. Figures written to figures/\n")
if (PLOT_ONLY) cat("Note: PLOT_ONLY = TRUE, no statistics computed. Set to FALSE for full analysis.\n")
