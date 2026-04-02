# Cilia Circular Statistics

R scripts for comparing the distribution of primary cilia angles between control and electrically stimulated conditions across multiple timepoints, using circular statistical methods.

## Background

Primary cilia angles are **directional data** — 0° and 360° are the same point, and arithmetic wraps around. Standard statistical tests (t-test, ANOVA) are therefore inappropriate. This repository uses the `circular` R package with permutation-based testing.

In the absence of a directional cue, primary cilia show no preferred orientation, producing a uniform circular distribution. Electrical stimulation (ES) is hypothesised to induce planar cell polarity (PCP), shifting cilia toward a concentrated directional distribution.

## Experimental design

| Factor | Levels |
|---|---|
| Condition | Control, ES (electrical stimulation) |
| Timepoint | T=0, T=4, T=8, T=12 |
| ES intensity | 100 mV/mm, 25 mV/mm |
| Angle range | −180° to +180° |

## Repository structure

```
cilia-circular-stats/
├── R/
│   └── cilia_circular_analysis.R    # Main analysis script
├── data/
│   ├── 100mvmm_Cilia_angles_25_9_24.csv
│   └── 25mvmm_Cilia_angles_30_9_24.csv
├── figures/                          # Created on first run
│   ├── rose_plots_100mVmm.pdf
│   ├── rose_plots_25mVmm.pdf
│   └── rose_plots_combined.pdf
├── docs/
│   └── methods.md
└── README.md
```

## Requirements

R (≥ 4.0) with the following packages. Install from the Linux terminal:

```bash
# Option 1 — system package manager
sudo apt install r-cran-circular r-cran-ggplot2 r-cran-gridextra

# Option 2 — from CRAN
sudo apt install r-base-dev gfortran
Rscript -e "install.packages(c('circular', 'ggplot2', 'gridExtra'), repos='https://cloud.r-project.org')"
```

## Usage

```r
setwd("/path/to/cilia-circular-stats")
source("R/cilia_circular_analysis.R")
```

Or from the terminal:

```bash
cd cilia-circular-stats
Rscript R/cilia_circular_analysis.R
```

The script runs both datasets sequentially. Permutations (9,999 per timepoint per dataset) take a few minutes in total.

## Output

**Console** (for each dataset):
- Sample sizes per group
- Rayleigh test results (uniformity check)
- Permutation Watson U² results (Control vs ES at each timepoint)
- Bonferroni and FDR-corrected p-values
- Summary statistics (mean direction, resultant length ρ, circular SD)

**figures/**:
- `rose_plots_100mVmm.pdf` — 2 × 4 rose plots for 100 mV/mm dataset
- `rose_plots_25mVmm.pdf` — 2 × 4 rose plots for 25 mV/mm dataset
- `rose_plots_combined.pdf` — all datasets in a single figure

## Adding new datasets

Add a new entry to the `datasets` list near the bottom of the script:

```r
list(
  path  = "data/your_new_file.csv",
  label = "50 mV/mm",
  file  = "figures/rose_plots_50mVmm.pdf"
)
```

The CSV must follow the same format as the existing files (see below).

## Data format

One column per condition-timepoint. Columns may have unequal numbers of rows. Angles in degrees, −180° to +180°.

```
Control T=0 | Control T=4 | Control T=8 | Control T=12 | ES T=0 | ES T=4 | ES T=8 | ES T=12
```

## Statistical approach

See [`docs/methods.md`](docs/methods.md) for full details.

