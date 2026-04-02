# Statistical Methods

## Why circular statistics?

Cilia angles are **circular** (directional) data. Unlike linear measurements,
angular data wrap around: 0° and 360° are the same direction, and −180° and
+180° are identical. Standard means, standard deviations, t-tests, and ANOVA
all assume unbounded linear data and give incorrect results on raw angles.

---

## Data preparation

Raw angles (−180° to +180°) are converted to the 0°–360° range:

```r
x_360 <- ifelse(x < 0, x + 360, x)
```

They are then wrapped into a `circular` object specifying degrees and 2π
modulus, so all subsequent functions treat the values as directional.

---

## Test selection

### Step 1 — Rayleigh test of uniformity

Run first on every group to determine whether data are concentrated around a
mean direction or uniformly spread around the circle.

- **p < 0.05**: data are concentrated — a mean direction exists
- **p ≥ 0.05**: data appear uniform — no preferred direction

**Result:** Control groups were uniform at all four timepoints (p = 0.07–0.88).
ES groups were strongly concentrated (p < 0.0001). Because the two conditions
have fundamentally different distributional shapes, the Watson-Williams test
(which assumes equal concentration and von Mises distributions in both groups)
is not appropriate.

### Step 2 — Watson U² test with permutation p-values

The Watson U² statistic is a non-parametric, distribution-free measure of
whether two circular samples come from the same population. It makes no
assumption about the shape of the underlying distribution and handles the
case where one group is uniform and the other is concentrated.

**Why permutation?** The `circular` R package computes Watson U² p-values
from a lookup table with a floor of p < 0.001 — it does not return a numeric
p-value. To obtain exact, numeric p-values we use a permutation approach:

1. Compute the observed U² statistic from the real group labels
2. Randomly shuffle the group labels 9,999 times, recomputing U² each time
3. The permutation p-value is the proportion of shuffled statistics ≥ the
   observed statistic, with a +1 correction to avoid p = 0:

```
p = (number of permuted U² >= observed U² + 1) / (n_perm + 1)
```

The minimum resolvable p-value with 9,999 permutations is 0.0001. All four
observed U² values (0.69–1.15) exceeded every single permuted statistic,
giving p = 0.0001 at all timepoints. The true p-value is ≤ 0.0001.

---

## Multiple comparison correction

Four simultaneous tests were conducted (one per timepoint). Two corrections
are reported:

| Method | Controls | Notes |
|---|---|---|
| Bonferroni | Family-wise error rate (FWER) | Multiplies each p by 4 |
| Benjamini-Hochberg (FDR) | Expected false discovery rate | Less conservative |

With p = 0.0001 at all timepoints, both corrected p-values (0.0004 and 0.0001
respectively) remain highly significant. Correction does not change the
conclusions but is reported for completeness.

---

## Summary statistics

For each condition-timepoint group:

| Statistic | Description |
|---|---|
| Mean direction (μ) | Circular mean — direction of the resultant vector |
| Resultant length (ρ) | 0 = uniform, 1 = all angles identical; measures concentration |
| Circular SD | Increases as ρ decreases |

Control ρ values near 0 confirm uniform distribution. ES ρ values
substantially above 0 confirm a consistent directional preference.

---

## Software

- R (≥ 4.0)
- `circular` package (Agostinelli & Lund, 2023)

## References

- Mardia, K. V. & Jupp, P. E. (2000). *Directional Statistics*. Wiley.
- Pewsey, A., Neuhäuser, M. & Ruxton, G. D. (2013). *Circular Statistics in R*. Oxford University Press.
- Agostinelli, C. & Lund, U. (2023). R package `circular`. https://cran.r-project.org/package=circular
