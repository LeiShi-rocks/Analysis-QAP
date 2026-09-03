# Analysis-QAP

Replication code for the paper on permutation tests for matrix regression (Quadratic Assignment Procedure / Mantel test).

---

## Repository structure

```
Analysis-QAP/
├── QAPpro.R          # Core method implementation
├── Fig1.R            # Replication: Figure 1
├── Fig2.R            # Replication: Figure 2
├── Table4.R          # Replication: Table 4
├── results/
│   ├── Fig1/         # Output figures and simulation data for Figure 1
│   ├── Fig2/         # Output figures and simulation data for Figure 2
│   └── DepressionData/   # Real-data application output (Table 4)
└── legacy/           # Archived earlier scripts and figures
```

### `QAPpro.R` — Method implementation

Contains the core functions used throughout all replication scripts:

- **`dyadicLM()`** — Fits a dyadic linear regression and returns coefficient estimates and a cluster-robust variance matrix (CR0).
- **`QAPpro()`** — Runs a QAP permutation test with three permutation modes (`permute-X`, `permute-Y`, `permute-e`), returns permutation p-values alongside asymptotic ones, and optionally plots the permutation distribution.
- **`super.pop()`** (defined in `Fig1.R`) — Computes the U-statistic-based Mantel correlation and its variance estimator.

---

## Replication scripts

### `Fig1.R` — Permutation distribution of QAP

Simulates the super-population and permutation distributions of the (non-studentized and studentized) Mantel correlation under two kernels:

- **Conservative case** — Walsh averages kernel with circular uniform data
- **Anti-conservative case** — Walsh averages kernel with heavy-tailed (`sinh`) data

Outputs saved to `results/Fig1/`:

| File | Description |
|---|---|
| `super_non_student.png` | Super-population distribution, non-studentized |
| `super_student.png` | Super-population distribution, studentized |
| `perm_non_student.png` | Permutation distribution, non-studentized (conservative) |
| `perm_student.png` | Permutation distribution, studentized (conservative) |
| `super_non_student_anti.png` | Super-population distribution, non-studentized (anti-conservative) |
| `super_student_anti.png` | Super-population distribution, studentized (anti-conservative) |
| `perm_non_student_anti.png` | Permutation distribution, non-studentized (anti-conservative) |
| `perm_student_anti.png` | Permutation distribution, studentized (anti-conservative) |
| `record_super_conservative.rds` | Simulation draws, super-population (conservative) |
| `record_perm_conservative.rds` | Simulation draws, permutation (conservative) |
| `record_super_anti.rds` | Simulation draws, super-population (anti-conservative) |
| `record_perm_anti.rds` | Simulation draws, permutation (anti-conservative) |

**Runtime:** ~30–40 minutes (two loops of MC = 2,000, n = 500 matrices).

---

### `Fig2.R` — Sampling and permutation distributions for MRQAP

Simulates the sampling distribution and permutation distribution of the MRQAP coefficient estimator and Wald statistic under a Walsh averages kernel with correlated predictors (ρ = 0.5).

Outputs saved to `results/Fig2/`:

| File | Description |
|---|---|
| `sampling_NS_MAQAP.png` | Sampling distribution, non-studentized |
| `sampling_S_MAQAP.png` | Sampling distribution, studentized |
| `perm_NS_MAQAP.png` | Permutation distribution, non-studentized |
| `perm_S_MAQAP.png` | Permutation distribution, Wald statistic vs χ²(1) |
| `record_stat_MRQAP.rds` | Simulation draws, sampling distribution |
| `record_perm_MRQAP.rds` | Simulation draws, permutation distribution |

**Runtime:** ~40 minutes (two loops of MC = 2,000, n = 250 matrices, with `dyadicLM` per iteration).

---

### `Table4.R` — Real-data application: depression and social isolation

Applies `QAPpro()` to the depression network dataset to test the depression-isolation and depression-homophily hypotheses. The main analysis uses `mode = "permute-X"` with 2,000 permutations.

Outputs saved to `results/DepressionData/`:

| File | Description |
|---|---|
| `Table4_results.csv` | Numerical output: coefficients and normal, non-studentized QAP, and studentized QAP p-values. The first two rows reproduce Table 4; the remaining rows are the interaction-term analyses discussed in the text. |
| `DPLevel.png` | Permutation distribution for depression mean coefficient |
| `DPSimilarity.png` | Permutation distribution for depression similarity coefficient |

**Runtime:** ~1–2 minutes.

---

## How to run replications

### Requirements

R (≥ 4.0) with the following packages:

```r
install.packages(c("dplyr", "tidyverse", "ggplot2", "mvtnorm",
                   "latex2exp", "clubSandwich", "car"))
```

### Running the scripts

The scripts use relative paths. Set the working directory to your local clone of the repository root (for example, `setwd("/path/to/Analysis-QAP")`) before running them:

```bash
Rscript Fig1.R     # ~30–40 min
Rscript Fig2.R     # ~40 min
Rscript Table4.R   # ~1–2 min
```

Or, from an R session started in the repository root:

```r
source("Fig1.R")
source("Fig2.R")
source("Table4.R")
```

### Reproducibility

All stochastic computations use fixed seed 2024. In `Table4.R`, the main QAP analysis
passes `user.seed = 2024` to `QAPpro()`, so its 2,000 permutations and the
saved `Table4_results.csv` file are exactly reproducible.

### Regenerating figures only (without re-running simulations)

`Fig1.R` and `Fig2.R` save intermediate simulation results as `.rds` files in `results/`. If you only need to regenerate figures after tweaking plot code, load the saved data directly:

```r
# Example: regenerate Fig1 permutation figures
record_perm  <- readRDS("results/Fig1/record_perm_conservative.rds")
record_super <- readRDS("results/Fig1/record_super_conservative.rds")
```

Then re-run only the visualization blocks at the bottom of each script.
