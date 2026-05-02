# GSCA-MM Replication Package

Replication code for **Ortu and Frigau (2026)**, *Structured Component
Regression for Covariate Effects in Mixture-of-Multinomials Models*,
adapted to the generalized `gscamm` R package.

## Layout

```
replication_package/
├── simulations/
│   ├── code/         R scripts to reproduce Section 4 (Table 1)
│   ├── results/      RData / CSV outputs
│   └── figures/      PNG plots (RMSE, coverage, perplexity)
└── real_world_application/
    ├── code/         R scripts for Section 5 (US presidential speeches)
    ├── data/         input corpus + metadata
    ├── results/      fitted models, coefficient tables
    └── figures/      log-odds plots, etc.
```

## Reproducing the simulations

```r
# from the project root
source("replication_package/simulations/code/00_setup.R")
source("replication_package/simulations/code/01_run_pilot.R")    # quick (~minutes)
source("replication_package/simulations/code/02_run_full.R")      # full R=100 (~hours)
source("replication_package/simulations/code/03_make_table1.R")
source("replication_package/simulations/code/04_make_figures.R")
```

The pilot fits R=10 replicates per scenario for sanity checking and
debugging; the full run reproduces the paper's R=100 Monte Carlo
(Section 4.3) with the addition of bootstrap-based coverage for the
GSCA-MM path coefficients (Remark 2).

## Methods compared

| Label              | Description                                                          |
|--------------------|----------------------------------------------------------------------|
| `gscamm`           | GSCA-MM with simplex-space GSCA step (paper default), plug-in WLS    |
| `gscamm_boot`      | GSCA-MM with noise-augmented basic-bootstrap CIs (recommended)       |
| `lda_alr`          | Latent Dirichlet Allocation followed by ALR-WLS regression           |
| `stm`              | Structural Topic Model with prevalence formula                       |

The `gscamm_boot` row reports identical point estimates to `gscamm`
(after first-order bootstrap bias correction the point may shift
slightly) but uses bootstrap CIs. This matches the proposal in Remark 2
of the paper and is the recommended inferential workflow.
