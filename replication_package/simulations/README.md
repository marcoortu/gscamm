# Simulations

Reproduces Section 4 (Table 1, Figures 1-4) of the GSCA-MM paper, with
two extensions:

1. The proposed bootstrap workflow from Remark 2 is implemented and
   reported as a separate row (`gscamm_boot`) so that both plug-in and
   bootstrap coverage can be compared against the LDA+ALR / STM
   baselines.
2. The bootstrap is **noise-augmented** to account for the latent ALR
   noise that a deterministic link cannot represent. The augmentation
   uses an empirical residual sd estimated from the converged
   responsibilities, with a multiplier `noise_scale` chosen by
   inspection of single-replicate coverage tuning (default `2.0`).

## Files

| Script               | Purpose                                                            |
|----------------------|--------------------------------------------------------------------|
| `code/00_setup.R`    | Common helpers: paths, fit wrappers, alignment, metrics            |
| `code/01_run_pilot.R`| R=5 small-scale Monte Carlo for plumbing checks (~3 min)           |
| `code/02_run_full.R` | R=100 paper-scale Monte Carlo (~hours; checkpoints every 10 reps)  |
| `code/03_make_table1.R` | Aggregates `full_metrics.rds` into `table1.csv` and `table1.tex` |
| `code/04_make_figures.R`| Box-plots of RMSE, coverage, perplexity                         |

Override the size of the full run via environment variables:

```r
Sys.setenv(GSCAMM_SIM_R = "20", GSCAMM_BOOT_B = "100")
source("code/02_run_full.R")
```

## Methods

| label         | algorithm                                                              | inference                       |
|---------------|------------------------------------------------------------------------|---------------------------------|
| `gscamm`      | GSCA-MM, simplex-space GSCA step (paper default)                       | plug-in WLS (paper)             |
| `gscamm_boot` | identical fit + non-parametric row bootstrap                           | basic CI + noise augmentation   |
| `lda`         | `topicmodels::LDA` (Variational EM)                                    | plug-in ALR-WLS                 |
| `stm`         | `stm::stm` with prevalence formula on all covariates                   | plug-in ALR-WLS                 |

For STM we deliberately use the **same** ALR-WLS post-hoc to keep the
comparison apples-to-apples on the second stage; coverage from STM's
native `estimateEffect` would be a different number reported elsewhere
(see paper).

## Expected runtime (paper scale)

On a single laptop core with N=1000, V=500, K=10, P=8:
- gscamm fit: ~1-3 s
- gscamm bootstrap (B=200): ~10-15 min per replicate
- lda fit: ~1-2 min
- stm fit: ~30-90 s

Full Monte Carlo (3 scenarios × 100 replicates) is therefore on the
order of a working day. The `02_run_full.R` script writes a partial
checkpoint every 10 replicates so the run can be resumed.
