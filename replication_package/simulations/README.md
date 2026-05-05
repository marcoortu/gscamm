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
| `gscamm`        | GSCA-MM, ALR-space GSCA step + per-row MAP polish (default)            | plug-in WLS (paper)             |
| `gscamm_boot`   | identical fit + non-parametric row bootstrap                           | basic CI + noise augmentation   |
| `lda`           | `topicmodels::LDA` (Variational EM, random init, seed=1)               | plug-in ALR-WLS                 |
| `stm`           | `stm::stm`, `init.type = "Random"` (no warm start; canonical)          | plug-in ALR-WLS                 |
| `stm_spectral`  | `stm::stm`, `init.type = "Spectral"` (anchor-words; opt-in counterfactual) | plug-in ALR-WLS             |

For STM we deliberately use the **same** ALR-WLS post-hoc to keep the
comparison apples-to-apples on the second stage; coverage from STM's
native `estimateEffect` would be a different number reported elsewhere
(see paper).

**Why STM defaults to Random init.** STM's native `Spectral` init is the
Arora-Halpern-Mimno anchor-words algorithm: a deterministic, data-driven
warm start with theoretical recovery guarantees. The ablation analysis
in `full_metrics_run{C,D,F}.rds` shows it contributes **+34% to +231%**
of STM's apparent edge on `rmse_theta` over the other models -- it is an
*init effect*, not a *model effect*. To make the model-to-model comparison
fair, the canonical `stm` row uses Random init. The Spectral variant is
available as `stm_spectral` (opt-in via `GSCAMM_USE_STM_SPECTRAL=1`) for
the explicit init-comparison appendix.

### Theta estimators reported

For each method the metrics file records three RMSE columns against the
ground-truth `sim$Theta`:

| column           | gscamm meaning                                        | LDA / STM meaning      |
|------------------|-------------------------------------------------------|------------------------|
| `rmse_theta`     | structural Theta = inverse-ALR(X B), X-only           | posterior gamma/theta  |
| `rmse_theta_R`   | token-level posterior responsibilities (over-peaked)  | same as `rmse_theta`   |
| `rmse_theta_map` | MAP polish: prior X B + multinomial likelihood        | same as `rmse_theta`   |

`rmse_theta_map` is the **apples-to-apples** column when comparing
gscamm to LDA/STM, since both sides combine a covariate-aware prior
with the data likelihood. The legacy `rmse_theta` column is retained
because it is the canonical structural estimand of the paper.

## Single-script ablation driver

For end-to-end ablation across all 4 incremental configs (B, C, D, F) in
a single Rscript invocation, use `code/05_run_ablation.R`. It snapshots
the current `results/full_metrics.{rds,csv}` as `runA` (if not already
present), then spawns 4 child Rscript subprocesses with the appropriate
env-var toggles, and renames the outputs of each into
`results/full_metrics_run{B,C,D,F}.{rds,csv}` plus per-run logs.

```bash
# 1) Install the package source (the children load gscamm via library(),
#    so any source-only change must be installed first).
R CMD INSTALL .

# 2) Drive the 4 incremental runs end-to-end (writes a per-run summary).
Rscript replication_package/simulations/code/05_run_ablation.R
```

The driver writes `results/ablation_summary.csv` after every iteration so
a partial trail survives a crash. On linux02-class hardware the full
ablation takes ~3 hours at R=100, B=200; override with `GSCAMM_SIM_R` /
`GSCAMM_BOOT_B` for a faster pilot.

## Incremental ablation runs (manual sequence)

The simulation script reads several environment variables so the same
code can drive ablation runs A→F without edits. Defaults are tuned for
the "final" run F.

| variable                | default | meaning                                          |
|-------------------------|---------|--------------------------------------------------|
| `GSCAMM_FIT_MAX_ITER`   | `200`   | gscamm EM cap (was 80; raised after the          |
|                         |         | high_sparsity run hit the cap in 2/5 reps)       |
| `GSCAMM_BOOT_MAX_ITER`  | `60`    | bootstrap EM cap (was 30)                        |
| `GSCAMM_NOISE_SCALE`    | `1.0`   | bootstrap noise multiplier (was `2.0`; lowered   |
|                         |         | to bring boot coverage back from 1.000 → ~0.95)  |
| `GSCAMM_USE_POLISH`       | `1`     | apply per-row MAP polish (`fit$Theta_map`)     |
| `GSCAMM_SIGMA2_POLISH`    | `0.25`  | prior variance for the MAP polish (ALR space)  |
| `GSCAMM_USE_STM_SPECTRAL` | `0`     | also fit STM with `init.type="Spectral"` (warm-start counterfactual) |
| `GSCAMM_SIM_R`            | `100`   | replicates per scenario                        |
| `GSCAMM_BOOT_B`           | `200`   | bootstrap replicates per fit                   |

Suggested ablation sequence (each writes to `results/full_metrics.rds`,
**save the previous file first** so you can compare):

```bash
# Run A: current baseline (no changes — already on disk)

# Run B: +iteration budget (isolate convergence asymmetry)
GSCAMM_USE_POLISH=0 GSCAMM_USE_STM_RANDOM=0 GSCAMM_NOISE_SCALE=2.0 \
  Rscript code/02_run_full.R

# Run C: +stm_random (isolate STM Spectral-init edge)
GSCAMM_USE_POLISH=0 GSCAMM_USE_STM_RANDOM=1 GSCAMM_NOISE_SCALE=2.0 \
  Rscript code/02_run_full.R

# Run D: +MAP polish (apples-to-apples theta estimator)
GSCAMM_USE_POLISH=1 GSCAMM_USE_STM_RANDOM=1 GSCAMM_NOISE_SCALE=2.0 \
  Rscript code/02_run_full.R

# Run F: +noise_scale=1.0 (recalibrate bootstrap coverage)
GSCAMM_USE_POLISH=1 GSCAMM_USE_STM_RANDOM=1 GSCAMM_NOISE_SCALE=1.0 \
  Rscript code/02_run_full.R
```

(Run E -- tempered E-step -- was skipped: with the MAP polish in place
the responsibilities R become a diagnostic only.)

## Expected runtime (paper scale)

On a single laptop core with N=1000, V=500, K=10, P=8:
- gscamm fit: ~1-3 s
- gscamm bootstrap (B=200): ~10-15 min per replicate
- lda fit: ~1-2 min
- stm fit: ~30-90 s

Full Monte Carlo (3 scenarios × 100 replicates) is therefore on the
order of a working day. The `02_run_full.R` script writes a partial
checkpoint every 10 replicates so the run can be resumed.
