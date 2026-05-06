# GSCA-MM

**Generalized Structured Component Analysis for Mixture-of-Multinomials Models**

`gscamm` is the R package accompanying the paper

> Ortu, M. and Frigau, L. (2026). *Structured Component Regression for
> Covariate Effects in Mixture-of-Multinomials Models.*

## About the article

The paper introduces **GSCA-MM**, a structured component regression framework
for mixture-of-multinomials data with unit-level covariates. GSCA-MM
generalizes the GSCA-TM topic model to arbitrary mixture-of-multinomials
applications such as microbiome composition analysis, market basket studies,
and topic modeling.

In the GSCA-MM model, covariate effects on the latent mixture weights are
encoded through a **path coefficient matrix** within a generalized structured
component analysis formulation. Estimation alternates:

1. **EM updates** for the component-category distributions (Phi) and the
   observation-component weights (Theta);
2. A **ridge-penalized log-ratio projection** of the posterior
   responsibilities onto the covariate space, performed in ALR coordinates
   (the GSCA step), with a structural zero on the reference component to
   enforce identifiability of the path coefficients (B).

This decouples latent structure recovery from structural parameter
estimation. Post-estimation inference on covariate effects is provided through
**additive log-ratio weighted least squares (ALR-WLS)**, with a noise-augmented
non-parametric bootstrap variant (Remark 2) recommended for confidence
intervals.

The paper (Section 4) compares GSCA-MM against LDA + ALR-WLS and the
Structural Topic Model (STM) on a Monte Carlo design with three scenarios
(`baseline`, `high_covariate`, `high_sparsity`), reporting RMSE, coverage,
interval width, and perplexity. Section 5 illustrates the framework on a
real-world corpus of US presidential speeches.

## About the R package

`gscamm` provides a clean implementation of the framework. Main entry points:

| Function                           | Purpose                                                          |
|------------------------------------|------------------------------------------------------------------|
| `fit_gscamm()`                     | Fits the EM-GSCA algorithm (Algorithm 1 of the paper)            |
| `covariate_effects()`              | Plug-in ALR-WLS inference on covariate effects (Section 2.4)     |
| `bootstrap_covariate_effects()`    | Noise-augmented bootstrap CIs for path coefficients (Remark 2)   |
| `simulate_gscamm()`                | Synthetic data generator mirroring the simulation design         |
| `search_optimal_components()`      | Component-count selection by perplexity                          |
| `align_components()`, `perplexity()` | Alignment / model evaluation utilities                         |

S3 methods `print()`, `summary()`, `coef()`, `predict()`, and `plot()` are
provided for `gscamm` and `gscamm_effects` objects.

## Installation

The package is hosted on GitHub. Install it with:

```r
# install dependencies first
install.packages(c("Matrix", "clue", "remotes"))

# install gscamm from GitHub
remotes::install_github("marcoortu/gscamm")
```

To install from a local clone of this repository:

```r
# from the parent directory of the cloned repo
remotes::install_local("gscamm")
# or, from inside the repo:
# devtools::install()
```

The simulation comparators additionally require:

```r
install.packages(c("topicmodels", "stm", "slam"))
```

## Reproducing Table 1 and all figures

The replication code lives in `replication_package/simulations/`. From the
project root, in an R session:

```r
# 1. Common setup: paths, fit wrappers, alignment, metrics
source("replication_package/simulations/code/00_setup.R")

# 2. (Optional) quick pilot with R = 5 replicates (~3 min) to check plumbing
source("replication_package/simulations/code/01_run_pilot.R")

# 3. Full Monte Carlo: 3 scenarios x R = 100 replicates (paper scale, ~hours).
#    Writes a checkpoint every 10 replicates so the run can be resumed.
source("replication_package/simulations/code/02_run_full.R")

# 4. Aggregate results into Table 1 (CSV + LaTeX)
source("replication_package/simulations/code/03_make_table1.R")

# 5. Produce all figures (RMSE, coverage, perplexity box-plots)
source("replication_package/simulations/code/04_make_figures.R")
```

Outputs are written to:

- `replication_package/simulations/results/` — `full_metrics.rds`, `table1.csv`,
  `table1.tex`
- `replication_package/simulations/figures/` — PNG plots

To run a smaller version of the full Monte Carlo, override the defaults via
environment variables before sourcing `02_run_full.R`:

```r
Sys.setenv(GSCAMM_SIM_R = "20", GSCAMM_BOOT_B = "100")
source("replication_package/simulations/code/02_run_full.R")
```

The full design (`replication_package/simulations/code/00_setup.R`,
`DEFAULT_DESIGN`) uses N = 1000 observations, V = 500 categories, K = 10
components, P = 8 covariates, R = 100 Monte Carlo replicates, and a bootstrap
of B = 200 resamples per replicate.

### Running from the shell

The same pipeline can be launched non-interactively:

```sh
Rscript replication_package/simulations/code/00_setup.R
Rscript replication_package/simulations/code/02_run_full.R
Rscript replication_package/simulations/code/03_make_table1.R
Rscript replication_package/simulations/code/04_make_figures.R
```

## Methods compared in Table 1

| Label          | Algorithm                                                     | Inference                       |
|----------------|---------------------------------------------------------------|---------------------------------|
| `gscamm`       | GSCA-MM, ALR-space log-ratio projection GSCA step (default)   | plug-in ALR-WLS                 |
| `gscamm_boot`  | Same fit + non-parametric row bootstrap                       | basic CI + noise augmentation   |
| `lda`          | `topicmodels::LDA` (Variational EM)                           | plug-in ALR-WLS                 |
| `stm`          | `stm::stm` with prevalence formula (`init.type = "Random"`)   | plug-in ALR-WLS                 |

## Citation

If you use this package or replicate the simulation study, please cite:

```
@article{ortu_frigau_2026_gscamm,
  author  = {Ortu, Marco and Frigau, Luca},
  title   = {Structured Component Regression for Covariate Effects in
             Mixture-of-Multinomials Models},
  year    = {2026}
}
```

## License

MIT License. Copyright (c) 2026 Marco Ortu and Luca Frigau. See `LICENSE`.

## Authors

- Marco Ortu — University of Cagliari (<marco.ortu@unica.it>),
  ORCID [0000-0003-4191-5058](https://orcid.org/0000-0003-4191-5058)
- Luca Frigau — ORCID [0000-0002-6316-4040](https://orcid.org/0000-0002-6316-4040)

## Issues

Bug reports and feature requests: <https://github.com/marcoortu/gscamm/issues>.
