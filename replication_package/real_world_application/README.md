# Real-world application

GSCA-MM applied to a corpus of US presidential and vice-presidential
campaign-cycle speeches (2019-2021), modernized from the legacy GSCA-TM
example under `old_gscatm/`.

## Pipeline

| Script | Role |
|---|---|
| `code/00_setup.R`         | Paths, helpers (preprocessing, ALR-WLS, paper theme), method palette |
| `code/01_preprocess.R`    | Read CSV → clean text → DFM → `data/processed_*` |
| `code/02_search_K.R`      | Grid over K=2..18 with held-out perplexity → `figures/search_K.png` |
| `code/03_fit.R`           | Fit GSCA-MM at the chosen K (default K=6) → `results/fit.rds` |
| `code/04_topic_terms.R`   | Top-N terms barplot panel + word cloud auxiliary |
| `code/05_covariate_effects.R` | Forest plot, heatmap, predicted prevalence-by-year |
| `code/06_robustness.R`    | LDA + STM (Random init) comparators; method-comparison plots |
| `code/07_make_tables.R`   | LaTeX tables (top terms, GSCA-MM effects, robustness summary) |
| `code/08_topic_network.R` | Topic-similarity graph (Jaccard on top-30 vocab) |

Run from the project root:

```bash
R CMD INSTALL .
Rscript replication_package/real_world_application/code/01_preprocess.R
Rscript replication_package/real_world_application/code/02_search_K.R
GSCAMM_APP_K=6 Rscript replication_package/real_world_application/code/03_fit.R
Rscript replication_package/real_world_application/code/04_topic_terms.R
Rscript replication_package/real_world_application/code/05_covariate_effects.R
Rscript replication_package/real_world_application/code/06_robustness.R
Rscript replication_package/real_world_application/code/07_make_tables.R
Rscript replication_package/real_world_application/code/08_topic_network.R
```

Total runtime on a single core: ~7 min for the K-search (17 fits at K=2..18),
~30 s for the final fit, ~2 min for the LDA/STM robustness fits, and
sub-second for everything else.

## Why this rewrite differs from `old_gscatm/`

### Cleaner NLP

The old pipeline kept noisy boilerplate vocabulary in the topic
distributions (audience-reaction tags, speaker-turn markers, ubiquitous
filler words). The new pipeline:

1. **Pre-tokenisation regex pass** that strips `[applause]`, `(laughter)`,
   `THE PRESIDENT:`, `Q.`, URLs, and emails.
2. Quanteda tokens with **`remove_punct + remove_symbols + remove_numbers
   + remove_url + split_hyphens`** all enabled.
3. **Smart-list stopwords** (~570 words) plus a custom list of
   speaker names, political-speech filler, and temporal deictics
   (defined in `00_setup.R`).
4. **Lemmatisation** via `textstem::lemmatize_words` so that
   "running"/"runs"/"ran" collapse to a single token.
5. **Significant bigrams** (top-50 collocations with lambda >= 3 and
   count >= 20) compounded into single tokens (e.g. `joe_biden`,
   `puerto_rico`, `executive_order`).
6. **DFM trim**: `min_termfreq = 20`, `min_docfreq = 1%`, `max_docfreq = 60%`.
7. **Document filter**: `>= 100` post-trim tokens (preserves the
   four-speaker balance; threshold 200 over-weighted Trump's long rallies).

### Identifiable covariate matrix

Old code used `model.matrix(~ . - 1)` which produced rank-deficient
factor encodings (one column per level). The new code uses proper
reference contrasts:

- `year`: 2019 reference -> dummies `year2020`, `year2021`
- `party`: republican reference -> `partydemocratic`
- `candidate_type`: president reference -> `candidate_typevp`

Four covariates total, all interpretable.

### K selection by held-out perplexity

`search_optimal_components(..., holdout = 0.10)` reports both in-sample
and held-out perplexity over K=2..18. The held-out minimum gives K* = 3,
with K=6 a close second (1474 vs 1481 in this dataset). We pick K=6 as
the operating point: held-out perplexity essentially identical, but the
five topic contrasts (vs two for K=3) give the covariate analysis enough
granularity. The legacy paper's K=9 sits well above the elbow in the
held-out curve (1546).

### Cleaner figures

| Figure | Purpose |
|---|---|
| `search_K.png`                 | Perplexity vs K, in-sample + held-out, with marker on K* |
| `topic_terms_barplot.png`      | Ranked horizontal bars per topic (within-facet sort), top 12 terms |
| `topic_terms_wordcloud.png`    | Auxiliary word cloud panel for visual scanning |
| `prevalence_by_year.png`       | Predicted mean theta by year, marginalised over party x VP |
| `effects_forest.png`           | Forest plot of ALR coefficients, faceted by covariate |
| `effects_heatmap.png`          | Covariate x topic heatmap with hierarchical clustering on topics |
| `method_comparison_effects.png`| Per-(covariate, topic) coefficient comparison, three colors |
| `method_comparison_signif_count.png` | Per-method significant-effect counts |
| `topic_network.png`            | Topic similarity graph (Jaccard) for redundancy check |

All figures share a single `theme_paper()` and a stable method palette
defined in `00_setup.R`.

### Honest comparator setup

STM uses `init.type = "Random"` (matching the simulation choice). The
spectral anchor-words init is a strong warm-start that systematically
gives STM an inflated edge in apparent recovery; using Random init
isolates the model from the warm-start effect.

## Files

```
data/
  processed_dfm.rds              N x V sparse count matrix (post-clean)
  processed_X.rds                N x P numeric covariate matrix
  processed_df.rds               cleaned data frame, N rows
  processed_collocations.csv     bigrams kept after collocation detection
  preprocess_summary.txt         human-readable preprocessing report

results/
  search_K.csv                   K-search results
  fit.rds                        the gscamm fit object
  fit_summary.txt                fit statistics
  topic_top_terms.csv            top-N terms per topic, long format
  topic_top_terms_wide.csv       same, wide format
  effects_long.csv               GSCA-MM ALR coefficients
  effects_lda.csv, effects_stm.csv, effects_all.csv  robustness fits
  tab_topic_terms.tex / .csv     LaTeX top-terms table
  tab_effects.tex / .csv         LaTeX GSCA-MM effects table
  tab_robustness.tex / .csv      LaTeX cross-method comparison

figures/
  search_K.png
  topic_terms_barplot.png
  topic_terms_wordcloud.png
  prevalence_by_year.png
  effects_forest.png
  effects_heatmap.png
  method_comparison_effects.png
  method_comparison_signif_count.png
  topic_network.png
```
