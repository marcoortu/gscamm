## ---------------------------------------------------------------------------
## 07_make_tables.R -- LaTeX tables for the paper.
##
## Produces three appendix-ready tables:
##   tab_topic_terms.tex     top-N terms per topic (one column per topic)
##   tab_effects.tex         GSCA-MM ALR coefficients, SE, 95% CI, BH p-adj
##   tab_robustness.tex      same coefficients juxtaposed across methods
##                           (GSCA-MM, LDA+ALR, STM)
##
## Reads:  results/{topic_top_terms_wide.csv, effects_long.csv, effects_all.csv}
## Writes: results/tab_*.tex (and the equivalent tab_*.csv)
## ---------------------------------------------------------------------------

.find_setup <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- regmatches(args, regexpr("(?<=--file=).*", args, perl = TRUE))
  here <- if (length(m)) dirname(normalizePath(m[1])) else getwd()
  if (!file.exists(file.path(here, "00_setup.R")))
    here <- file.path(getwd(),
                      "replication_package/real_world_application/code")
  file.path(here, "00_setup.R")
}
source(.find_setup())

suppressPackageStartupMessages({
  library(xtable)
})

## ---------------------------------------------------------------------------
## 1. Topic-terms table: top-N terms per topic, one row per rank.
## ---------------------------------------------------------------------------
tt <- read.csv(file.path(RESULT_DIR, "topic_top_terms_wide.csv"),
               check.names = FALSE)
print(
  xtable(tt,
         caption = paste0("Top terms per topic, ranked by phi_kv. ",
                          "GSCA-MM fit with K = ", ncol(tt),
                          " components on the cleaned campaign-speech corpus."),
         label   = "tab:topic-terms"),
  file = file.path(RESULT_DIR, "tab_topic_terms.tex"),
  include.rownames = TRUE, booktabs = TRUE)

## ---------------------------------------------------------------------------
## 2. GSCA-MM effects: covariate, topic, estimate, SE, CI, p, p.adj.
## ---------------------------------------------------------------------------
co <- read.csv(file.path(RESULT_DIR, "effects_long.csv"))
co <- co[co$covariate != "(Intercept)", , drop = FALSE]

.cov_label <- c(
  "year2020"                       = "Year 2020 vs 2019",
  "year2021"                       = "Year 2021 vs 2019",
  "partydemocratic"                = "Democratic vs Republican",
  "candidate_typevice_president"   = "VP vs President"
)
co <- co %>%
  mutate(
    Covariate = .cov_label[covariate],
    Topic     = component_label,
    Estimate  = round(estimate,  3),
    SE        = round(std.error, 3),
    `95% CI`  = sprintf("[%5.2f, %5.2f]",
                        round(conf.low,  2), round(conf.high, 2)),
    `p-adj BH` = format.pval(p.adj, digits = 2, eps = 1e-3)
  ) %>%
  arrange(Covariate, Topic) %>%
  select(Covariate, Topic, Estimate, SE, `95% CI`, `p-adj BH`)

write.csv(co, file.path(RESULT_DIR, "tab_effects.csv"), row.names = FALSE)
print(
  xtable(co,
         caption = paste0("GSCA-MM covariate effects on topic prevalence ",
                          "(ALR scale). Plug-in WLS standard errors with ",
                          "BH-adjusted p-values."),
         label   = "tab:gscamm-effects"),
  file = file.path(RESULT_DIR, "tab_effects.tex"),
  include.rownames = FALSE, booktabs = TRUE,
  sanitize.text.function = function(x) x)

## ---------------------------------------------------------------------------
## 3. Robustness table: GSCA-MM vs LDA vs STM coefficients side by side.
## ---------------------------------------------------------------------------
ea <- read.csv(file.path(RESULT_DIR, "effects_all.csv"))
ea <- ea[ea$covariate != "(Intercept)", , drop = FALSE]
ea$Covariate <- .cov_label[ea$covariate]
.n_topics_m1 <- length(unique(ea$component))
ea$Topic     <- factor(ea$component,
                       levels = paste0("comp", seq_len(.n_topics_m1)),
                       labels = paste0("Topic ", seq_len(.n_topics_m1)))

ea_summary <- ea %>%
  mutate(value = sprintf("%5.2f%s",
                         round(estimate, 2),
                         ifelse(p.adj < 0.05, "*", " "))) %>%
  select(method, Covariate, Topic, value) %>%
  tidyr::pivot_wider(names_from = method, values_from = value) %>%
  arrange(Covariate, Topic)

write.csv(ea_summary, file.path(RESULT_DIR, "tab_robustness.csv"),
          row.names = FALSE)
print(
  xtable(ea_summary,
         caption = paste0("ALR coefficients across methods (* indicates ",
                          "BH-adjusted p < 0.05). Topics aligned to GSCA-MM ",
                          "via the Hungarian algorithm."),
         label   = "tab:robustness"),
  file = file.path(RESULT_DIR, "tab_robustness.tex"),
  include.rownames = FALSE, booktabs = TRUE,
  sanitize.text.function = function(x) x)

cat("LaTeX tables saved to:\n")
for (f in c("tab_topic_terms.tex", "tab_effects.tex",
            "tab_robustness.tex"))
  cat("  ", file.path(RESULT_DIR, f), "\n")
