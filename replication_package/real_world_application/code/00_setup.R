## ---------------------------------------------------------------------------
## 00_setup.R -- common settings, paths, helper functions, and pre-processing
## utilities used across the real-world-application scripts (01_preprocess.R,
## 02_search_K.R, 03_fit.R, 04_topic_terms.R, 05_covariate_effects.R,
## 06_robustness.R, 07_make_tables.R).
##
## Run from the project root.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(gscamm)
  library(quanteda)
  library(quanteda.textstats)
  library(stopwords)
  ## NOTE: do NOT library(textstem) -- it pulls in koRpus which masks
  ## quanteda::tokens(). We instead call textstem::lemmatize_words()
  ## qualified, which loads only what we need.
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(Matrix)
  library(ggplot2)
})

## ---------------------------------------------------------------------------
## Paths -- resolved by walking up from the entry-point script (top-level
## Rscript --file=... argument) when available, falling back to the
## conventional layout under getwd() otherwise. The fallback assumes the
## script is run from the project root.
## ---------------------------------------------------------------------------
.find_setup_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- regmatches(args, regexpr("(?<=--file=).*", args, perl = TRUE))
  if (length(m)) {
    entry <- normalizePath(m[1], mustWork = FALSE)
    candidate <- dirname(entry)
    ## entry-point may live in code/; if 00_setup.R is here, return it.
    if (file.exists(file.path(candidate, "00_setup.R"))) return(candidate)
  }
  ## fallback: assume we are in the project root or in any subdir that
  ## contains the canonical replication_package layout.
  cwd <- getwd()
  cand <- file.path(cwd, "replication_package/real_world_application/code")
  if (file.exists(file.path(cand, "00_setup.R"))) return(cand)
  ## last resort: search upwards
  d <- cwd
  while (d != dirname(d)) {
    cand <- file.path(d, "replication_package/real_world_application/code")
    if (file.exists(file.path(cand, "00_setup.R"))) return(cand)
    d <- dirname(d)
  }
  stop("Could not locate replication_package/real_world_application/code/")
}
HERE <- .find_setup_dir()

REPL_ROOT  <- normalizePath(file.path(HERE, ".."))    ## .../real_world_application
RAW_CSV    <- normalizePath(file.path(REPL_ROOT,
              "old_gscatm/usa_election_speeches_2024.csv"))
DATA_DIR   <- file.path(REPL_ROOT, "data")
RESULT_DIR <- file.path(REPL_ROOT, "results")
FIGURE_DIR <- file.path(REPL_ROOT, "figures")
for (d in c(DATA_DIR, RESULT_DIR, FIGURE_DIR))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

## ---------------------------------------------------------------------------
## Custom stopwords for political speeches
## ---------------------------------------------------------------------------
## Speaker names (the four candidates plus first names): they tag who is
## speaking but carry no topical information.
SPEAKER_NAMES <- c("trump", "donald", "biden", "joe",
                   "harris", "kamala", "pence", "mike")

## Filler words that dominate spoken political English without topical
## content. Carefully curated: words that, if left in, end up at the top of
## several topics' phi vectors and noise the interpretation.
POLITICAL_FILLER <- c(
  "going", "got", "gotta", "really", "people", "know", "think", "lot",
  "like", "look", "looking", "say", "said", "saying", "want", "wanted",
  "go", "get", "gets", "getting", "make", "making", "made",
  "right", "now", "well", "yes", "no", "good", "great", "greatest",
  "hello", "ladies", "gentlemen", "thank", "thanks", "thanking",
  "applause", "laughter", "audience", "cheers", "okay", "ok",
  "one", "two", "three", "first", "second", "many", "much", "lots",
  "every", "everyone", "everybody", "anybody", "somebody", "nobody",
  "thing", "things", "way", "ways", "back", "still", "even",
  "ever", "never", "always", "sometimes", "really",
  "let", "lets", "must", "may", "might", "shall", "would",
  "us", "u.s", "yeah", "huh", "uh", "um"
)

## Temporal deictics: anchor a speech in time but are non-topical (and
## already partially captured by the year covariate).
TEMPORAL_DEICTICS <- c("today", "tonight", "yesterday", "tomorrow",
                       "year", "years", "month", "months",
                       "week", "weeks", "day", "days",
                       "morning", "evening", "afternoon", "night")

CUSTOM_STOPWORDS <- unique(c(SPEAKER_NAMES, POLITICAL_FILLER,
                             TEMPORAL_DEICTICS))

## ---------------------------------------------------------------------------
## clean_speech() -- regex pre-tokenization cleanup of a raw speech string
##
## Drops audience markers, speaker turn tags, dialogue cues, URLs, and
## emails. Lowercasing is left to quanteda::tokens() to keep the function
## composable.
## ---------------------------------------------------------------------------
clean_speech <- function(text) {
  if (is.na(text) || !nzchar(text)) return("")
  ## audience reaction tags: [applause], (laughter), {cheers}, etc.
  text <- str_replace_all(text, "\\[[^\\]]*\\]", " ")
  text <- str_replace_all(text, "\\([^\\)]*\\)", " ")
  text <- str_replace_all(text, "\\{[^\\}]*\\}", " ")
  ## speaker turn markers: "THE PRESIDENT:", "Q.", "MR. SMITH:"
  text <- str_replace_all(text, "^[A-Z][A-Z\\.\\s]+:\\s*", " ")
  text <- str_replace_all(text, "(?m)^[A-Z][A-Z\\.\\s]+:\\s*", " ")
  text <- str_replace_all(text, "(?m)^Q\\.\\s+", " ")
  text <- str_replace_all(text, "(?m)^A\\.\\s+", " ")
  ## URLs and emails
  text <- str_replace_all(text, "https?://\\S+", " ")
  text <- str_replace_all(text, "www\\.\\S+", " ")
  text <- str_replace_all(text,
    "\\b[[:alnum:]._%+-]+@[[:alnum:].-]+\\.[A-Za-z]{2,}\\b", " ")
  ## collapse whitespace
  text <- str_squish(text)
  text
}

## ---------------------------------------------------------------------------
## build_token_pipeline() -- end-to-end tokens object from raw text
##
## Stages:
##   1. clean_speech() regex pre-pass
##   2. quanteda::tokens with strict filtering (punct, symbols, numbers, urls)
##   3. lowercase
##   4. remove smart-list stopwords + custom (speakers / filler / deictics)
##   5. drop tokens shorter than min_nchar
##   6. lemmatize via textstem::lemmatize_words
##
## Bigram detection / compounding is handled separately in build_dfm() so the
## collocation statistics can be computed on the cleaned-and-lemmatized
## token stream.
## ---------------------------------------------------------------------------
build_token_pipeline <- function(text_vec, min_nchar = 3L) {
  cleaned <- vapply(text_vec, clean_speech, character(1), USE.NAMES = FALSE)
  toks <- tokens(cleaned,
                 remove_punct   = TRUE,
                 remove_symbols = TRUE,
                 remove_numbers = TRUE,
                 remove_url     = TRUE,
                 split_hyphens  = TRUE)
  toks <- tokens_tolower(toks)
  toks <- tokens_remove(toks, stopwords::stopwords("en", source = "smart"))
  toks <- tokens_remove(toks, CUSTOM_STOPWORDS)
  toks <- tokens_select(toks, min_nchar = min_nchar)
  ## lemmatize: build a lookup from unique types and use tokens_replace
  types <- unique(unlist(as.list(toks), use.names = FALSE))
  if (length(types)) {
    lemmas <- textstem::lemmatize_words(types)
    keep_map <- types != lemmas
    if (any(keep_map))
      toks <- tokens_replace(toks,
                             pattern = types[keep_map],
                             replacement = lemmas[keep_map],
                             valuetype = "fixed")
  }
  toks
}

## ---------------------------------------------------------------------------
## detect_and_compound_bigrams() -- find significant bigrams via
## quanteda.textstats::textstat_collocations and compound them into single
## tokens (e.g., "joe biden" -> "joe_biden").
##
## Returns the compounded tokens object and the table of selected bigrams
## (for later inspection / paper appendix).
## ---------------------------------------------------------------------------
detect_and_compound_bigrams <- function(toks,
                                        min_count   = 20L,
                                        top_n       = 50L,
                                        min_lambda  = 3.0) {
  cands <- textstat_collocations(toks, size = 2L, min_count = min_count)
  if (!nrow(cands)) {
    return(list(tokens = toks,
                collocations = cands[integer(0), , drop = FALSE]))
  }
  cands <- cands[cands$lambda >= min_lambda, , drop = FALSE]
  cands <- cands[order(-cands$lambda), , drop = FALSE]
  if (nrow(cands) > top_n) cands <- cands[seq_len(top_n), , drop = FALSE]
  if (!nrow(cands))
    return(list(tokens = toks,
                collocations = cands[integer(0), , drop = FALSE]))
  toks_cmp <- tokens_compound(toks, pattern = cands, concatenator = "_")
  list(tokens = toks_cmp, collocations = cands)
}

## ---------------------------------------------------------------------------
## build_dfm() -- assemble a trimmed DFM from a tokens object.
##
## Parameters mirror the choices justified in the manuscript: drop terms
## below min_termfreq, below min_docfreq fraction, or above max_docfreq
## fraction (boilerplate).
## ---------------------------------------------------------------------------
build_dfm <- function(toks,
                      min_termfreq = 20L,
                      min_docfreq  = 0.01,
                      max_docfreq  = 0.6) {
  d <- dfm(toks)
  d <- dfm_trim(d,
                min_termfreq = min_termfreq,
                min_docfreq  = min_docfreq,
                max_docfreq  = max_docfreq,
                docfreq_type = "prop")
  d
}

## ---------------------------------------------------------------------------
## build_covariate_matrix() -- model.matrix with proper reference contrasts.
##
## Old code used ~ . - 1 which gave one column per factor level (rank
## deficient). We use ~ year + party + candidate_type with reference levels
## chosen for interpretability:
##   year:           reference 2019  -> dummies year2020, year2021
##   party:          reference republican -> dummy partydemocratic
##   candidate_type: reference president  -> dummy candidate_typevp
## ---------------------------------------------------------------------------
build_covariate_matrix <- function(df,
                                   year_ref       = "2019",
                                   party_ref      = "republican",
                                   ctype_ref      = "president") {
  dat <- data.frame(
    year           = factor(as.character(df$year),
                            levels = c(year_ref,
                              setdiff(sort(unique(as.character(df$year))),
                                      year_ref))),
    party          = factor(df$party,
                            levels = c(party_ref,
                              setdiff(unique(df$party), party_ref))),
    candidate_type = factor(df$candidate_type,
                            levels = c(ctype_ref,
                              setdiff(unique(df$candidate_type), ctype_ref)))
  )
  X <- stats::model.matrix(~ year + party + candidate_type, data = dat)
  X <- X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE]
  X
}

## ---------------------------------------------------------------------------
## .alr_wls() -- ALR-weighted-least-squares estimator of covariate effects on
## a topic mixture matrix. Used for the LDA and STM comparators where the
## first stage produces theta directly (no structural B). Mirrors the helper
## used by the simulation pipeline (00_setup.R in simulations/).
## ---------------------------------------------------------------------------
.alr_wls <- function(Theta, X_std, ref = ncol(Theta), level = 0.95) {
  Theta <- pmax(Theta, .Machine$double.eps)
  N <- nrow(Theta); K <- ncol(Theta); P <- ncol(X_std)
  D <- cbind(`(Intercept)` = 1, X_std)
  non_ref <- setdiff(seq_len(K), ref)
  K1 <- length(non_ref)
  B  <- matrix(NA_real_, P, K1)
  se <- matrix(NA_real_, P, K1)
  z  <- stats::qnorm(1 - (1 - level) / 2)
  for (j in seq_along(non_ref)) {
    k <- non_ref[j]
    y <- log(Theta[, k] / Theta[, ref])
    w <- Theta[, k] * Theta[, ref]
    f <- stats::lm.wfit(D, y, w)
    bhat <- f$coefficients
    bhat[is.na(bhat)] <- 0
    e <- as.numeric(f$residuals)
    XtWX <- crossprod(D, w * D)
    XtWX_inv <- tryCatch(solve(XtWX), error = function(...) NULL)
    if (is.null(XtWX_inv)) next
    meat <- crossprod(D, (w^2 * e^2) * D)
    Vmat <- XtWX_inv %*% meat %*% XtWX_inv
    B[, j]  <- bhat[-1L]
    se[, j] <- sqrt(pmax(diag(Vmat)[-1L], 0))
  }
  rownames(B)  <- colnames(X_std)
  rownames(se) <- colnames(X_std)
  colnames(B)  <- paste0("comp", non_ref)
  colnames(se) <- paste0("comp", non_ref)
  list(B = B, se = se,
       ci_lo = B - z * se, ci_hi = B + z * se)
}

## ---------------------------------------------------------------------------
## ggplot theme used by all paper figures. Clean monochrome, large enough
## fonts for two-column journal layout.
## ---------------------------------------------------------------------------
theme_paper <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background   = element_rect(fill = "grey95", colour = NA),
      strip.text         = element_text(face = "bold"),
      plot.title         = element_text(face = "bold"),
      plot.subtitle      = element_text(colour = "grey30"),
      legend.position    = "bottom",
      legend.title       = element_text(face = "bold")
    )
}

## ---------------------------------------------------------------------------
## Method palette used across robustness comparison figures.
## ---------------------------------------------------------------------------
METHOD_PALETTE <- c(
  "GSCA-MM" = "#1E3A8A",   ## deep navy
  "LDA+ALR" = "#10B981",   ## green
  "STM"     = "#EF4444"    ## red
)
