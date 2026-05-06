## ---------------------------------------------------------------------------
## 01_preprocess.R -- read the raw speech corpus, clean and tokenize each
## speech, detect significant bigrams, build a trimmed document-feature
## matrix, and emit the processed objects consumed by all later scripts.
##
## Reads:  old_gscatm/usa_election_speeches_2024.csv
## Writes: data/processed_dfm.rds       (sparse N x V count matrix)
##         data/processed_X.rds          (N x P numeric covariate matrix)
##         data/processed_df.rds         (cleaned data frame, N rows)
##         data/processed_collocations.csv  (top bigrams kept)
##         data/preprocess_summary.txt   (human-readable run summary)
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

t0 <- Sys.time()
.log_lines <- character(0)
.log <- function(...) {
  msg <- sprintf("[%s] %s",
                 format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                 paste0(..., collapse = ""))
  cat(msg, "\n", sep = "")
  .log_lines <<- c(.log_lines, msg)
}

## ---------------------------------------------------------------------------
## 1. Read raw CSV and standardize covariates
## ---------------------------------------------------------------------------
.log("Reading raw CSV: ", RAW_CSV)
df_raw <- read.csv(RAW_CSV, stringsAsFactors = FALSE)
.log("  raw rows: ", nrow(df_raw))

## Drop rows with missing essentials. The CSV already carries the
## party / candidate_type / year columns produced by the legacy ETL script;
## we re-validate them defensively.
df_raw <- df_raw %>%
  mutate(
    year = suppressWarnings(as.integer(year)),
    party = tolower(trimws(party)),
    candidate_type = tolower(gsub("\\s+", "_", trimws(candidate_type)))
  ) %>%
  filter(
    !is.na(RawText), nzchar(RawText),
    !is.na(party), party %in% c("republican", "democratic"),
    !is.na(candidate_type),
    candidate_type %in% c("president", "vice_president"),
    !is.na(year)
  )
.log("  rows after covariate validation: ", nrow(df_raw))

## ---------------------------------------------------------------------------
## 2. Pre-tokenization cleaning + tokenization + lemmatization
## ---------------------------------------------------------------------------
.log("Pre-cleaning text and tokenizing (this can take a moment)...")
toks <- build_token_pipeline(df_raw$RawText, min_nchar = 3L)
.log("  total tokens after lemmatization: ",
     sum(ntoken(toks)), "  unique types: ",
     length(unique(unlist(as.list(toks), use.names = FALSE))))

## ---------------------------------------------------------------------------
## 3. Significant bigram detection and compounding
## ---------------------------------------------------------------------------
.log("Detecting significant bigrams (min_count=20, lambda>=3, top 50)...")
bg <- detect_and_compound_bigrams(toks, min_count = 20L,
                                  top_n = 50L, min_lambda = 3.0)
toks <- bg$tokens
.log("  bigrams compounded: ", nrow(bg$collocations))
if (nrow(bg$collocations)) {
  preview <- head(bg$collocations[, c("collocation", "count", "lambda")], 15)
  .log("  top 15 bigrams kept (collocation / count / lambda):")
  for (i in seq_len(nrow(preview)))
    .log(sprintf("    %-30s %6d  %5.2f",
                 preview$collocation[i], preview$count[i], preview$lambda[i]))
}

## ---------------------------------------------------------------------------
## 4. Build trimmed DFM
## ---------------------------------------------------------------------------
.log("Building DFM (min_termfreq=20, min_docfreq=1%, max_docfreq=60%)...")
dfm_full <- build_dfm(toks, min_termfreq = 20L,
                      min_docfreq = 0.01, max_docfreq = 0.6)
.log("  DFM shape: ", nrow(dfm_full), " docs x ", ncol(dfm_full), " types")
.log("  total non-zero entries: ", length(dfm_full@x))

## ---------------------------------------------------------------------------
## 5. Document-length filtering: drop very short, cap very long
## ---------------------------------------------------------------------------
doc_len_post_dfm <- rowSums(dfm_full)
.log("  doc length (post-trim) summary:")
ls_str <- capture.output(summary(as.numeric(doc_len_post_dfm)))
for (l in ls_str) .log("    ", l)

## MIN_DOC_LEN chosen as the lowest cap that preserves balance across all
## four speakers without introducing selection bias toward Trump's longer
## rally-style speeches. At 200 the speaker-specific drop rates were
## ~30/70/71/13% (Trump/Biden/Harris/Pence); at 100 they fall to
## ~27/44/25/8%, with Biden and Harris both retaining > 90 documents.
MIN_DOC_LEN <- 100L
keep_docs <- doc_len_post_dfm >= MIN_DOC_LEN
.log(sprintf("  dropping %d docs with < %d post-trim tokens (keep %d)",
             sum(!keep_docs), MIN_DOC_LEN, sum(keep_docs)))

## Trim the corresponding rows in df and dfm
dfm_kept <- dfm_full[keep_docs, ]
df_kept  <- df_raw[keep_docs, , drop = FALSE]

## Drop any term that became empty (zero column) after dropping short docs
col_keep <- colSums(dfm_kept) > 0
dfm_kept <- dfm_kept[, col_keep]
.log("  final DFM shape: ", nrow(dfm_kept), " docs x ",
     ncol(dfm_kept), " types")

## ---------------------------------------------------------------------------
## 6. Covariate matrix with proper reference contrasts
## ---------------------------------------------------------------------------
X <- build_covariate_matrix(df_kept,
                            year_ref  = "2019",
                            party_ref = "republican",
                            ctype_ref = "president")
.log("  covariate matrix: ", nrow(X), " x ", ncol(X),
     "  cols: ", paste(colnames(X), collapse = ", "))

## ---------------------------------------------------------------------------
## 7. Save processed objects as sparse / efficient formats
## ---------------------------------------------------------------------------
W_sparse <- as(dfm_kept, "dgCMatrix")
saveRDS(W_sparse,    file.path(DATA_DIR, "processed_dfm.rds"))
saveRDS(X,           file.path(DATA_DIR, "processed_X.rds"))
saveRDS(df_kept,     file.path(DATA_DIR, "processed_df.rds"))
write.csv(bg$collocations,
          file.path(DATA_DIR, "processed_collocations.csv"),
          row.names = FALSE)

elapsed <- as.numeric(Sys.time() - t0, units = "secs")
.log(sprintf("Preprocessing complete in %.1f s", elapsed))
.log("Artifacts written to: ", DATA_DIR)

writeLines(.log_lines, file.path(DATA_DIR, "preprocess_summary.txt"))
