## ---------------------------------------------------------------------------
## 04_topic_terms.R -- build the topic-term tables and figures.
##
## Produces two views of each topic's vocabulary distribution:
##   (a) ranked bar charts of top-N terms per topic, panel-faceted
##       (the primary, paper-ready figure -- precise, readable, journal-safe)
##   (b) a word-cloud panel as auxiliary supplement (size = phi_kv)
## And a long table of top-N terms per topic for the LaTeX appendix.
##
## Reads:  results/fit.rds, data/processed_dfm.rds
## Writes: figures/topic_terms_barplot.png
##         figures/topic_terms_wordcloud.png
##         results/topic_top_terms.csv
##         results/topic_top_terms_wide.csv
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
  library(quanteda.textplots)
})

fit <- readRDS(file.path(RESULT_DIR, "fit.rds"))
W   <- as.matrix(readRDS(file.path(DATA_DIR, "processed_dfm.rds")))
vocab <- colnames(W)

K <- fit$K
Phi <- fit$Phi
rownames(Phi) <- paste0("Topic ", seq_len(K))
colnames(Phi) <- vocab

## ---------------------------------------------------------------------------
## Top-N terms per topic
## ---------------------------------------------------------------------------
TOP_N <- 12L
top_terms <- do.call(rbind, lapply(seq_len(K), function(k) {
  ord <- order(Phi[k, ], decreasing = TRUE)[seq_len(TOP_N)]
  data.frame(
    topic = paste0("Topic ", k),
    rank  = seq_len(TOP_N),
    term  = vocab[ord],
    phi   = Phi[k, ord],
    stringsAsFactors = FALSE
  )
}))
write.csv(top_terms,
          file.path(RESULT_DIR, "topic_top_terms.csv"), row.names = FALSE)

## Wide format for a quick eyeball table
top_wide <- do.call(cbind, lapply(seq_len(K), function(k) {
  data.frame(top_terms$term[top_terms$topic == paste0("Topic ", k)],
             stringsAsFactors = FALSE)
}))
colnames(top_wide) <- paste0("Topic_", seq_len(K))
write.csv(top_wide,
          file.path(RESULT_DIR, "topic_top_terms_wide.csv"),
          row.names = FALSE)

## ---------------------------------------------------------------------------
## Figure (a): faceted bar chart of top-N terms per topic.
## Same term can appear in multiple topics, so we build a per-topic unique
## label (term + sep + topic) for the ordering, and strip the suffix at the
## axis text. This is the tidytext::reorder_within idiom inlined.
## ---------------------------------------------------------------------------
top_terms$topic <- factor(top_terms$topic,
                          levels = paste0("Topic ", seq_len(K)))

reorder_within <- function(x, by, within, fun = mean, sep = "___") {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}
strip_suffix <- function(x, sep = "___") gsub(paste0(sep, ".+$"), "", x)

top_terms <- top_terms %>%
  mutate(term_ord = reorder_within(term, phi, topic))

p_bar <- ggplot(top_terms,
                aes(x = phi, y = term_ord, fill = topic)) +
  geom_col(width = 0.7) +
  facet_wrap(~ topic, scales = "free_y", ncol = 3) +
  scale_y_discrete(labels = strip_suffix) +
  scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.85, guide = "none") +
  labs(title = sprintf("Top-%d terms by topic", TOP_N),
       subtitle = sprintf("Bars are component-category probabilities phi_kv (K = %d)",
                          K),
       x = expression(phi[kv]),
       y = NULL) +
  theme_paper() +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 9))

n_rows_facet <- ceiling(K / 3)
ggsave(file.path(FIGURE_DIR, "topic_terms_barplot.png"),
       plot = p_bar,
       width = 10, height = 2.4 * n_rows_facet, dpi = 150,
       limitsize = FALSE)
cat(sprintf("Bar plot saved -> %s\n",
            file.path(FIGURE_DIR, "topic_terms_barplot.png")))

## ---------------------------------------------------------------------------
## Figure (b): wordcloud panel (auxiliary)
##
## We assemble a single PNG with K subplots arranged in a 3-wide grid. Each
## subplot is a wordcloud of top-50 terms per topic, sized by sqrt(phi_kv)
## so the visual weight remains comparable across topics.
## ---------------------------------------------------------------------------
TOP_N_WC <- 50L
n_cols_wc <- 3L
n_rows_wc <- ceiling(K / n_cols_wc)

png(file.path(FIGURE_DIR, "topic_terms_wordcloud.png"),
    width = n_cols_wc * 600, height = n_rows_wc * 480, res = 130)
op <- par(mfrow = c(n_rows_wc, n_cols_wc),
          mar = c(0.5, 0.5, 2.0, 0.5))
on.exit({ par(op); dev.off() })

set.seed(2026)
pal <- viridis::viridis(8, begin = 0.15, end = 0.85)

for (k in seq_len(K)) {
  ord <- order(Phi[k, ], decreasing = TRUE)[seq_len(TOP_N_WC)]
  freqs <- Phi[k, ord]
  ## wordcloud requires named numeric vector; we scale by sqrt to compress
  ## the dynamic range and reduce the dominance of the top-1 term
  freqs_scaled <- sqrt(freqs)
  names(freqs_scaled) <- vocab[ord]
  ## use base wordcloud since quanteda's textplot_wordcloud expects a dfm
  tryCatch(
    wordcloud::wordcloud(words = names(freqs_scaled),
                          freq = freqs_scaled,
                          colors = pal,
                          random.order = FALSE,
                          rot.per = 0.15,
                          scale = c(3.2, 0.7),
                          min.freq = 0),
    error = function(e) {
      ## fallback: text-only display
      plot.new(); text(0.5, 0.5, paste(names(freqs_scaled)[1:8], collapse = "\n"),
                       cex = 0.8)
    })
  title(main = paste0("Topic ", k), cex.main = 1.2, line = 0.5)
}
cat(sprintf("Wordcloud saved -> %s\n",
            file.path(FIGURE_DIR, "topic_terms_wordcloud.png")))
