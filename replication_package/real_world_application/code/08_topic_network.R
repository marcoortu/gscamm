## ---------------------------------------------------------------------------
## 08_topic_network.R -- topic-similarity network. Each node is a topic;
## edges connect topics whose top-N vocabulary distributions overlap (Jaccard
## similarity above a threshold). Useful for spotting redundant or
## co-occurring themes in the topic solution.
##
## Reads:  results/fit.rds, results/topic_top_terms.csv
## Writes: figures/topic_network.png
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
  library(igraph)
  library(ggraph)
  library(tidygraph)
})

fit <- readRDS(file.path(RESULT_DIR, "fit.rds"))
tt  <- read.csv(file.path(RESULT_DIR, "topic_top_terms.csv"),
                stringsAsFactors = FALSE)

K <- fit$K
TOP_N <- 30L

## Per-topic top-N term sets
topics_top <- split(tt$term, tt$topic)
topics_top <- lapply(topics_top, function(v) head(v, TOP_N))

## Pairwise Jaccard similarity
.jaccard <- function(a, b) {
  ix <- length(intersect(a, b)); un <- length(union(a, b))
  if (un == 0) 0 else ix / un
}
M <- matrix(0, K, K, dimnames = list(names(topics_top), names(topics_top)))
for (i in seq_len(K)) for (j in seq_len(K)) {
  if (i == j) next
  M[i, j] <- .jaccard(topics_top[[i]], topics_top[[j]])
}

## Build edge list above a similarity floor (keeps the figure readable)
edges <- as.data.frame(as.table(M), stringsAsFactors = FALSE) %>%
  rename(from = Var1, to = Var2, weight = Freq) %>%
  filter(from < to, weight > 0.05)
cat(sprintf("Edges with Jaccard > 0.05: %d\n", nrow(edges)))

## Node attributes: top term per topic, and total prevalence in the corpus
prev <- colMeans(fit$Theta_map)
nodes <- data.frame(
  name = paste0("Topic ", seq_len(K)),
  top_term = vapply(topics_top, function(v) v[1], character(1)),
  prevalence = prev,
  stringsAsFactors = FALSE
)

g <- as_tbl_graph(igraph::graph_from_data_frame(d = edges, vertices = nodes,
                                                directed = FALSE))

p <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(width = weight, alpha = weight), colour = "grey50",
                 show.legend = FALSE) +
  geom_node_point(aes(size = prevalence), colour = "#1E3A8A") +
  geom_node_text(aes(label = paste0(name, "\n[", top_term, "]")),
                 repel = TRUE, size = 3.5,
                 box.padding = 0.4) +
  scale_edge_width(range = c(0.5, 2.5)) +
  scale_edge_alpha(range = c(0.4, 0.9)) +
  scale_size_continuous(range = c(4, 10), name = "Mean theta") +
  labs(title = "Topic similarity network",
       subtitle = paste0("Edges = Jaccard overlap of top-", TOP_N,
                          " terms; node size = mean topic prevalence"),
       x = NULL, y = NULL) +
  theme_paper() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank())

ggsave(file.path(FIGURE_DIR, "topic_network.png"),
       plot = p, width = 8, height = 6, dpi = 150)
cat(sprintf("Topic network -> %s\n",
            file.path(FIGURE_DIR, "topic_network.png")))
