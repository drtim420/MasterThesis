# ============================================================
# 02_simulate_and_evaluate.R
# Modular pipeline: simulate DAG data + return amat & weights
#                   evaluate any DAG on any dataset via CI tests
# ============================================================

suppressPackageStartupMessages({
  library(dagitty)
  library(igraph)
  library(dplyr)
  library(tibble)
  library(readr)
})

# -------------------------------
# Config
# -------------------------------
set.seed(123)

data_path <- "../data"   # for Sachs real data (cd3cd28.xls)
out_dir   <- "../results"
alpha     <- 0.05

# ensure output dir
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -------------------------------
# Sachs graph helpers (same node order as your main script)
# -------------------------------

nms <- c("Raf","Mek","PLCg","PIP2","PIP3","Erk","Akt","PKA","PKC","p38","JNK")

get_sachs_amat <- function(){
  A <- matrix(c(
    0,0,0,0,0,0,0,1,1,0,0,
    1,0,0,0,0,0,0,1,1,0,0,
    0,0,0,0,0,0,0,0,0,0,0,
    0,0,1,0,1,0,0,0,0,0,0,
    0,0,1,0,0,0,0,0,0,0,0,
    0,1,0,0,0,0,0,1,0,0,0,
    0,0,0,0,0,1,0,1,0,0,0,
    0,0,0,0,0,0,0,0,1,0,0,
    0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,1,1,0,0,
    0,0,0,0,0,0,0,1,1,0,0
  ), ncol = 11)
  colnames(A) <- nms; rownames(A) <- nms
  A
}

adj2dag <- function(adj){
  nodes <- rownames(adj)
  s <- "dag {"
  for (i in seq_len(nrow(adj))) for (j in seq_len(ncol(adj)))
    if (adj[i,j]==1) s <- paste(s, nodes[i], "->", nodes[j], ";")
  dagitty::dagitty(paste0(s, "}"))
}

# -------------------------------
# (A) Simulation-only API
# -------------------------------
# Returns: list(data, amat, weights)
# weights is a matrix W[parent, child] with linear SEM coefficients
# -------------------------------

simulate_dag_data <- function(amat, n, seed = 1,
                              w_mean = 0.8, w_sd = 0.2,
                              noise_sd = 1.0){
  stopifnot(all(rownames(amat) == colnames(amat)))
  set.seed(seed)
  vars <- rownames(amat)
  g <- igraph::graph_from_adjacency_matrix(amat, mode = "directed")
  stopifnot(igraph::is_dag(g))
  order <- igraph::as_ids(igraph::topo_sort(g, mode = "out"))
  
  E <- matrix(rnorm(n * length(vars), 0, noise_sd), nrow = n, ncol = length(vars))
  colnames(E) <- vars
  X <- matrix(0, n, length(vars)); colnames(X) <- vars
  W <- matrix(0, nrow = length(vars), ncol = length(vars), dimnames = list(vars, vars))
  
  for (v in order){
    parents <- vars[which(amat[, v] == 1L)]
    if (!length(parents)){
      X[, v] <- E[, v]
    } else {
      w  <- rnorm(length(parents), w_mean, w_sd)
      mu <- as.matrix(X[, parents, drop = FALSE]) %*% w
      X[, v] <- as.numeric(mu) + E[, v]
      W[parents, v] <- w
    }
  }
  
  list(
    data    = as.data.frame(X)[, rownames(amat), drop = FALSE],
    amat    = amat,
    weights = W
  )
}

# -------------------------------
# (B) CI enumeration and testing API
# - enumerate_all_cis(): brute-force minimal separating sets via dagitty::dseparated
# - run_ci_tests_from_list(): test arbitrary CI list via COMETs
# -------------------------------

# Enumerate all implied CIs (optionally empty set; minimal-only)
enumerate_all_cis <- function(amat,
                              include_empty = TRUE,
                              max_Z_size = NULL,
                              minimal_only = TRUE){
  g <- adj2dag(amat)
  nodes <- rownames(amat)
  p <- length(nodes)
  
  subsets <- function(vec, maxk = NULL){
    out <- list(character(0))
    if (length(vec)==0) return(out)
    for (k in 1:length(vec)){
      if (!is.null(maxk) && k > maxk) break
      out <- c(out, combn(vec, k, simplify = FALSE))
    }
    out
  }
  
  cis_list <- list()
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      X <- nodes[i]; Y <- nodes[j]
      rest <- setdiff(nodes, c(X, Y))
      allZ <- subsets(rest, maxk = max_Z_size)
      
      pair_sets <- list()
      for (Z in allZ){
        if (length(Z)==0 && !include_empty) next
        if (dagitty::dseparated(g, X = X, Y = Y, Z = Z)){
          pair_sets[[length(pair_sets)+1]] <- Z
        }
      }
      
      if (minimal_only && length(pair_sets) > 1){
        keep <- rep(TRUE, length(pair_sets))
        for (a in seq_along(pair_sets)){
          Za <- pair_sets[[a]]
          for (b in seq_along(pair_sets)){
            if (a==b) next
            Zb <- pair_sets[[b]]
            if (length(Zb) < length(Za) && all(Zb %in% Za)){
              keep[a] <- FALSE; break
            }
          }
        }
        pair_sets <- pair_sets[keep]
      }
      
      if (length(pair_sets)){
        for (Z in pair_sets){
          cis_list[[length(cis_list)+1]] <- list(X = X, Y = Y, Z = Z)
        }
      }
    }
  }
  cis_list
}

# Pretty CI table for manual browsing/picking
make_ci_table <- function(amat,
                          include_empty = TRUE,
                          max_Z_size = NULL,
                          minimal_only = TRUE){
  cis_list <- enumerate_all_cis(amat, include_empty, max_Z_size, minimal_only)
  tibble(
    ci_id = seq_along(cis_list),
    X     = vapply(cis_list, function(ci) ci$X, character(1)),
    Y     = vapply(cis_list, function(ci) ci$Y, character(1)),
    Z_str = vapply(cis_list, function(ci){
      if (length(ci$Z)) paste(sort(ci$Z), collapse = ", ") else "∅"
    }, character(1)),
    ci_obj = cis_list
  )
}

# -------------------------------
# COMETs wrapper (requires comets package installed)
# -------------------------------

need_pkg <- function(p) if (!requireNamespace(p, quietly = TRUE))
  stop(sprintf("Please install '%s': install.packages('%s')", p, p))
need_pkg("comets")

# General CI-test runner on a provided CI list
run_ci_tests_from_list <- function(cis_list, dat,
                                   tests = c("gcm","pcm"), alpha = 0.05){
  res <- lapply(tests, function(tst){
    pv <- vapply(seq_along(cis_list), function(k){
      ci <- cis_list[[k]]
      rhs <- paste0(
        paste0(ci$Y, collapse = "+"),
        " | ",
        if (length(ci$Z)) paste0(ci$Z, collapse = "+") else "1"
      )
      fm <- reformulate(rhs, response = ci$X)
      tryCatch(comets::comets(fm, dat, test = tst, coin = TRUE)$p.value,
               error = function(e) NA_real_)
    }, numeric(1))
    
    tibble(
      test = tst,
      CI = vapply(cis_list, function(ci){
        if (length(ci$Z))
          paste(ci$X, "_||_", paste(ci$Y, collapse = "+"), "|", paste(ci$Z, collapse = ", "))
        else
          paste(ci$X, "_||_", paste(ci$Y, collapse = "+"), "|", "∅")
      }, character(1)),
      p.value = pv,
      adj.p.value = p.adjust(pv, "holm"),
      rejected = p.adjust(pv, "holm") < alpha
    )
  }) |> dplyr::bind_rows()
  res
}

# Evaluate a DAG (amat) on a dataset (data): enumerate-all-CIs -> test
evaluate_dag_on_data <- function(amat, data,
                                 tests = c("gcm","pcm"), alpha = 0.05,
                                 include_empty = TRUE,
                                 max_Z_size = NULL,
                                 minimal_only = TRUE,
                                 tag = NULL){
  cis_list <- enumerate_all_cis(amat, include_empty, max_Z_size, minimal_only)
  res <- run_ci_tests_from_list(cis_list, data, tests, alpha)
  
  summary <- res %>%
    summarize(
      n_tests = n(),
      n_rejected = sum(rejected),
      min_adj_p = min(adj.p.value, na.rm = TRUE),
      median_adj_p = stats::median(adj.p.value, na.rm = TRUE)
    )
  if (!is.null(tag)) summary$tag <- tag
  
  list(results = res,
       rejected = dplyr::filter(res, rejected),
       summary = summary)
}

# -------------------------------
# Convenience helper: does a CI string involve a pair (a,b)?
# -------------------------------

involves_pair <- function(ci_str, a, b){
  parts <- strsplit(ci_str, " _\\|\\|_ ", perl = TRUE)[[1]]
  if (length(parts) < 2) return(FALSE)
  X <- trimws(parts[1])
  Y <- trimws(strsplit(parts[2], "\\|", perl = TRUE)[[1]][1])
  s <- sort(c(X, Y))
  identical(s, sort(c(a, b)))
}

# -------------------------------
# Example usage (uncomment to run)
# -------------------------------

# ## (1) SIMULATE once from (possibly modified) DAG
# amat_true <- get_sachs_amat()
# # Example: add Raf -> PIP3 to true generator
# amat_true["Raf","PIP3"] <- 1L
# sim <- simulate_dag_data(amat_true, n = 1000, seed = 777)
# # Save simulation set for reuse
# saveRDS(sim, file.path(out_dir, "sim_set_Raf_to_PIP3_n1000_seed777.rds"))
# write_csv(as_tibble(sim$data), file.path(out_dir, "sim_data.csv"))
# write_csv(as_tibble(as.data.frame(sim$weights), row_names = TRUE),
#           file.path(out_dir, "sim_weights_matrix.csv"))

# ## (2) EVALUATE any DAG on that data (e.g., original Sachs)
# amat_model <- get_sachs_amat()
# eval1 <- evaluate_dag_on_data(amat_model, sim$data, tests = c("gcm","pcm"), alpha = alpha,
#                               include_empty = TRUE, minimal_only = TRUE, tag = "Sachs on sim")
# print(eval1$summary)
# write_csv(eval1$results, file.path(out_dir, "eval_Sachs_on_sim_allCIs.csv"))
# write_csv(eval1$rejected, file.path(out_dir, "eval_Sachs_on_sim_rejected.csv"))

# ## (3) Focus on a manual pair and show its rejected CIs
# pair_X <- "Raf"; pair_Y <- "PIP3"
# pair_rej <- eval1$rejected %>%
#   dplyr::filter(vapply(CI, involves_pair, logical(1), a = pair_X, b = pair_Y)) %>%
#   arrange(adj.p.value)
# print(pair_rej)

# -------------------------------
# (C) Plotting helpers
# -------------------------------

suppressPackageStartupMessages({
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("ggraph", quietly = TRUE)
  requireNamespace("scales", quietly = TRUE)
  requireNamespace("tidyr", quietly = TRUE)
})

# Nice circular DAG plot with optional highlighted edge
plot_dag_circle <- function(amat,
                            order = c("PIP2","PLCg","Mek","Raf","JNK","p38","PKC","PKA","Akt","Erk","PIP3"),
                            anchor = "PIP2",
                            mirror_y = TRUE,
                            highlight_edge = NULL,   # c(from, to) or NULL
                            title = "Graph") {
  stopifnot(all(rownames(amat) %in% order), all(order %in% rownames(amat)))
  dag <- adj2dag(amat)
  k <- length(order)
  theta <- seq(0, 2*pi, length.out = k + 1)[1:k]; names(theta) <- order
  theta <- theta + ((pi/2) - theta[anchor])
  x <- cos(theta); y <- sin(theta); if (mirror_y) x <- -x
  pos <- data.frame(name = order, x = as.numeric(x), y = as.numeric(y), row.names = order)
  
  el <- dagitty::edges(dag)
  g  <- igraph::graph_from_data_frame(el, directed = TRUE, vertices = data.frame(name = order))
  
  p <- ggraph::ggraph(g, layout = "manual", x = pos$x, y = pos$y) +
    ggraph::geom_edge_link(arrow = grid::arrow(length = grid::unit(3, "mm"), type = "closed"),
                           end_cap = ggraph::circle(3.2, "mm"),
                           width = 0.6, alpha = 0.8) +
    ggraph::geom_node_point(size = 12, shape = 21, fill = "black", color = "black", stroke = 0) +
    ggraph::geom_node_text(ggplot2::aes(label = name), color = "white", size = 3.8, fontface = "bold") +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0, face = "bold"))
  
  if (!is.null(highlight_edge)) {
    he <- data.frame(from = highlight_edge[1], to = highlight_edge[2])
    # draw a bold segment on top
    p <- p + ggplot2::geom_segment(
      data = transform(merge(merge(he, pos, by.x = "from", by.y = "name"), pos, by.x = "to", by.y = "name"),
                       xend = x.y, yend = y.y),
      ggplot2::aes(x = x.x, y = y.x, xend = xend, yend = yend),
      inherit.aes = FALSE, linewidth = 1.6
    )
  }
  p
}

# Heatmap of edge weights (W[parent, child])
plot_weights_heatmap <- function(W, title = "Edge weights (true coefficients)") {
  df <- as.data.frame(W)
  df$parent <- rownames(df)
  df_long <- tidyr::pivot_longer(df, -parent, names_to = "child", values_to = "w")
  ggplot2::ggplot(df_long, ggplot2::aes(parent, child, fill = w)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title, x = "Parent", y = "Child", fill = "w") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

# Heatmap of adjacency matrix (0/1)
plot_adjacency_heatmap <- function(amat, title = "Adjacency (parent→child)") {
  df <- as.data.frame(amat)
  df$parent <- rownames(df)
  df_long <- tidyr::pivot_longer(df, -parent, names_to = "child", values_to = "a")
  ggplot2::ggplot(df_long, ggplot2::aes(parent, child, fill = factor(a))) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_manual(values = c("0" = "#f0f0f0", "1" = "#313695")) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title, x = "Parent", y = "Child", fill = "edge") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

# Scatter to visualize dependence for a selected pair (raw)
plot_pair_scatter <- function(dat, x, y, title = NULL) {
  if (is.null(title)) title <- paste0(x, " vs ", y)
  ggplot2::ggplot(dat, ggplot2::aes_string(x = x, y = y)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(title = title, x = x, y = y)
}

# Bar summary of rejections per test type
plot_rejection_summary <- function(results, title = "CI rejections by test") {
  stopifnot(all(c("test","rejected") %in% names(results)))
  agg <- results %>% group_by(test) %>% summarize(n_rejected = sum(rejected), .groups = "drop")
  ggplot2::ggplot(agg, ggplot2::aes(test, n_rejected)) +
    ggplot2::geom_col() +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(title = title, x = "Test", y = "# rejections (Holm)")
}

# ECDF of adjusted p-values per test
plot_pvalue_ecdf <- function(results, title = "ECDF of adjusted p-values") {
  ggplot2::ggplot(results, ggplot2::aes(adj.p.value, color = test)) +
    ggplot2::stat_ecdf(geom = "step") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(title = title, x = "adjusted p-value", y = "ECDF", color = "Test")
}

# -------------------------------
# Example plotting workflow (uncomment to run)
# -------------------------------
# ## Simulate data and plot DAG & weights
# amat_true <- get_sachs_amat(); amat_true["Raf","PIP3"] <- 1L
# sim <- simulate_dag_data(amat_true, n = 1000, seed = 777)
# p1 <- plot_dag_circle(amat_true, highlight_edge = c("Raf","PIP3"), title = "True DAG (Raf→PIP3)")
# p2 <- plot_adjacency_heatmap(amat_true)
# p3 <- plot_weights_heatmap(sim$weights)
# p4 <- plot_pair_scatter(sim$data, "Raf", "PIP3", title = "Dependence after adding Raf→PIP3")
# ggplot2::ggsave(file.path(out_dir, "dag_true_circle.png"), p1, width = 6, height = 6, dpi = 150)
# ggplot2::ggsave(file.path(out_dir, "adjacency_heatmap.png"), p2, width = 6, height = 5, dpi = 150)
# ggplot2::ggsave(file.path(out_dir, "weights_heatmap.png"), p3, width = 7, height = 6, dpi = 150)
# ggplot2::ggsave(file.path(out_dir, "raf_pip3_scatter.png"), p4, width = 6, height = 5, dpi = 150)
# 
# ## Evaluate Sachs DAG on the same data and visualize outcomes
# eval1 <- evaluate_dag_on_data(get_sachs_amat(), sim$data, tests = c("gcm","pcm"), alpha = 0.05,
#                               include_empty = TRUE, minimal_only = TRUE, tag = "Sachs on sim")
# p5 <- plot_rejection_summary(eval1$results)
# p6 <- plot_pvalue_ecdf(eval1$results)
# ggplot2::ggsave(file.path(out_dir, "rejections_by_test.png"), p5, width = 6, height = 4, dpi = 150)
# ggplot2::ggsave(file.path(out_dir, "adj_p_ecdf.png"), p6, width = 6, height = 4, dpi = 150)

# -------------------------------
# (D) Rejected-CI chords overlay on circular DAG
# -------------------------------

# Parse "X _||_ Y | Z" -> c(X, Y)
parse_ci_pair_xy <- function(ci_str){
  parts <- strsplit(ci_str, " _\\|\\|_ ", perl = TRUE)[[1]]
  parts
}

# Plot DAG and overlay dashed red chords for rejected CI pairs
plot_dag_with_rejected_cis <- function(amat, results,
                                       order = c("PIP2","PLCg","Mek","Raf","JNK","p38","PKC","PKA","Akt","Erk","PIP3"),
                                       anchor = "PIP2",
                                       mirror_y = TRUE,
                                       title = "Rejected CIs (red, dashed)"){
  stopifnot(all(c("CI","rejected") %in% names(results)))
  dag <- adj2dag(amat)
  # positions
  k <- length(order)
  theta <- seq(0, 2*pi, length.out = k + 1)[1:k]; names(theta) <- order
  theta <- theta + ((pi/2) - theta[anchor])
  x <- cos(theta); y <- sin(theta); if (mirror_y) x <- -x
  pos <- data.frame(name = order, x = as.numeric(x), y = as.numeric(y), row.names = order)
  
  # rejected CI pairs (unique, unordered)
  rej <- subset(results, isTRUE(rejected))
  if (nrow(rej) == 0){
    warning("No rejected CIs provided; plotting base DAG only.")
  }
  pairs <- do.call(rbind, lapply(rej$CI, function(s){
    xy <- parse_ci_pair_xy(s)
    if (any(is.na(xy))) return(NULL)
    data.frame(X = xy[1], Y = xy[2], stringsAsFactors = FALSE)
  }))
  if (!is.null(pairs) && nrow(pairs)){ 
    pairs <- as.data.frame(unique(t(apply(pairs, 1, function(v) sort(v)))))
    names(pairs) <- c("X","Y")
  } else {
    pairs <- data.frame(X=character(), Y=character())
  }
  
  # chords
  chords <- merge(pairs, transform(pos, name = rownames(pos)), by.x = "X", by.y = "name", all.x = TRUE)
  chords <- merge(chords, transform(pos, name = rownames(pos)), by.x = "Y", by_y = "name", all.x = TRUE, suffixes = c("", ".end"))
  chords <- chords[is.finite(chords$x) & is.finite(chords$y) & is.finite(chords$x.end) & is.finite(chords$y.end), , drop = FALSE]
  
  # base graph
  el <- dagitty::edges(dag)
  g  <- igraph::graph_from_data_frame(el, directed = TRUE, vertices = data.frame(name = order))
  p <- ggraph::ggraph(g, layout = "manual", x = pos$x, y = pos$y) +
    ggraph::geom_edge_link(arrow = grid::arrow(length = grid::unit(3, "mm"), type = "closed"),
                           end_cap = ggraph::circle(3.2, "mm"), width = 0.6, alpha = 0.8) +
    ggraph::geom_node_point(size = 12, shape = 21, fill = "black", color = "black", stroke = 0) +
    ggraph::geom_node_text(ggplot2::aes(label = name), color = "white", size = 3.8, fontface = "bold") +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0, face = "bold"))
  
  if (nrow(chords)){
    p <- p + ggplot2::geom_segment(
      data = chords,
      ggplot2::aes(x = x, y = y, xend = x.end, yend = y.end),
      inherit.aes = FALSE,
      linewidth = 0.9,
      linetype = "dashed",
      color = "red",
      alpha = 0.9
    )
  }
  p
}

# -------------------------------
# Example (uncomment)
# -------------------------------
# amat_true <- get_sachs_amat(); amat_true["Raf","PIP3"] <- 1L
# sim <- simulate_dag_data(amat_true, n = 1000, seed = 777)
# eval1 <- evaluate_dag_on_data(get_sachs_amat(), sim$data, tests = c("gcm","pcm"), alpha = 0.05,
#                               include_empty = TRUE, minimal_only = TRUE, tag = "Sachs on sim")
# p_rej <- plot_dag_with_rejected_cis(get_sachs_amat(), eval1$results,
#                                     title = "Sachs DAG — rejected CIs after Raf→PIP3 data")
# ggplot2::ggsave(file.path(out_dir, "dag_rejected_cis.png"), p_rej, width = 7, height = 7, dpi = 150)

# ============================================================
# End of file
# ============================================================
