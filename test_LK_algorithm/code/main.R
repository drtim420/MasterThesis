# run_kook_on_data_folder.R
# Runs Kook/COMETs on ALL datasets found in ../data (csv/xlsx/xls).
# Data/DAG generation stays in data_simulation_with_DAG.R.

set.seed(1)
alpha <- 0.05
save  <- TRUE
data_path <- "../data"
tests <- c("gcm","pcm")

library(comets)
library(tidyverse)
library(readxl)
library(dagitty)
library(tidyr)

# bring node list + helpers (NO generation here)
source("data_simulation_with_DAG.R")   # provides nms, get_dag(), adj2dag(), etc.
use_saved_dag <- FALSE   # <- FALSE = test ALL datasets against get_dag() (the original)

# ---- discover all datasets (csv/xlsx/xls) ----
files_csv  <- list.files(data_path, pattern = "\\.csv$",  full.names = TRUE, ignore.case = TRUE)
files_xlsx <- list.files(data_path, pattern = "\\.xlsx$", full.names = TRUE, ignore.case = TRUE)
files_xls  <- list.files(data_path, pattern = "\\.xls$",  full.names = TRUE, ignore.case = TRUE)
files <- unique(c(files_csv, files_xlsx, files_xls))
if (length(files) == 0) stop("No datasets found in ", data_path)

# ---- loader: reads data + finds matching DAG (basename_dagitty.txt) ----
load_dataset_with_dag <- function(file_path) {
  ext <- tools::file_ext(file_path) |> tolower()
  dat <- switch(ext,
                "csv"  = readr::read_csv(file_path, show_col_types = FALSE),
                "xlsx" = readxl::read_xlsx(file_path),
                "xls"  = readxl::read_xls(file_path),
                stop("Unsupported extension: ", ext)
  )
  
  # Keep only variables we know; reorder to nms if possible
  if (!all(nms %in% colnames(dat))) {
    warning("Skipping (missing required columns): ", basename(file_path))
    return(NULL)
  }
  dat <- dat[, nms, drop = FALSE]
  colnames(dat) <- nms
  
  base <- tools::file_path_sans_ext(basename(file_path))
  dir  <- dirname(file_path)
  dag_file <- file.path(dir, paste0(base, "_dagitty.txt"))
  
  # choose hypothesis DAG
  if (use_saved_dag) {
    if (file.exists(dag_file)) {
      dag <- dagitty::dagitty(paste(readLines(dag_file), collapse = " "))
      amat <- matrix(0, length(nms), length(nms), dimnames = list(nms, nms))
      eds  <- dagitty::edges(dag)
      if (nrow(eds)) for (k in seq_len(nrow(eds))) amat[eds[k,1], eds[k,2]] <- 1
    } else {
      amat <- get_dag()  # fallback to original if no saved DAG
    }
  } else {
    amat <- get_dag()    # <-- ORIGINAL DAG as hypothesis for ALL datasets
  }
  list(data = dat, amat = amat, dataset = basename(file_path))
  
}

# ---- run COMETs on one dataset ----
run_ci_tests <- function(amat, dat, tests = c("gcm","pcm")) {
  dag  <- adj2dag(amat)
  cis  <- impliedConditionalIndependencies(dag)
  tcis <- cis[vapply(cis, function(x) length(x$Z) > 0, logical(1))]
  
  run_one <- function(tst) {
    pv <- vapply(tcis, function(ci) {
      fm <- reformulate(
        paste0(paste0(ci$Y, collapse = "+"), "|", paste0(ci$Z, collapse = "+")),
        response = ci$X
      )
      comets(fm, dat, test = tst, coin = TRUE)$p.value
    }, numeric(1))
    tibble(test = tst, CI = sapply(tcis, paste), p.value = pv,
           adj.p.value = p.adjust(pv, "holm"))
  }
  bind_rows(lapply(tests, run_one))
}

# ---- batch runner over ALL files in ../data ----
batched <- list()
for (f in files) {
  obj <- try(load_dataset_with_dag(f), silent = TRUE)
  if (inherits(obj, "try-error") || is.null(obj)) next
  res <- run_ci_tests(obj$amat, obj$data, tests = tests) %>%
    mutate(dataset = obj$dataset)
  batched[[length(batched)+1]] <- res
}
if (length(batched) == 0) stop("No valid datasets with required columns found.")

res <- bind_rows(batched)

# ---- summaries & reporting ----
summary_tbl <- res %>%
  mutate(viol = adj.p.value < alpha) %>%
  group_by(dataset, test) %>%
  summarise(n_CI = n(), n_viol = sum(viol), prop_viol = mean(viol), .groups = "drop")

cat("Datasets tested:", length(unique(res$dataset)), "\n")
print(summary_tbl %>% arrange(dataset, test))

viol <- res %>% filter(adj.p.value < alpha) %>% arrange(adj.p.value)
if (nrow(viol) == 0) {
  message("No CI violations at Holm-adjusted alpha = ", alpha, ".")
} else {
  print(viol, n = Inf)
}

# optional: top-10 smallest raw p per test
res %>% arrange(p.value) %>% group_by(test) %>% slice_head(n = 10) %>% ungroup() %>% print(n = Inf)

# ---- save outputs ----
if (save) {
  dir.create("../results", showWarnings = FALSE)
  readr::write_csv(res, "../results/all-tests_multi.csv")
  readr::write_csv(summary_tbl, "../results/summary_multi.csv")
  if (nrow(viol) > 0) readr::write_csv(viol, "../results/violations_multi.csv")
}








# new controlled approach

# ============================
# Changed-edge detection (independent block)
# ============================

message("\n=== Changed-edge detection block (independent) ===")
dir.create("../results", showWarnings = FALSE)

# -- helpers (local to this block) --

.read_dataset_only <- function(file_path) {
  ext <- tolower(tools::file_ext(file_path))
  dat <- switch(ext,
                "csv"  = readr::read_csv(file_path, show_col_types = FALSE),
                "xlsx" = readxl::read_xlsx(file_path),
                "xls"  = readxl::read_xls(file_path),
                stop("Unsupported extension: ", ext)
  )
  if (!all(nms %in% colnames(dat))) return(NULL)
  dat <- dat[, nms, drop = FALSE]; colnames(dat) <- nms
  dat
}

.adj_from_dagitty_file <- function(path, nms) {
  if (!file.exists(path)) return(NULL)
  dag  <- dagitty::dagitty(paste(readLines(path), collapse = " "))
  amat <- matrix(0, length(nms), length(nms), dimnames = list(nms, nms))
  eds  <- dagitty::edges(dag)
  if (nrow(eds)) for (k in seq_len(nrow(eds))) amat[eds[k,1], eds[k,2]] <- 1
  amat
}

.edges_df <- function(amat) {
  el <- which(amat == 1, arr.ind = TRUE)
  if (!nrow(el)) return(tibble::tibble(from=character(), to=character()))
  tibble::tibble(from = rownames(amat)[el[,1]], to = colnames(amat)[el[,2]])
}

.diff_edges <- function(before, after) {
  added   <- dplyr::anti_join(after,  before, by = c("from","to"))
  removed <- dplyr::anti_join(before, after,  by = c("from","to"))
  flips <- dplyr::inner_join(
    dplyr::rename(added,  from2=from, to2=to),
    dplyr::rename(removed,from2=to,   to2=from),
    by = c("from2","to2")
  )
  flipped <- if (nrow(flips)) dplyr::transmute(flips, from = to2, to = from2) else added[0,]
  added   <- dplyr::anti_join(added,   flipped, by=c("from","to"))
  removed <- dplyr::anti_join(removed, flipped, by=c("from","to"))
  list(added=added, removed=removed, flipped=flipped)
}

# FIXED: no dagitty.ci objects in tibble (store only X/Y/Z/CI)
.run_ci_tests_rich <- function(amat, dat, tests = c("gcm","pcm")) {
  dag  <- adj2dag(amat)
  cis  <- impliedConditionalIndependencies(dag)
  tcis <- cis[vapply(cis, function(x) length(x$Z) > 0, logical(1))]
  
  make_ci_tbl <- function() tibble::tibble(
    X      = vapply(tcis, function(ci) ci$X, character(1)),
    Y_list = lapply(tcis, function(ci) ci$Y),
    Z_list = lapply(tcis, function(ci) ci$Z),
    CI     = vapply(tcis, paste, character(1))
  )
  
  run_one <- function(tst) {
    pv <- vapply(tcis, function(ci) {
      fm <- reformulate(
        paste0(paste0(ci$Y, collapse = "+"), "|", paste0(ci$Z, collapse = "+")),
        response = ci$X
      )
      comets(fm, dat, test = tst, coin = TRUE)$p.value
    }, numeric(1))
    out <- make_ci_tbl()
    out$test <- tst
    out$p.value <- pv
    out$adj.p.value <- p.adjust(pv, "holm")
    out
  }
  
  dplyr::bind_rows(lapply(tests, run_one))
}

.edge_detection_scores <- function(amat_true, amat_alt, res_one, alpha = 0.05) {
  chg <- .diff_edges(.edges_df(amat_true), .edges_df(amat_alt))
  changed <- dplyr::bind_rows(
    dplyr::mutate(chg$flipped, type="flipped"),
    dplyr::mutate(chg$added,   type="added"),
    dplyr::mutate(chg$removed, type="removed")
  )
  if (!nrow(changed)) return(dplyr::mutate(changed, score=0L, n_related=0L, edge_label=character()))
  
  viol <- res_one %>% dplyr::filter(adj.p.value < alpha)
  
  ci_mentions <- function(ci_row, nodes) {
    x_hit <- ci_row$X %in% nodes
    y_hit <- any(unlist(ci_row$Y_list) %in% nodes)
    x_hit || y_hit
  }
  
  changed %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      edge_label = paste0(from,"→",to),
      n_related  = sum(vapply(seq_len(nrow(viol)), function(i)
        ci_mentions(viol[i,], c(from,to)), logical(1))),
      score = n_related
    ) %>% dplyr::ungroup()
}

.plot_edge_detection_graph <- function(amat_true, scores, title="Changed edges (colored by detection score)") {
  g <- igraph::graph_from_adjacency_matrix(amat_true, mode="directed")
  igraph::E(g)$color <- "grey80"; igraph::E(g)$lwd <- 1
  if (nrow(scores)) {
    for (k in seq_len(nrow(scores))) {
      f <- scores$from[k]; t <- scores$to[k]
      idx1 <- igraph::get.edge.ids(g, c(f, t))
      idx2 <- igraph::get.edge.ids(g, c(t, f))
      col <- if (scores$score[k] > 0) "#d62728" else "#ff9896"
      lw  <- 1 + 2*log1p(scores$score[k])
      if (idx1 != 0) { igraph::E(g)$color[idx1] <- col; igraph::E(g)$lwd[idx1] <- lw }
      if (idx2 != 0) { igraph::E(g)$color[idx2] <- col; igraph::E(g)$lwd[idx2] <- lw }
    }
  }
  plot(g,
       main = title,
       vertex.color = "white",
       vertex.frame.color = "grey30",
       vertex.label.color = "black",
       edge.arrow.size = 0.3)
}

.plot_edge_detection_bars <- function(scores, title="Detection score per changed edge") {
  if (!nrow(scores)) { message("No changed edges."); return(NULL) }
  df <- scores %>% dplyr::mutate(edge = edge_label)
  gg <- ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(edge, score), y = score, fill = type)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "changed edge", y = "detection score\n(# rejected CIs mentioning endpoints)",
                  title = title) +
    ggplot2::theme_minimal()
  gg
}

# -- Discover datasets again, but ignore logs (edges/added/removed/flipped/CI files) --
files2_csv  <- list.files(data_path, pattern = "\\.csv$",  full.names = TRUE, ignore.case = TRUE)
files2_xlsx <- list.files(data_path, pattern = "\\.xlsx$", full.names = TRUE, ignore.case = TRUE)
files2_xls  <- list.files(data_path, pattern = "\\.xls$",  full.names = TRUE, ignore.case = TRUE)
files2 <- unique(c(files2_csv, files2_xlsx, files2_xls))
bad_pat <- "(_edges|_added|_removed|_flipped|_ci_only_in_)"
files2 <- files2[!grepl(bad_pat, basename(files2), ignore.case = TRUE)]
if (!length(files2)) stop("Changed-edge block: no datasets found in ", data_path)

# -- Run detection per dataset --
for (f in files2) {
  dat <- try(.read_dataset_only(f), silent = TRUE)
  if (inherits(dat, "try-error") || is.null(dat)) {
    message("Skipping (missing required columns): ", basename(f))
    next
  }
  ds_name  <- basename(f)
  ds_dir   <- dirname(f)
  base_only <- tools::file_path_sans_ext(ds_name)
  
  # TRUE & PERTURBED DAGs: look for files your simulation block created
  true_txt <- list.files(ds_dir, pattern = "_TRUE_dagitty\\.txt$",      full.names = TRUE, ignore.case = TRUE)
  alt_txt  <- list.files(ds_dir, pattern = "_PERTURBED_dagitty\\.txt$", full.names = TRUE, ignore.case = TRUE)
  sidecar  <- file.path(ds_dir, paste0(base_only, "_dagitty.txt"))
  
  amat_true <- if (length(true_txt)) .adj_from_dagitty_file(true_txt[1], nms) else NULL
  if (is.null(amat_true)) amat_true <- get_dag()  # fallback to original DAG
  
  amat_alt  <- if (length(alt_txt))  .adj_from_dagitty_file(alt_txt[1],  nms) else NULL
  if (is.null(amat_alt) && file.exists(sidecar)) amat_alt <- .adj_from_dagitty_file(sidecar, nms)
  
  if (is.null(amat_alt)) {
    message("No perturbed DAG found for ", ds_name, " (no _PERTURBED_dagitty.txt or sidecar). Skipping plots.")
    next
  }
  
  # Hypothesis = TRUE/original DAG; data = from perturbed DAG
  res_rich <- .run_ci_tests_rich(amat_true, dat, tests = tests)
  
  # Score changed edges
  scores <- .edge_detection_scores(amat_true, amat_alt, res_rich, alpha = alpha)
  
  # Save plots
  png(file.path("../results", paste0(base_only, "_edge_detection_graph.png")), width = 1200, height = 900)
  .plot_edge_detection_graph(amat_true, scores,
                             title = paste0(base_only, ": changed edges colored by detection score"))
  dev.off()
  
  p <- .plot_edge_detection_bars(scores,
                                 title = paste0(base_only, ": detection score per changed edge"))
  if (!is.null(p)) {
    ggplot2::ggsave(filename = file.path("../results", paste0(base_only, "_edge_detection_bars.png")),
                    plot = p, width = 10, height = 6, dpi = 150)
  }
  
  # Save a CSV with the scores
  readr::write_csv(scores, file.path("../results", paste0(base_only, "_edge_detection_scores.csv")))
  
  message("Detection plots + scores saved for: ", ds_name)
}

message("=== Changed-edge detection block finished ===\n")
