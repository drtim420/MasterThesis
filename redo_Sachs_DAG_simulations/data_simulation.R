# ============================================================
# - Build Kook-style CI set: one CI per missing edge with the
#   smallest NON-EMPTY separating set Z (searched up to size 3)
# - Test with COMETs (GCM/PCM) + Holm
# - Simulate faithful data from the same DAG and compare
# ============================================================

set.seed(123)


data_path <- "../data"   
out_dir   <- "../results"
alpha     <- 0.05
tests     <- c("gcm","pcm")
max_k_Z   <- 3           # search Z up to size 3 for minimal non-empty d-sep


need_pkg <- function(p) if (!requireNamespace(p, quietly = TRUE))
  stop(sprintf("Please install '%s': install.packages('%s')", p, p))
for (p in c("readxl","dagitty","comets","igraph","MASS","dplyr","tibble","readr","utils"))
  need_pkg(p)

suppressPackageStartupMessages({
  library(readxl); library(dagitty); library(comets)
  library(igraph); library(MASS);   library(dplyr)
  library(tibble); library(readr)
})

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


if (requireNamespace("ranger", quietly = TRUE)) {
  
  if ("set_default_learner" %in% getNamespaceExports("comets")) {
    try(comets::set_default_learner("ranger"), silent = TRUE)
  }
  options(comets.default_learner = "ranger")
}


nms <- c("Raf","Mek","PLCg","PIP2","PIP3","Erk","Akt","PKA","PKC","p38","JNK")

# Sachs DAG 
get_sachs_amat <- function() {
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
  ), ncol = 11, byrow = TRUE)
  colnames(A) <- nms; rownames(A) <- nms
  A
}


adj2dag <- function(adj) {
  nodes <- rownames(adj); s <- "dag {"
  for (i in seq_len(nrow(adj))) for (j in seq_len(ncol(adj)))
    if (adj[i,j]==1) s <- paste(s, nodes[i], "->", nodes[j], ";")
  dagitty::dagitty(paste0(s, "}"))
}


read_sachs_observational <- function(path = "../data") {
  file <- file.path(path, "cd3cd28.xls")
  if (!file.exists(file)) stop("Missing file: ", file)
  dat <- readxl::read_xls(file)
  # log-transform numeric columns like lukas kook
  for (j in seq_along(dat)) if (is.numeric(dat[[j]])) dat[[j]] <- log(dat[[j]])
  colnames(dat) <- nms
  dat
}

# Find minimal NON-EMPTY separating set Z 
# Returns a character vector Z or NULL if none found within size <= max_k.
minimal_nonempty_sep <- function(dag, X, Y, universe, max_k = 3) {
  # Try |Z| = 1, then 2, ..., up to max_k
  for (k in 1:max_k) {
    if (length(universe) < k) break
    combs <- utils::combn(universe, k, simplify = FALSE)
    for (Z in combs) {
      if (dagitty::dseparated(dag, X, Y, Z)) return(Z)
    }
  }
  NULL
}


kook_ci_set <- function(amat, max_k = 3) {
  dag    <- adj2dag(amat)
  nodes  <- colnames(amat)
  pairs  <- utils::combn(nodes, 2, simplify = FALSE)
  out    <- list()
  for (pr in pairs) {
    a <- pr[1]; b <- pr[2]
    # skip if edge exists in either direction
    if (amat[a,b] == 1L || amat[b,a] == 1L) next
    Z <- minimal_nonempty_sep(dag, a, b, setdiff(nodes, c(a,b)), max_k = max_k)
    if (is.null(Z)) {
      # try the other direction of the query (d-sep is symmetric, but we keep the search symmetric for safety)
      Z <- minimal_nonempty_sep(dag, b, a, setdiff(nodes, c(a,b)), max_k = max_k)
    }
    if (!is.null(Z)) {
      # pick deterministic regression direction: alphabetic
      X <- sort(c(a,b))[1]
      Y <- setdiff(c(a,b), X)
      out[[length(out)+1]] <- list(X = X, Y = Y, Z = Z)
    }
  }
  # Return as tibble
  if (!length(out)) return(tibble(X=character(), Y=character(), Z_list=list(), CI=character()))
  tb <- tibble::tibble(
    X      = vapply(out, `[[`, character(1), "X"),
    Y      = vapply(out, `[[`, character(1), "Y"),
    Z_list = lapply(out, `[[`, "Z")
  )
  tb$CI <- mapply(function(x,y,z) paste0(x," _||_ ", y, " | ", paste(z, collapse = "+")),
                  tb$X, tb$Y, tb$Z_list, USE.NAMES = FALSE)
  tb
}


run_ci_tests_kook <- function(amat, dat, tests = c("gcm","pcm"), alpha = 0.05, max_k = 3) {
  cis_tbl <- kook_ci_set(amat, max_k = max_k)
  if (!nrow(cis_tbl)) {
    return(tibble(test=character(), X=character(), Y=character(), Z_list=list(),
                  CI=character(), p.value=double(), adj.p.value=double(), rejected=logical()))
  }
  run_one <- function(tst) {
    pv <- vapply(seq_len(nrow(cis_tbl)), function(i) {
      X <- cis_tbl$X[i]; Y <- cis_tbl$Y[i]; Z <- cis_tbl$Z_list[[i]]
      fm <- reformulate(paste0(paste(Y, collapse="+"), "|", paste(Z, collapse="+")), response = X)
      comets(fm, dat, test = tst, coin = TRUE)$p.value
    }, numeric(1))
    tibble(test = tst,
           X = cis_tbl$X, Y = cis_tbl$Y, Z_list = cis_tbl$Z_list, CI = cis_tbl$CI,
           p.value = pv,
           adj.p.value = p.adjust(pv, "holm"),
           rejected = p.adjust(pv, "holm") < alpha)
  }
  bind_rows(lapply(tests, run_one))
}

# Simulation linear SEM 
simulate_linear_from_dag <- function(amat, n, seed = 1,
                                     w_mean = 0.8, w_sd = 0.2,
                                     noise_sd = 1.0) {
  set.seed(seed)
  stopifnot(all(rownames(amat) == colnames(amat)))
  vars <- rownames(amat)
  g <- igraph::graph_from_adjacency_matrix(amat, mode="directed")
  stopifnot(igraph::is_dag(g))
  order <- igraph::as_ids(igraph::topo_sort(g, mode = "out"))
  
  # independent Gaussian noise
  E <- matrix(rnorm(n * length(vars), 0, noise_sd), nrow = n, ncol = length(vars))
  colnames(E) <- vars
  
  X <- matrix(0, n, length(vars)); colnames(X) <- vars
  for (v in order) {
    parents <- vars[which(amat[, v] == 1L)]
    if (!length(parents)) {
      X[, v] <- E[, v]
    } else {
      w  <- rnorm(length(parents), w_mean, w_sd)            # fixed weights per parent
      mu <- as.matrix(X[, parents, drop = FALSE]) %*% w     # linear mechanism
      X[, v] <- as.numeric(mu) + E[, v]
    }
  }
  as.data.frame(X)[, nms, drop = FALSE]
}





#  RUN 
cat("\n=== (A) REAL data: test Kook-style CI set ===\n")
amat     <- get_sachs_amat()
dat_real <- read_sachs_observational(data_path)

res_real <- run_ci_tests_kook(amat, dat_real, tests = tests, alpha = alpha, max_k = max_k_Z)
viol_real <- res_real %>%
  dplyr::filter(rejected) %>%
  dplyr::arrange(adj.p.value)


readr::write_csv(res_real, file.path(out_dir, "real_all_CIs_KOOK.csv"))
readr::write_csv(viol_real, file.path(out_dir, "real_rejected_CIs_KOOK.csv"))

cat("Total CIs tested (Kook set):", nrow(res_real), "\n")
cat("Rejected after Holm @ alpha =", alpha, ":", nrow(viol_real), "\n")
if (nrow(viol_real)) {
  cat("\nTop real-data rejections:\n")
  print(viol_real %>% select(test, CI, adj.p.value) %>% arrange(adj.p.value) %>% slice_head(n = 10))
} else {
  cat("\nNo CI was rejected on real data after correction (with this CI set).\n")
}

cat("\n=== (B) SIMULATE faithful data from the SAME DAG ===\n")
set.seed(999)
sim_dat <- simulate_linear_from_dag(
  amat, n = nrow(dat_real), seed = 999,
  w_mean = 0.8, w_sd = 0.2, noise_sd = 1.0
)

res_sim <- run_ci_tests_kook(amat, sim_dat, tests = tests, alpha = alpha, max_k = max_k_Z)
viol_sim <- res_sim %>%
  dplyr::filter(.data$rejected) %>%
  dplyr::arrange(.data$`adj.p.value`)


readr::write_csv(res_sim,  file.path(out_dir, "sim_all_CIs_KOOK.csv"))
readr::write_csv(viol_sim, file.path(out_dir, "sim_rejected_CIs_KOOK.csv"))

cat("Simulated: Total CIs tested (Kook set):", nrow(res_sim), "\n")
cat("Simulated: Rejected after Holm @ alpha =", alpha, ":", nrow(viol_sim), "\n")

cat("\n=== (C) SIDE-BY-SIDE comparison (same CI set, same tests) ===\n")
cmp <- dplyr::full_join(
  res_real %>% dplyr::select(test, CI, adj.p.value) %>% dplyr::rename(adjp_real = adj.p.value),
  res_sim  %>% dplyr::select(test, CI, adj.p.value) %>% dplyr::rename(adjp_sim  = adj.p.value),
  by = c("test","CI")
) %>%
  dplyr::mutate(
    rej_real = !is.na(adjp_real) & adjp_real < alpha,
    rej_sim  = !is.na(adjp_sim)  & adjp_sim  < alpha
  ) %>%
  dplyr::arrange(dplyr::desc(rej_real), adjp_real)


readr::write_csv(cmp, file.path(out_dir, "compare_real_vs_sim_all_CIs_KOOK.csv"))

only_real <- cmp %>%
  dplyr::filter(rej_real & !rej_sim)

readr::write_csv(only_real, file.path(out_dir, "compare_only_real_rejections_KOOK.csv"))

cat("CIs rejected on REAL but not on SIM (Kook set):", nrow(only_real), "\n")
if (nrow(only_real)) {
  print(only_real %>% select(test, CI, adjp_real, adjp_sim))
} else {
  cat("(None — try a different seed or larger n if you expect near-zero sim rejections.)\n")
}












######### consensus vs sachs DAG analysis



# Consensus DAG 
cons_amat <- matrix(c(
  0,0,0,0,0,0,0,1,1,0,0,
  1,0,0,0,0,0,0,1,1,0,0,
  0,0,0,0,1,0,0,0,0,0,0,
  0,0,1,0,1,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,
  0,1,0,0,0,0,0,1,0,0,0,
  0,0,0,0,1,0,0,1,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,
  0,0,1,1,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,1,1,0,0,
  0,0,0,0,0,0,0,1,1,0,0
), ncol = 11, byrow = TRUE)
dimnames(cons_amat) <- list(nms, nms)

# Sachs 
sachs_amat <- matrix(c(
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
), ncol = 11, byrow = TRUE)
dimnames(sachs_amat) <- list(nms, nms)


sachs_improved_amat <- sachs_amat


# Helpers 
edges_df <- function(amat) {
  which(amat == 1, arr.ind = TRUE) |>
    as.data.frame() |>
    transform(from = rownames(amat)[row], to = colnames(amat)[col]) |>
    dplyr::select(from, to)
}
undir_pair <- function(a, b) paste(sort(c(a, b)), collapse = " -- ")
get_pair_from_ci <- function(ci) {
  # 1) Keep the left side BEFORE the conditioning bar " | " (with spaces)
  left <- sub("\\s\\|\\s.*$", "", ci)  # e.g., "Akt _||_ Erk"
  
  # 2) Extract X and Y around the literal "_||_"
  m <- regexec("^\\s*(.*?)\\s*_\\|\\|_\\s*(.*?)\\s*$", left)
  parts <- regmatches(left, m)[[1]]
  if (length(parts) >= 3) {
    x <- trimws(parts[2]); y <- trimws(parts[3])
    return(undir_pair(x, y))  # "X -- Y"
  } else {
    return(NA_character_)
  }
}



E0 <- edges_df(sachs_amat)
E1 <- edges_df(sachs_improved_amat)
E0$key <- paste(E0$from, "->", E0$to)
E1$key <- paste(E1$from, "->", E1$to)

changed_keys <- base::setdiff(
  base::union(E0$key, E1$key),
  base::intersect(E0$key, E1$key)
)

changed_dir   <- dplyr::bind_rows(E0, E1) |>
  dplyr::filter(paste(from, "->", to) %in% changed_keys)
changed_pairs <- unique(undir_pair(changed_dir$from, changed_dir$to))

cat("\n[DIFF Sachs orig vs improved] #changed directed edges:", length(changed_keys),
    " | #changed unordered pairs:", length(changed_pairs), "\n")

# Run Kook/Lukas battery on CONSENSUS 
res_cons <- run_ci_tests_kook(cons_amat, dat_real, tests = tests, alpha = alpha, max_k = max_k_Z)
viol_cons <- res_cons |>
  dplyr::filter(.data$rejected) |>
  dplyr::arrange(.data$`adj.p.value`)
viol_cons$pair <- vapply(viol_cons$CI, get_pair_from_ci, character(1))

cat("[CONSENSUS] Total CIs:", nrow(res_cons),
    " | Rejected (Holm @", alpha, "):", nrow(viol_cons), "\n")

# Overlap: consensus rejections vs Sachs changed pairs 
viol_cons$matches_changed_edge <- viol_cons$pair %in% changed_pairs
matches <- viol_cons |>
  dplyr::filter(.data$matches_changed_edge) |>
  dplyr::select(test, CI, adj.p.value, pair)

if (nrow(matches)) {
  cat("\nMATCH FOUND: Rejected CI(s) in CONSENSUS involve the SAME variable pair as a changed edge in improved SACHS:\n")
  print(matches)
} else if (nrow(viol_cons)) {
  cat("\nNo overlap: rejected CI pairs in CONSENSUS do not coincide with changed Sachs edge pairs.\n")
  print(viol_cons |> dplyr::select(test, CI, adj.p.value) |> dplyr::slice_head(n = 5))
} else {
  cat("\nNo rejected CIs in CONSENSUS after Holm; nothing to compare.\n")
}

# Visuals 
plot_dag <- function(amat, title, highlight_pairs = character(0)) {
  g   <- igraph::graph_from_adjacency_matrix(amat, mode = "directed")
  lay <- igraph::layout_with_sugiyama(g)$layout
  vnm <- igraph::V(g)$name
  
  V(g)$color        <- "grey85"
  V(g)$frame.color  <- "grey40"
  E(g)$color        <- "grey55"
  E(g)$arrow.size   <- 0.35
  
  # highlight nodes that appear in any highlighted pair
  if (length(highlight_pairs)) {
    hp_nodes <- unique(unlist(strsplit(highlight_pairs, " -- ")))
    hit <- match(hp_nodes, vnm)
    V(g)$color[hit[!is.na(hit)]] <- "tomato"
  }
  
  plot(g, layout = lay, vertex.size = 22,
       vertex.label.color = "black", main = title)
  
  # dashed connectors between highlighted pairs
  if (length(highlight_pairs)) {
    for (p in highlight_pairs) {
      ab  <- strsplit(p, " -- ")[[1]]
      idx <- match(ab, vnm)           # <- use vertex names, not rownames(lay)
      if (any(is.na(idx)) || length(idx) != 2) next
      coords <- lay[idx, , drop = FALSE]
      if (nrow(coords) == 2) {
        segments(coords[1,1], coords[1,2], coords[2,1], coords[2,2], lty = 2, lwd = 2)
      }
    }
    legend("topleft",
           legend = c("Pair(s) of interest","Other nodes"),
           pch = 21, pt.bg = c("tomato","grey85"), bty = "n", pt.cex = 1.4)
  }
}


pairs_to_show <- unique(c(matches$pair, head(viol_cons$pair, 1)))  # show matches; fallback to top rejection

op <- par(mfrow = c(1,3), mar = c(1,1,3,1))
plot_dag(sachs_amat,          "Sachs (Original)",  highlight_pairs = changed_pairs)
plot_dag(cons_amat,           "Consensus (Rejected CI pairs)", highlight_pairs = pairs_to_show)
par(op)


run_ci_both_dirs(dat_real, left = "Akt", right = "Erk", Z = c("Raf"), tests = c("gcm","pcm"), alpha = 0.05)



# 1) After you build `viol_cons`, collect the rejected (unordered) pairs:
rej_pairs <- unique(viol_cons$pair)  # e.g., "Akt -- Erk", ...

# 2) Replace your plot_dag() with this version that can draw rejected pairs in red
plot_dag <- function(amat, title,
                     highlight_pairs = character(0),
                     rejected_pairs  = character(0),
                     draw_bidirectional = TRUE) {
  g   <- igraph::graph_from_adjacency_matrix(amat, mode = "directed")
  lay <- igraph::layout_with_sugiyama(g)$layout
  vnm <- igraph::V(g)$name
  
  V(g)$color       <- "grey85"
  V(g)$frame.color <- "grey40"
  E(g)$color       <- "grey55"
  E(g)$arrow.size  <- 0.35
  
  # Highlight nodes from highlight_pairs (tomato)
  if (length(highlight_pairs)) {
    hp_nodes <- unique(unlist(strsplit(highlight_pairs, " -- ")))
    idx <- match(hp_nodes, vnm)
    V(g)$color[idx[!is.na(idx)]] <- "tomato"
  }
  
  # Also highlight nodes that appear in rejected_pairs (tomato)
  if (length(rejected_pairs)) {
    rp_nodes <- unique(unlist(strsplit(rejected_pairs, " -- ")))
    idx <- match(rp_nodes, vnm)
    V(g)$color[idx[!is.na(idx)]] <- "tomato"
  }
  
  plot(g, layout = lay, vertex.size = 22,
       vertex.label.color = "black", main = title)
  
  # Draw dashed connectors for highlight_pairs (grey dashed)
  if (length(highlight_pairs)) {
    for (p in highlight_pairs) {
      ab  <- strsplit(p, " -- ")[[1]]
      idx <- match(ab, vnm)
      if (any(is.na(idx)) || length(idx) != 2) next
      coords <- lay[idx, , drop = FALSE]
      segments(coords[1,1], coords[1,2], coords[2,1], coords[2,2],
               lty = 2, lwd = 2)
    }
  }
  
  # Draw REJECTED CI pairs in RED (these are non-edges implied by Consensus but rejected by data)
  if (length(rejected_pairs)) {
    for (p in rejected_pairs) {
      ab  <- strsplit(p, " -- ")[[1]]
      idx <- match(ab, vnm)
      if (any(is.na(idx)) || length(idx) != 2) next
      coords <- lay[idx, , drop = FALSE]
      # thick red line
      segments(coords[1,1], coords[1,2], coords[2,1], coords[2,2],
               col = "red", lwd = 3)
      # optional: draw red arrows both ways to suggest a missing connection
      if (draw_bidirectional) {
        arrows(coords[1,1], coords[1,2], coords[2,1], coords[2,2],
               col = "red", lwd = 2, length = 0.08)
        arrows(coords[2,1], coords[2,2], coords[1,1], coords[1,2],
               col = "red", lwd = 2, length = 0.08)
      }
    }
    legend("topleft",
           legend = c("Rejected CI pair (non-edge)", "Other nodes"),
           pch = 21, pt.bg = c("tomato","grey85"),
           text.col = c("red","black"),
           col = c("red","grey40"), bty = "n", pt.cex = 1.4)
  }
}

# keep Sachs plots as-is, and pass rejected_pairs to the Consensus plot
pairs_to_show <- unique(c(matches$pair, head(viol_cons$pair, 1)))

op <- par(mfrow = c(1,3), mar = c(1,1,3,1))
plot_dag(sachs_amat,          "Sachs (Original)",  highlight_pairs = character(0))
plot_dag(cons_amat,           "Consensus (Rejected CIs)",
         rejected_pairs = rej_pairs, draw_bidirectional = TRUE)
par(op)

print(rej_pairs)  
print(nms)
