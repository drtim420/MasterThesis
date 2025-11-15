# ============================================================
# 01_test_vs_sim_sachs.R
# Test DAG vs REAL Sachs data: list rejected implied CIs
# Simulate faithful data from the DAG and re-test same CIs
#    
# ============================================================

set.seed(123)


data_path <- "../data"   # must contain cd3cd28.xls
out_dir   <- "../results"
alpha     <- 0.05
tests     <- c("gcm","pcm")

# Dependencies 
need_pkg <- function(p) if (!requireNamespace(p, quietly = TRUE))
  stop(sprintf("Please install '%s': install.packages('%s')", p, p))
for (p in c("readxl","dagitty","comets","igraph","dplyr","tibble","readr"))
  need_pkg(p)

suppressPackageStartupMessages({
  library(readxl); library(dagitty); library(comets)
  library(igraph); library(dplyr);   library(tibble); library(readr)
})

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Sachs nodes taken from the Paper implementation of lukas kook
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
  ), ncol = 11)
  colnames(A) <- nms; rownames(A) <- nms
  A
}


adj2dag <- function(adj) {
  nodes <- rownames(adj); s <- "dag {"
  for (i in seq_len(nrow(adj))) for (j in seq_len(ncol(adj)))
    if (adj[i,j]==1) s <- paste(s, nodes[i], "->", nodes[j], ";")
  dagitty::dagitty(paste0(s, "}"))
}




################################################################################
## now we plot the DAG, reproducing the DAG plot style from the lukas kook paper
################################################################################
library(ggdag)
# Paper-style circle with anchor at 12 and optional y-axis mirror
plot_sachs_paper_circle <- function(
    dag,
    order  = c("PIP2","PLCg","Mek","Raf","JNK","p38","PKC","PKA","Akt","Erk","PIP3"),
    anchor = "PIP2",
    mirror_y = TRUE,        
    title  = "Sachs graph"
){
  stopifnot(setequal(order, names(dag)))
  stopifnot(anchor %in% order)
  
  # angles and rotation so anchor is at 12 o'clock
  k <- length(order)
  theta <- seq(0, 2*pi, length.out = k + 1)[1:k]
  names(theta) <- order
  shift <- (pi/2) - theta[anchor]     # rotate so anchor -> top
  theta <- theta + as.numeric(shift)
  
  # base circle coords
  x <- cos(theta); y <- sin(theta)
  if (isTRUE(mirror_y)) x <- -x       # mirror across y-axis (left-right flip)
  
  pos <- data.frame(name = order, x = x, y = y, row.names = order)
  
  # nice plot with ggraph if present; else base dagitty
  if (requireNamespace("ggraph", quietly = TRUE) &&
      requireNamespace("igraph", quietly = TRUE) &&
      requireNamespace("ggplot2", quietly = TRUE)) {
    
    el <- dagitty::edges(dag)
    g  <- igraph::graph_from_data_frame(el, directed = TRUE,
                                        vertices = data.frame(name = order))
    p <- ggraph::ggraph(g, layout = "manual", x = pos$x, y = pos$y) +
      ggraph::geom_edge_link(
        arrow   = grid::arrow(length = grid::unit(3, "mm"), type = "closed"),
        end_cap = ggraph::circle(3.2, "mm"),
        width   = 0.6
      ) +
      ggraph::geom_node_point(size = 12, shape = 21, fill = "black", color = "black", stroke = 0) +
      ggraph::geom_node_text(ggplot2::aes(label = name), color = "white", size = 3.8, fontface = "bold") +
      ggplot2::coord_fixed() +
      ggplot2::labs(title = title) +
      ggplot2::theme_void(base_size = 12) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0, face = "bold"))
    print(p); invisible(p)
    
  } else {
    dagitty::coordinates(dag) <- pos[, c("x","y")]
    op <- par(mar = c(1,1,2,1)); on.exit(par(op), add = TRUE)
    plot(dag, main = title)
    invisible(NULL)
  }
}


amat <- get_sachs_amat()
dag  <- adj2dag(amat)

paper_order <- c("PIP2","PLCg","Mek","Raf","JNK","p38","PKC","PKA","Akt","Erk","PIP3")
plot_sachs_paper_circle(dag, order = paper_order, anchor = "PIP2", mirror_y = TRUE,
                        title = "Sachs graph")

################################################################################



# find the conditional independencies 

cis_with_Z <- function(amat) {
  g <- adj2dag(amat)
  cis <- dagitty::impliedConditionalIndependencies(g)
  cis[vapply(cis, function(x) length(x$Z) > 0, logical(1))]
}


dag <- adj2dag(get_sachs_amat())    
cis <- impliedConditionalIndependencies(dag)
length(Filter(function(ci) length(ci$Z) > 0, cis))  # expect 22 




################################################################################
# plot the CI's
################################################################################
# Robust: set circular coords as a LIST with named numeric x,y vectors
ensure_circular_coords <- function(dag,
                                   order,
                                   anchor = order[1],
                                   mirror_y = TRUE) {
  nodes <- names(dag)
  if (!setequal(order, nodes)) {
    stop("`order` must be a permutation of these nodes: ",
         paste(sort(nodes), collapse = ", "))
  }
  if (!(anchor %in% order)) stop("`anchor` must be in `order`.")
  
  k <- length(order)
  theta <- seq(0, 2*pi, length.out = k + 1)[1:k]
  names(theta) <- order
  
  # rotate so `anchor` is at 12 o'clock
  shift <- (pi/2) - theta[anchor]
  theta <- theta + as.numeric(shift)
  
  x <- cos(theta); y <- sin(theta)
  if (mirror_y) x <- -x
  
  # Reorder x,y to dag's node order and ensure numeric + names
  xv <- setNames(as.numeric(x[nodes]), nodes)
  yv <- setNames(as.numeric(y[nodes]), nodes)
  
  if (!all(is.finite(xv)) || !all(is.finite(yv))) stop("Non-finite coords computed.")
  
  # >>> Key: assign as a LIST with named vectors
  dagitty::coordinates(dag) <- list(x = xv, y = yv)
  dag
}

# now we want to highlight the CI pairs we test
# Plot DAG + overlay tested CI pairs as red dashed chords
plot_dag_with_tested_cis <- function(dag, cis_list,
                                     title = "Sachs DAG — tested CIs (red, dashed)",
                                     node_size = 11) {
  # Unique unordered tested pairs
  pair_key <- vapply(cis_list, function(ci) paste(sort(c(ci$X, ci$Y)), collapse = "||"), "")
  pair_key <- unique(pair_key)
  pairs <- do.call(rbind, lapply(pair_key, function(k) {
    uv <- strsplit(k, "\\|\\|")[[1]]
    data.frame(X = uv[1], Y = uv[2], stringsAsFactors = FALSE)
  }))
  
  pos <- as.data.frame(dagitty::coordinates(dag))
  pos$name <- rownames(pos)
  
  # chords df
  chords <- merge(pairs, pos, by.x = "X", by.y = "name", all.x = TRUE)
  chords <- merge(chords, pos, by.x = "Y", by.y = "name", all.x = TRUE,
                  suffixes = c("", ".end"))
  
  # sanity: drop any rows with missing coords (should be none now)
  chords <- chords[is.finite(chords$x) & is.finite(chords$y) &
                     is.finite(chords$x.end) & is.finite(chords$y.end), , drop = FALSE]
  
  if (requireNamespace("ggraph", quietly = TRUE) &&
      requireNamespace("igraph", quietly = TRUE) &&
      requireNamespace("ggplot2", quietly = TRUE)) {
    
    el <- dagitty::edges(dag)
    vdf <- data.frame(name = pos$name, stringsAsFactors = FALSE)
    g  <- igraph::graph_from_data_frame(el, directed = TRUE, vertices = vdf)
    
    layout_df <- pos[match(igraph::V(g)$name, pos$name), c("x","y")]
    layout_df <- as.data.frame(lapply(layout_df, as.numeric))
    
    p <- ggraph::ggraph(g, layout = layout_df) +
      ggraph::geom_edge_link(
        arrow   = grid::arrow(length = grid::unit(3, "mm"), type = "closed"),
        end_cap = ggraph::circle(3.2, "mm"),
        width   = 0.6
      ) +
      ggraph::geom_node_point(size = node_size, shape = 21,
                              fill = "black", color = "black", stroke = 0) +
      ggraph::geom_node_text(ggplot2::aes(label = name),
                             color = "white", size = 3.8, fontface = "bold") +
      ggplot2::geom_segment(
        data = chords,
        ggplot2::aes(x = x, y = y, xend = x.end, yend = y.end),
        inherit.aes = FALSE,
        linetype = "dashed",
        linewidth = 0.8,
        color = "red"
      ) +
      ggplot2::coord_fixed(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), expand = TRUE) +
      ggplot2::labs(title = title) +
      ggplot2::theme_void(base_size = 12) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0, face = "bold"))
    
    print(p); invisible(p)
    
  } else {
    op <- par(mar = c(1,1,2,1)); on.exit(par(op), add = TRUE)
    plot(dag, main = title, cex = 1.1)
    segments(chords$x, chords$y, chords$x.end, chords$y.end,
             lty = 2, lwd = 2, col = "red")
    invisible(NULL)
  }
}


amat <- get_sachs_amat()
dag  <- adj2dag(amat)

paper_order <- c("PIP2","PLCg","Mek","Raf","JNK","p38","PKC","PKA","Akt","Erk","PIP3")
dag_circ <- ensure_circular_coords(dag, order = paper_order, anchor = "PIP2", mirror_y = TRUE)

str(dagitty::coordinates(dag_circ))  # should show numeric x,y (not logical NA)
# 'data.frame': 11 obs. of  2 variables:
#  $ x: num ...
#  $ y: num ...

# Now overlay the CI chords:
tcis <- cis_with_Z(amat)
plot_dag_with_tested_cis(dag_circ, tcis,
                         title = "Sachs DAG — tested CI pairs (red, dashed)")


# Set robust circular coords (PIP2 at 12 o’clock; mirror to match your figure)
dag_circ <- ensure_circular_coords(dag, order = paper_order, anchor = "PIP2", mirror_y = TRUE)

# Get the CI list you test
tcis <- cis_with_Z(amat)

# Plot with red dashed CI chords
plot_dag_with_tested_cis(dag_circ, tcis,
                         title = "Sachs DAG — tested CI pairs (red, dashed)")


str(dagitty::coordinates(dag_circ))   # should be data.frame with x,y for all nodes

## end of plots
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 



# Read original dataset 
read_sachs_observational <- function(path = "../data") {
  file <- file.path(path, "cd3cd28.xls")
  if (!file.exists(file)) stop("Missing file: ", file)
  dat <- readxl::read_xls(file)
  # log-transform numeric columns (Sachs convention)
  for (j in seq_along(dat)) if (is.numeric(dat[[j]])) dat[[j]] <- log(dat[[j]])
  colnames(dat) <- nms
  dat
}

# ------------- CI testing ---------------
run_ci_tests <- function(amat, dat, tests = c("gcm", "pcm"), alpha = 0.05) {
  tcis <- cis_with_Z(amat)
  
  res <- lapply(tests, function(tst) {
    pv <- vapply(seq_along(tcis), function(k) {
      ci <- tcis[[k]]
      fm <- reformulate(
        paste0(paste0(ci$Y, collapse = "+"), "|", paste0(ci$Z, collapse = "+")),
        response = ci$X
      )
      comets(fm, dat, test = tst, coin = TRUE)$p.value
    }, numeric(1))
    
    tibble(
      test = tst,
      CI = vapply(tcis, paste, character(1)),
      p.value = pv,
      adj.p.value = p.adjust(pv, "holm"),
      rejected = p.adjust(pv, "holm") < alpha
    )
  }) |> dplyr::bind_rows()
  
  return(res)
}

## now we are simulating data based on the DAG given 
####################################################
# Simulation: linear 
simulate_linear_from_dag <- function(amat, n, seed = 1,
                                     w_mean = 0.8, w_sd = 0.2,
                                     noise_sd = 1.0) {
  # amat = 
  # n = 10

  set.seed(seed)
  stopifnot(all(rownames(amat) == colnames(amat)))
  vars <- rownames(amat)
  g <- igraph::graph_from_adjacency_matrix(amat, mode = "directed")
  stopifnot(igraph::is_dag(g))
  order <- igraph::as_ids(igraph::topo_sort(g, mode = "out"))
  
  E <- matrix(rnorm(n * length(vars), 0, noise_sd), nrow = n, ncol = length(vars))
  colnames(E) <- vars
  
  X <- matrix(0, n, length(vars)); colnames(X) <- vars
  for (v in order) {
    parents <- vars[which(amat[, v] == 1L)]
    if (!length(parents)) {
      X[, v] <- E[, v]
    } else {
      w  <- rnorm(length(parents), w_mean, w_sd)              # fixed weights per parent
      mu <- as.matrix(X[, parents, drop = FALSE]) %*% w       # linear mechanism
      X[, v] <- as.numeric(mu) + E[, v]
    }
  }
  as.data.frame(X)[, nms, drop = FALSE]
}





# now after defining all the functions we can run the experiments
# RUN 

cat("\n=== (A) REAL data: test all implied CIs of the DAG ===\n")
amat     <- get_sachs_amat()
dat_real <- read_sachs_observational(data_path)

res_real <- run_ci_tests(amat, dat_real, tests = tests, alpha = alpha)
viol_real <- res_real %>% dplyr::filter(rejected) %>% dplyr::arrange(adj.p.value)

readr::write_csv(res_real, file.path(out_dir, "real_all_CIs.csv"))
readr::write_csv(viol_real, file.path(out_dir, "real_rejected_CIs.csv"))

cat("Total CIs tested:", nrow(res_real), "\n")
cat("Rejected after Holm @ alpha =", alpha, ":", nrow(viol_real), "\n")
if (nrow(viol_real)) {
  cat("\nTop real-data rejections:\n")
  print(
    viol_real %>%
      dplyr::select(test, CI, adj.p.value) %>%
      dplyr::arrange(adj.p.value) %>%
      dplyr::slice_head(n = 10)
  )
} else {
  cat("\nNo implied CI was rejected on real data after correction.\n")
}

cat("\n=== (B) SIMULATE faithful data from the SAME DAG ===\n")
set.seed(999)
sim_dat <- simulate_linear_from_dag(
  amat, n = nrow(dat_real), seed = 999,
  w_mean = 0.8, w_sd = 0.2, noise_sd = 1.0
)

res_sim <- run_ci_tests(amat, sim_dat, tests = tests, alpha = alpha)
viol_sim <- res_sim %>% dplyr::filter(rejected) %>% dplyr::arrange(adj.p.value)

readr::write_csv(res_sim,  file.path(out_dir, "sim_all_CIs.csv"))
readr::write_csv(viol_sim, file.path(out_dir, "sim_rejected_CIs.csv"))

cat("Simulated: Total CIs tested:", nrow(res_sim), "\n")
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

readr::write_csv(cmp, file.path(out_dir, "compare_real_vs_sim_all_CIs.csv"))

only_real <- cmp %>% dplyr::filter(rej_real & !rej_sim)
readr::write_csv(only_real, file.path(out_dir, "compare_only_real_rejections.csv"))

cat("CIs rejected on REAL but not on SIM:", nrow(only_real), "\n")
if (nrow(only_real)) {
  print(only_real %>% dplyr::select(test, CI, adjp_real, adjp_sim))
} else {
  cat("(None — try a different seed or larger n if you expect near-zero sim rejections.)\n")
}







##### update 22.10
# showing the edge weights


library(ggplot2)
simulate_linear_from_dag <- function(amat, n, seed = 1,
                                     w_mean = 0.8, w_sd = 0.2,
                                     noise_sd = 1.0,
                                     return_weights = FALSE) {
  set.seed(seed)
  stopifnot(all(rownames(amat) == colnames(amat)))
  vars <- rownames(amat)
  g <- igraph::graph_from_adjacency_matrix(amat, mode = "directed")
  stopifnot(igraph::is_dag(g))
  order <- igraph::as_ids(igraph::topo_sort(g, mode = "out"))
  
  E <- matrix(rnorm(n * length(vars), 0, noise_sd), nrow = n, ncol = length(vars))
  colnames(E) <- vars
  
  X <- matrix(0, n, length(vars)); colnames(X) <- vars
  W <- matrix(0, nrow = length(vars), ncol = length(vars),
              dimnames = list(vars, vars))  # W[parent, child]
  
  for (v in order) {
    parents <- vars[which(amat[, v] == 1L)]
    if (!length(parents)) {
      X[, v] <- E[, v]
    } else {
      w  <- rnorm(length(parents), w_mean, w_sd)     # fixed weights per parent
      mu <- as.matrix(X[, parents, drop = FALSE]) %*% w
      X[, v] <- as.numeric(mu) + E[, v]
      W[parents, v] <- w
    }
  }
  
  Xdf <- as.data.frame(X)[, rownames(amat), drop = FALSE]
  
  if (!return_weights) return(Xdf)
  
  edge_weights <- do.call(rbind, lapply(vars, function(child) {
    parents <- vars[which(amat[, child] == 1L)]
    if (!length(parents)) return(NULL)
    data.frame(from = parents, to = child, weight = W[parents, child],
               row.names = NULL, stringsAsFactors = FALSE)
  }))
  list(data = Xdf, weights = edge_weights, W = W)
}




plot_sachs_paper_circle <- function(
    dag,
    order  = c("PIP2","PLCg","Mek","Raf","JNK","p38","PKC","PKA","Akt","Erk","PIP3"),
    anchor = "PIP2",
    mirror_y = TRUE,
    title  = "Sachs graph"
){
  stopifnot(setequal(order, names(dag)))
  k <- length(order)
  theta <- seq(0, 2*pi, length.out = k + 1)[1:k]; names(theta) <- order
  theta <- theta + ((pi/2) - theta[anchor])
  x <- cos(theta); y <- sin(theta); if (mirror_y) x <- -x
  pos <- data.frame(name = order, x = as.numeric(x), y = as.numeric(y), row.names = order)
  
  el <- dagitty::edges(dag)
  g  <- igraph::graph_from_data_frame(el, directed = TRUE,
                                      vertices = data.frame(name = order))
  
  p <- ggraph::ggraph(g, layout = "manual", x = pos$x, y = pos$y) +
    ggraph::geom_edge_link(
      arrow   = grid::arrow(length = grid::unit(3, "mm"), type = "closed"),
      end_cap = ggraph::circle(3.2, "mm"),
      width   = 0.6
    ) +
    ggraph::geom_node_point(size = 12, shape = 21, fill = "black", color = "black", stroke = 0) +
    ggraph::geom_node_text(ggplot2::aes(label = name), color = "white", size = 3.8, fontface = "bold") +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0, face = "bold"))
  print(p); invisible(p)
}


library(ggraph); library(ggplot2); library(dplyr); library(igraph); library(dagitty)

amat <- get_sachs_amat()
dag  <- adj2dag(amat)
paper_order <- c("PIP2","PLCg","Mek","Raf","JNK","p38","PKC","PKA","Akt","Erk","PIP3")
dag_circ <- ensure_circular_coords(dag, order = paper_order, anchor = "PIP2", mirror_y = TRUE)

set.seed(999)
sim_true <- simulate_linear_from_dag(amat, n = 1000, seed = 999,
                                     w_mean = 0.8, w_sd = 0.2, noise_sd = 1.0,
                                     return_weights = TRUE)

plot_sachs_paper_circle(dag_circ, order = paper_order, anchor = "PIP2",
                        mirror_y = TRUE, title = "Sachs graph (mirrored)")

plot_sachs_with_edge_weights(dag_circ, edge_weights = sim_true$weights,
                             title = "Sachs DAG with TRUE simulation weights")





## now we are testing beyond lukas kook: what results do we get if we test every CI
# TEST ALL IMPLIED CIs (not just the reduced basis)


cat("\n=== TESTING ALL IMPLIED CI RELATIONS ===\n")

# Enumerate all CIs from the Sachs DAG 
# All implied CIs via brute-force d-separation (fixed) 
# include_empty: include unconditional independences (Z = ∅)
# max_Z_size  : cap |Z| to control combinatorial growth (NULL = no cap)
# minimal_only: keep only inclusion-minimal separating sets per (X,Y)
enumerate_all_cis <- function(amat,
                              include_empty = TRUE,
                              max_Z_size = NULL,
                              minimal_only = TRUE) {
  g <- adj2dag(amat)
  nodes <- rownames(amat)
  p <- length(nodes)
  
  # all subsets of vec, optionally size-capped
  subsets <- function(vec, maxk = NULL) {
    out <- list(character(0))
    if (length(vec) == 0) return(out)
    for (k in 1:length(vec)) {
      if (!is.null(maxk) && k > maxk) break
      out <- c(out, combn(vec, k, simplify = FALSE))
    }
    out
  }
  
  cis_list <- list()
  
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      X <- nodes[i]; Y <- nodes[j]
      rest <- setdiff(nodes, c(X, Y))
      allZ <- subsets(rest, maxk = max_Z_size)
      
      # gather all separating sets for this (X,Y)
      pair_sets <- list()
      for (Z in allZ) {
        if (length(Z) == 0 && !include_empty) next
        if (dagitty::dseparated(g, X = X, Y = Y, Z = Z)) {
          pair_sets[[length(pair_sets) + 1]] <- Z
        }
      }
      
      # keep only inclusion-minimal Z (optional)
      if (minimal_only && length(pair_sets) > 1) {
        keep <- rep(TRUE, length(pair_sets))
        for (a in seq_along(pair_sets)) {
          Za <- pair_sets[[a]]
          for (b in seq_along(pair_sets)) {
            if (a == b) next
            Zb <- pair_sets[[b]]
            if (length(Zb) < length(Za) && all(Zb %in% Za)) {
              keep[a] <- FALSE; break
            }
          }
        }
        pair_sets <- pair_sets[keep]
      }
      
      # append to global list
      if (length(pair_sets)) {
        for (Z in pair_sets) {
          cis_list[[length(cis_list) + 1]] <- list(X = X, Y = Y, Z = Z)
        }
      }
    }
  }
  
  cis_list
}


# Function to run CI tests on any CI list 
run_ci_tests_from_list <- function(cis_list, dat, tests = c("gcm","pcm"), alpha = 0.05) {
  res <- lapply(tests, function(tst) {
    pv <- vapply(seq_along(cis_list), function(k) {
      ci <- cis_list[[k]]
      
      # always provide two RHS parts; use "1" when Z is empty 
      rhs <- paste0(
        paste0(ci$Y, collapse = "+"),
        " | ",
        if (length(ci$Z)) paste0(ci$Z, collapse = "+") else "1"
      )
      
      fm <- reformulate(rhs, response = ci$X)
      
      # be defensive in case any formula fails
      tryCatch(
        comets(fm, dat, test = tst, coin = TRUE)$p.value,
        error = function(e) NA_real_
      )
    }, numeric(1))
    
    tibble::tibble(
      test = tst,
      CI = vapply(cis_list, function(ci) {
        if (length(ci$Z))
          paste(ci$X, "_||_", paste(ci$Y, collapse="+"), "|", paste(ci$Z, collapse=", "))
        else
          paste(ci$X, "_||_", paste(ci$Y, collapse="+"), "|", "∅")
      }, character(1)),
      p.value = pv,
      adj.p.value = p.adjust(pv, "holm"),
      rejected = p.adjust(pv, "holm") < alpha
    )
  }) |> dplyr::bind_rows()
  
  res
}


# Compare reduced vs. all implied CI sets 
cis_reduced <- cis_with_Z(amat)
res_reduced <- run_ci_tests_from_list(cis_reduced, dat_real, tests = tests, alpha = alpha)

cis_all <- enumerate_all_cis(amat, include_empty = TRUE, max_Z_size = NULL, minimal_only = TRUE)
res_all <- run_ci_tests_from_list(cis_all, dat_real, tests = tests, alpha = alpha)

summarize_ci_run <- function(tag, df) {
  df |> dplyr::summarize(
    tag = tag,
    n_tests = dplyr::n(),
    n_sig = sum(rejected),
    min_adj_p = min(adj.p.value, na.rm = TRUE),
    median_adj_p = stats::median(adj.p.value, na.rm = TRUE)
  )
}

summary_tbl <- dplyr::bind_rows(
  summarize_ci_run("Reduced (Paper-style)", res_reduced),
  summarize_ci_run("All implied", res_all)
)
print(summary_tbl)

# Save results 
readr::write_csv(res_reduced, file.path(out_dir, "real_reduced_CIs.csv"))
readr::write_csv(res_all, file.path(out_dir, "real_all_CIs_full.csv"))






################################################################################
### 
### adding an edge to the Sachs Graph 
### trying to identify which direction it could go
################################################################################




# Parse "X _||_ Y | Z1, Z2" -> c("X","Y")
parse_ci_pair <- function(ci_str) {
  parts <- strsplit(ci_str, " _\\|\\|_ ", perl = TRUE)[[1]]
  stopifnot(length(parts) >= 2)
  X <- trimws(parts[1])
  yz <- strsplit(parts[2], "\\|", perl = TRUE)[[1]]
  Y <- trimws(yz[1])
  c(X, Y)
}

# Unique unordered pairs from rejected CI strings
pairs_from_rejections <- function(ci_strings) {
  pairs <- t(vapply(ci_strings, parse_ci_pair, character(2)))
  pairs <- t(apply(pairs, 1, sort))
  pairs <- unique(pairs)
  colnames(pairs) <- c("u","v")
  as.data.frame(pairs, stringsAsFactors = FALSE)
}

# Add oriented edge if acyclic; else return NULL
add_oriented_edge_if_dag <- function(amat, u, v, direction = c("u2v","v2u")) {
  direction <- match.arg(direction)
  A <- amat
  A[u, v] <- 0; A[v, u] <- 0
  if (direction == "u2v") A[u, v] <- 1 else A[v, u] <- 1
  g <- igraph::graph_from_adjacency_matrix(A, mode = "directed")
  if (!igraph::is_dag(g)) return(NULL)
  A
}

# Run COMETs on a given adjacency
run_all_tests_for_amat <- function(amat, dat, tests = c("gcm","pcm"), alpha = 0.05) {
  res <- run_ci_tests(amat, dat, tests = tests, alpha = alpha)
  list(n_rejected = sum(res$rejected), results = res)
}

# Unexported v-structure helpers (replacement for dagitty::vstructures) 
skeleton_from_amat <- function(A) (A + t(A)) > 0

vstructures_from_amat <- function(A){
  stopifnot(all(rownames(A) == colnames(A)))
  nodes <- rownames(A); out <- list()
  for (z in nodes) {
    pars <- nodes[A[, z] == 1L]
    if (length(pars) >= 2) {
      for (i in 1:(length(pars)-1)) for (j in (i+1):length(pars)) {
        x <- pars[i]; y <- pars[j]
        if (A[x,y]==0L && A[y,x]==0L) out[[length(out)+1]] <- c(X=x,Z=z,Y=y)
      }
    }
  }
  if (!length(out)) data.frame(X=character(),Z=character(),Y=character())
  else as.data.frame(do.call(rbind,out), stringsAsFactors = FALSE)
}

markov_equivalent <- function(A1, A2){
  same_skel <- identical(skeleton_from_amat(A1), skeleton_from_amat(A2))
  norm_vs <- function(vs){
    if (!nrow(vs)) return(data.frame(X=character(),Y=character(),Z=character()))
    d <- t(apply(vs,1,function(r) c(sort(c(r["X"],r["Y"])), r["Z"])))
    d <- as.data.frame(d); names(d) <- c("X","Y","Z")
    d[order(d$X,d$Y,d$Z), , drop=FALSE]
  }
  same_vs <- identical(norm_vs(vstructures_from_amat(A1)), norm_vs(vstructures_from_amat(A2)))
  list(same_skeleton = same_skel,
       same_vstructures = same_vs,
       markov_equivalent = same_skel && same_vs)
}

# suggest & test orientations for each rejected CI pair
suggest_and_test_orientations <- function(viol_tbl, amat0, dat, tests = c("gcm","pcm"), alpha = 0.05) {
  stopifnot("CI" %in% names(viol_tbl))
  cand_pairs <- pairs_from_rejections(viol_tbl$CI)
  
  base <- run_all_tests_for_amat(amat0, dat, tests, alpha)
  base_nrej <- base$n_rejected
  
  summaries <- list(); details <- list()
  
  for (i in seq_len(nrow(cand_pairs))) {
    u <- cand_pairs$u[i]; v <- cand_pairs$v[i]
    A1 <- add_oriented_edge_if_dag(amat0, u, v, "u2v")
    A2 <- add_oriented_edge_if_dag(amat0, u, v, "v2u")
    
    s1 <- d1 <- s2 <- d2 <- NULL
    
    if (!is.null(A1)) {
      out1 <- run_all_tests_for_amat(A1, dat, tests, alpha)
      s1 <- data.frame(pair = paste(u,v,sep="—"),
                       orientation = sprintf("%s→%s", u, v),
                       valid = TRUE,
                       n_rejected = out1$n_rejected,
                       delta_vs_base = out1$n_rejected - base_nrej,
                       stringsAsFactors = FALSE)
      d1 <- transform(out1$results, pair = paste(u,v,sep="—"), orientation = sprintf("%s→%s", u, v))
    } else {
      s1 <- data.frame(pair = paste(u,v,sep="—"),
                       orientation = sprintf("%s→%s", u, v),
                       valid = FALSE,
                       n_rejected = NA_integer_,
                       delta_vs_base = NA_integer_,
                       stringsAsFactors = FALSE)
    }
    
    if (!is.null(A2)) {
      out2 <- run_all_tests_for_amat(A2, dat, tests, alpha)
      s2 <- data.frame(pair = paste(u,v,sep="—"),
                       orientation = sprintf("%s→%s", v, u),
                       valid = TRUE,
                       n_rejected = out2$n_rejected,
                       delta_vs_base = out2$n_rejected - base_nrej,
                       stringsAsFactors = FALSE)
      d2 <- transform(out2$results, pair = paste(u,v,sep="—"), orientation = sprintf("%s→%s", v, u))
    } else {
      s2 <- data.frame(pair = paste(u,v,sep="—"),
                       orientation = sprintf("%s→%s", v, u),
                       valid = FALSE,
                       n_rejected = NA_integer_,
                       delta_vs_base = NA_integer_,
                       stringsAsFactors = FALSE)
    }
    
    summaries[[i]] <- rbind(s1, s2)
    details[[i]]   <- dplyr::bind_rows(d1, d2)
  }
  
  summary_df <- do.call(rbind, summaries)
  
  winner <- summary_df %>%
    dplyr::filter(valid) %>%
    dplyr::group_by(pair) %>%
    dplyr::filter(n_rejected == min(n_rejected, na.rm = TRUE)) %>%
    dplyr::mutate(recommendation = dplyr::case_when(
      dplyr::n() > 1 ~ "Markov equivalent (direction not identifiable)",
      TRUE ~ paste("prefer", orientation)
    )) %>% dplyr::ungroup()
  
  list(
    baseline_rejections = base_nrej,
    orientation_summary = summary_df,
    recommendation      = winner,
    orientation_details = dplyr::bind_rows(details)
  )
}


# Orientation suggestions based on rejected CIs 
amat0 <- get_sachs_amat()
out <- suggest_and_test_orientations(
  viol_tbl = viol_real,
  amat0    = amat0,
  dat      = dat_real,
  tests    = tests,
  alpha    = alpha
)

out$baseline_rejections
out$orientation_summary %>% dplyr::arrange(pair, orientation)
out$recommendation       %>% dplyr::arrange(pair)
# Optionally inspect per-CI details after each orientation:
# out$orientation_details %>% dplyr::filter(rejected) %>% dplyr::arrange(pair, orientation, adj.p.value)




#################################################################################
## now we use the SACHS DAG and take away an edge.
## the data is from the simulations with the SACHS dag, so it should reject the CI
#################################################################################

## Simulate from the original DAG, then remove Akt—Erk and re-test

set.seed(4242)

# Base Sachs DAG and simulated data (same n as real data for apples-to-apples)
amat_base <- get_sachs_amat()
dat_real  <- read_sachs_observational(data_path)
sim_dat   <- simulate_linear_from_dag(
  amat_base,
  n     = nrow(dat_real),
  seed  = 4242,
  w_mean = 0.8, w_sd = 0.2, noise_sd = 1.0
)

# Remove the Akt—Erk edge (both directions to be safe)
amat_drop <- amat_base
amat_drop["Akt","Erk"] <- 0L
amat_drop["Erk","Akt"] <- 0L

# run test for CI's
res_sim_base <- run_ci_tests(amat_base, sim_dat, tests = tests, alpha = alpha)
res_sim_drop <- run_ci_tests(amat_drop, sim_dat, tests = tests, alpha = alpha)

viol_sim_base <- dplyr::filter(res_sim_base, rejected) %>% dplyr::arrange(adj.p.value)
viol_sim_drop <- dplyr::filter(res_sim_drop, rejected) %>% dplyr::arrange(adj.p.value)

cat("\n--- Simulated from ORIGINAL DAG ---\n")
cat("Total CIs tested:", nrow(res_sim_base), "\n")
cat("Rejected after Holm @ alpha =", alpha, ":", nrow(viol_sim_base), "\n")

cat("\n--- Same data, but DAG with Akt—Erk REMOVED ---\n")
cat("Total CIs tested:", nrow(res_sim_drop), "\n")
cat("Rejected after Holm @ alpha =", alpha, ":", nrow(viol_sim_drop), "\n")

if (nrow(viol_sim_drop)) {
  cat("\nTop rejections after removing Akt—Erk (should include Akt/Erk pairs):\n")
  print(
    viol_sim_drop %>%
      dplyr::select(test, CI, adj.p.value) %>%
      dplyr::arrange(adj.p.value) %>%
      dplyr::slice_head(n = 10)
  )
}

# Specifically check whether any rejected CI involves the pair (Akt, Erk)
involves_Akt_Erk <- function(ci_str) {
  # look for both names on the left side of the pipe
  # CI format: "X _||_ Y | Z1, Z2"
  # we check unordered pair in X,Y
  parts <- strsplit(ci_str, " _\\|\\|_ ", perl = TRUE)[[1]]
  if (length(parts) < 2) return(FALSE)
  X <- trimws(parts[1])
  Y <- trimws(strsplit(parts[2], "\\|", perl = TRUE)[[1]][1])
  s <- sort(c(X, Y))
  identical(s, sort(c("Akt", "Erk")))
}

viol_AktErk <- viol_sim_drop %>% dplyr::filter(vapply(CI, involves_Akt_Erk, logical(1)))

cat("\nRejected CIs that directly involve (Akt, Erk):\n")
if (nrow(viol_AktErk)) {
  print(viol_AktErk %>% dplyr::select(test, CI, adj.p.value) %>% dplyr::arrange(adj.p.value))
} else {
  cat("(none)\n")
}

### now some functions to help understand the changes
# CI diff helper: which CIs are newly implied (or removed) 
ci_strings <- function(amat){
  g <- adj2dag(amat)
  cis <- dagitty::impliedConditionalIndependencies(g)
  vapply(Filter(function(ci) length(ci$Z)>0, cis), paste, character(1))
}

ci_diff <- function(A_old, A_new){
  old <- ci_strings(A_old); new <- ci_strings(A_new)
  list(added = setdiff(new, old), removed = setdiff(old, new))
}


diff <- ci_diff(amat_base, amat_drop)

cat("\nNewly implied CIs after removing Akt—Erk:\n")
print(diff$added)

cat("\nCIs no longer implied (were implied before, not now):\n")
print(diff$removed)

# Join with your test results to see which of the NEW CIs were rejected
added_tbl <- tibble::tibble(CI = diff$added) %>%
  dplyr::left_join(res_sim_drop %>% dplyr::select(CI, test, adj.p.value, rejected),
                   by = "CI") %>%
  dplyr::arrange(dplyr::desc(rejected), adj.p.value)

cat("\nNewly implied CIs with test outcomes (after edge removal):\n")
print(added_tbl, n = nrow(added_tbl))




################################################################################
## additional edge Erk -> PIP3, then we simulate data from it, 
## we use old DAG without the edge and test 
## if we can find the (Erk, PIP3) CI statement
################################################################################


################################################################################
## Make all implied CIs visible
################################################################################

make_ci_table <- function(amat,
                          include_empty = TRUE,
                          max_Z_size    = NULL,
                          minimal_only  = TRUE) {
  cis_list <- enumerate_all_cis(
    amat,
    include_empty = include_empty,
    max_Z_size    = max_Z_size,
    minimal_only  = minimal_only
  )
  
  tibble::tibble(
    ci_id = seq_along(cis_list),
    X     = vapply(cis_list, function(ci) ci$X, character(1)),
    Y     = vapply(cis_list, function(ci) ci$Y, character(1)),
    Z_str = vapply(cis_list, function(ci) {
      if (length(ci$Z)) paste(sort(ci$Z), collapse = ", ") else "∅"
    }, character(1)),
    # keep full object in case you need it later
    ci_obj = cis_list
  )
}



# CI format from run_ci_tests_from_list:
#   "X _||_ Y1+Y2 | Z1, Z2"   or   "X _||_ Y | ∅"
involves_pair <- function(ci_str, a, b) {
  parts <- strsplit(ci_str, " _\\|\\|_ ", perl = TRUE)[[1]]
  if (length(parts) < 2) return(FALSE)
  X <- trimws(parts[1])
  Y <- trimws(strsplit(parts[2], "\\|", perl = TRUE)[[1]][1])  # Y part before '|'
  s <- sort(c(X, Y))
  identical(s, sort(c(a, b)))
}



################################################################################
## 2) Main experiment function:
##    - You manually choose X, Y, and direction ("X_to_Y" or "Y_to_X")
##    - We add that edge to the DAG
##    - Simulate data from modified DAG
##    - Test ALL original CIs on this new data
##    - Show which CIs involving (X, Y) are rejected
################################################################################

run_manual_edge_experiment <- function(amat_base,
                                       dat_template,
                                       X, Y,
                                       direction = c("X_to_Y", "Y_to_X"),
                                       tests = c("gcm", "pcm"),
                                       alpha = 0.05,
                                       sim_seed = 1234) {
  direction <- match.arg(direction)
  
  stopifnot(all(rownames(amat_base) == colnames(amat_base)))
  if (!(X %in% rownames(amat_base)) || !(Y %in% rownames(amat_base))) {
    stop("X and Y must be node names in the DAG.")
  }
  
  # Decide orientation
  if (direction == "X_to_Y") {
    from <- X; to <- Y
  } else {
    from <- Y; to <- X
  }
  
  message("Chosen pair: (", X, ", ", Y, "), direction: ", from, " → ", to)
  
  # Check we don't already have an edge
  if (amat_base[from, to] == 1L) {
    stop("Edge ", from, " → ", to, " already exists in the base DAG.")
  }
  if (amat_base[to, from] == 1L) {
    message("Note: there is already an edge ", to, " → ", from,
            " in the base DAG. You are not changing that here.")
  }
  
  # Build modified adjacency
  amat_mod <- amat_base
  amat_mod[from, to] <- 1L
  
  # Check acyclicity
  g_mod <- igraph::graph_from_adjacency_matrix(amat_mod, mode = "directed")
  if (!igraph::is_dag(g_mod)) {
    stop("Adding edge ", from, " → ", to, " creates a cycle. Choose another direction or pair.")
  }
  
  # Simulate data from modified DAG, same n as template data
  n <- nrow(dat_template)
  set.seed(sim_seed)
  sim_dat <- simulate_linear_from_dag(
    amat_mod,
    n        = n,
    seed     = sim_seed,
    w_mean   = 0.8,
    w_sd     = 0.2,
    noise_sd = 1.0
  )
  
  # Enumerate ALL implied CIs from ORIGINAL DAG
  cis_all_base <- enumerate_all_cis(
    amat_base,
    include_empty = TRUE,
    max_Z_size    = NULL,
    minimal_only  = TRUE
  )
  
  # Test those CIs on data from MODIFIED DAG
  res_mismatch <- run_ci_tests_from_list(
    cis_list = cis_all_base,
    dat      = sim_dat,
    tests    = tests,
    alpha    = alpha
  )
  
  viol_mismatch <- res_mismatch %>%
    dplyr::filter(rejected) %>%
    dplyr::arrange(adj.p.value)
  
  cat("\n--- Data from MODIFIED DAG (", from, "→", to,
      "), tested with ORIGINAL DAG CIs ---\n", sep = "")
  cat("Total CIs tested: ", nrow(res_mismatch), "\n")
  cat("Rejected after Holm @ alpha =", alpha, ":", nrow(viol_mismatch), "\n")
  
  if (nrow(viol_mismatch)) {
    cat("\nTop 10 rejections:\n")
    print(
      viol_mismatch %>%
        dplyr::select(test, CI, adj.p.value) %>%
        dplyr::slice_head(n = 10)
    )
  } else {
    cat("\nNo CI rejections.\n")
  }
  
  # Focus specifically on CIs involving (X, Y)
  viol_XY <- viol_mismatch %>%
    dplyr::filter(vapply(CI, involves_pair, logical(1), a = X, b = Y)) %>%
    dplyr::arrange(adj.p.value)
  
  cat("\nRejected CIs that involve (", X, ", ", Y, "):\n", sep = "")
  if (nrow(viol_XY)) {
    print(viol_XY %>% dplyr::select(test, CI, adj.p.value))
  } else {
    cat("(none for this pair)\n")
  }
  
  invisible(list(
    amat_mod        = amat_mod,
    sim_data        = sim_dat,
    all_results     = res_mismatch,
    all_rejections  = viol_mismatch,
    pair_rejections = viol_XY
  ))
}




# 0) Base objects as usual
amat_base <- get_sachs_amat()
dat_real  <- read_sachs_observational(data_path)

# 1) (Optional) Explore all CIs and choose one
ci_tbl <- make_ci_table(amat_base)
View(ci_tbl)  # or head(ci_tbl), or write_csv

# Example: manually pick a pair
X <- "Raf"
Y <- "PIP3"

# 2) Choose edge direction by hand:
#    "X_to_Y" means Raf → PIP3, "Y_to_X" means PIP3 → Raf
direction <- "X_to_Y"

# 3) Run the experiment
out <- run_manual_edge_experiment(
  amat_base    = amat_base,
  dat_template = dat_real,
  X            = X,
  Y            = Y,
  direction    = direction,
  tests        = tests,   # your c("gcm","pcm")
  alpha        = alpha,   # your 0.05
  sim_seed     = 2025
)

