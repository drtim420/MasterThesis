#data_simulation_with_DAG

# --- 1) Nodes (same order as your code) ---
nms <- c("Raf","Mek","PLCg","PIP2","PIP3","Erk","Akt","PKA","PKC","p38","JNK")

# --- 2) NEW DAG in matrix form (row -> col is an arrow) ---
get_dag <- function() {
  idx <- function(x) match(x, nms)
  M <- matrix(0, length(nms), length(nms), dimnames = list(nms, nms))
  add <- function(m, a, b) { m[idx(a), idx(b)] <- 1; m }
  
  M <- add(M,"PIP2","PLCg")
  M <- add(M,"PLCg","PIP3")
  M <- add(M,"PIP3","PKA"); M <- add(M,"PIP3","Akt")
  M <- add(M,"PLCg","PKC")
  M <- add(M,"PKC","Raf"); M <- add(M,"PKC","JNK"); M <- add(M,"PKC","p38")
  M <- add(M,"Raf","Mek")
  M <- add(M,"Mek","Erk")
  M <- add(M,"PKA","Erk"); M <- add(M,"PKA","Akt"); M <- add(M,"PKA","JNK")
  M <- add(M,"JNK","p38")
  M <- add(M,"Akt","Mek"); M <- add(M,"Akt","PKC")
  M
}

# (optional) dagitty helper
adj2dag <- function(adj) {
  nodes <- rownames(adj); s <- "dag {"
  for (i in seq_len(nrow(adj))) for (j in seq_len(ncol(adj)))
    if (adj[i,j]==1) s <- paste(s, nodes[i], "->", nodes[j], ";")
  dagitty::dagitty(paste0(s, "}"))
}

# --- 3) Simulator for observational data (nonlinear SEM on log scale) ---
simulate_from_dag <- function(amat, n = 1000, seed = 1L) {
  set.seed(seed)
  stopifnot(all(rownames(amat)==colnames(amat)))
  vars <- rownames(amat)
  
  g <- igraph::graph_from_adjacency_matrix(amat, mode="directed")
  if (!igraph::is_dag(g)) stop("DAG has a cycle.")
  order <- igraph::topo_sort(g, mode="out") |> igraph::as_ids()
  
  softplus <- function(x) log1p(exp(x))
  X <- matrix(0, n, length(vars)); colnames(X) <- vars
  
  for (v in order) {
    parents <- vars[which(amat[, v]==1L)]
    if (!length(parents)) {
      z <- rnorm(n,0,1)
    } else {
      w <- rnorm(length(parents), 0.9, 0.2)
      nonlin <- Reduce(`+`, Map(function(p, wi) wi * tanh(X[,p]), parents, w))
      inter <- 0
      if (length(parents) >= 2)
        for (i in seq_along(parents)) for (j in seq_len(i-1))
          inter <- inter + 0.15 * X[, parents[i]] * X[, parents[j]]
      z <- nonlin + inter + rnorm(n,0,0.5) + rnorm(1,0,0.2)
    }
    X[, v] <- log( softplus(z + rnorm(n,0,0.3)) )
  }
  as.data.frame(X)[, nms, drop=FALSE]
}

# --- 4) Generate & save like the paper (observational only) ---
generate_observational <- function(n = 1000, seed = 42, path = "../data") {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  amat <- get_dag()
  dat  <- simulate_from_dag(amat, n = n, seed = seed)
  # Save as CSV with Sachs-like basename
  readr::write_csv(dat, file.path(path, "cd3cd28.csv"))
  # (optional) also XLSX if desired:
  if (requireNamespace("writexl", quietly = TRUE))
    writexl::write_xlsx(dat, file.path(path, "cd3cd28.xlsx"))
  invisible(list(amat = amat, data = dat))
}

# --- run once ---
# install.packages(c("igraph","dagitty","readr"))
out <- generate_observational(n = 1000, seed = 42)
