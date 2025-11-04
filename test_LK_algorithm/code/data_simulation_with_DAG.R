# data_simulation_with_DAG.R
# Utilities to (1) build a DAG, (2) perturb it, and (3) simulate + save datasets.


requireNamespace("igraph", quietly = TRUE)
requireNamespace("dagitty", quietly = TRUE)
requireNamespace("readr", quietly = TRUE)



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
#
#plot graph
#
# Checks acyclicity using igraph::is_dag(...).
is_dag_mat <- function(A) {
  igraph::is_dag(igraph::graph_from_adjacency_matrix(A, mode = "directed"))
}

# makes a graphical dag representation needed for the CI search
adj2dag <- function(adj) {
  nodes <- rownames(adj); s <- "dag {"
  for (i in seq_len(nrow(adj))) for (j in seq_len(ncol(adj)))
    if (adj[i,j]==1) s <- paste(s, nodes[i], "->", nodes[j], ";")
  dagitty::dagitty(paste0(s, "}"))
}

# --- 2c) Perturb the DAG (add/remove/flip edges, keep acyclic) ---
perturb_dag <- function(amat, add = 0, remove = 0, flip = 0, seed = 1L) {
  set.seed(seed)
  stopifnot(is_dag_mat(amat))
  A <- amat
  
  # try to add one edge at a time (randomly) if graph stays a DAG
  try_add <- function(A) {
    Z <- which(A == 0, arr.ind = TRUE)
    Z <- Z[Z[,1] != Z[,2], , drop = FALSE]
    if (!nrow(Z)) return(A)
    for (k in sample(seq_len(nrow(Z)))) {
      i <- Z[k,1]; j <- Z[k,2]
      B <- A; B[i,j] <- 1
      if (is_dag_mat(B)) return(B)
    }
    A
  }
  # remove one random existing edge
  try_remove <- function(A) {
    O <- which(A == 1, arr.ind = TRUE)
    if (!nrow(O)) return(A)
    k <- sample(seq_len(nrow(O)), 1)
    A[O[k,1], O[k,2]] <- 0
    A
  }
  # flip u->v to v->u if still DAG (avoid creating 2-cycles)
  try_flip <- function(A) {
    O <- which(A == 1, arr.ind = TRUE)
    if (!nrow(O)) return(A)
    for (k in sample(seq_len(nrow(O)))) {
      i <- O[k,1]; j <- O[k,2]
      if (A[j,i] == 1) next
      B <- A; B[i,j] <- 0; B[j,i] <- 1
      if (is_dag_mat(B)) return(B)
    }
    A
  }
  
  for (k in seq_len(remove)) A <- try_remove(A)
  for (k in seq_len(flip))   A <- try_flip(A)
  for (k in seq_len(add))    A <- try_add(A)
  
  stopifnot(is_dag_mat(A))
  A
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

# --- 4) Save utilities (Sachs-like naming, Kook-friendly) ---
.sachs_base <- function(int = c("none","Akt","PIP2","Erk","PKC","PIP3")) {
  int <- match.arg(int)
  switch(int,
         "none" = "cd3cd28",
         "Akt"  = "cd3cd28+aktinhib",
         "PIP2" = "cd3cd28+psitect",
         "Erk"  = "cd3cd28+u0126",
         "PKC"  = "cd3cd28+g0076",
         "PIP3" = "cd3cd28+ly"
  )
}

write_dataset <- function(dat, path = "../data",
                          base = .sachs_base("none"),
                          suffix = NULL, write_xlsx = FALSE) {  # <- keep FALSE
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  name <- if (!is.null(suffix) && nzchar(suffix)) paste0(base, "_", suffix) else base
  readr::write_csv(dat, file.path(path, paste0(name, ".csv")))
  # no XLSX output
  invisible(file.path(path, paste0(name, ".csv")))
}


# --- 5) High-level generators ---

# (A) Generate one dataset from a (possibly perturbed) DAG
generate_dataset <- function(n = 1000, seed = 42,
                             add = 0, remove = 0, flip = 0,
                             path = "../data",
                             base = .sachs_base("none"),
                             suffix = NULL,
                             save_graph_files = TRUE) {
  amat0 <- get_dag()
  # plot dag
  amat  <- perturb_dag(amat0, add = add, remove = remove, flip = flip, seed = seed)
  dat   <- simulate_from_dag(amat, n = n, seed = seed)
  
  # save data
  write_dataset(dat, path = path, base = base, suffix = suffix)
  
  # optionally save graph artifacts for traceability
  if (save_graph_files) {
    # edge list
    el <- which(amat == 1, arr.ind = TRUE)
    if (nrow(el)) {
      edges <- data.frame(from = rownames(amat)[el[,1]], to = colnames(amat)[el[,2]])
    } else edges <- data.frame(from = character(0), to = character(0))
    readr::write_csv(edges, file.path(path, paste0(base, if (!is.null(suffix) && nzchar(suffix)) paste0("_", suffix), "_edges.csv")))
    
    # dagitty
    dag <- adj2dag(amat)
    writeLines(as.character(dag),
               con = file.path(path, paste0(base, if (!is.null(suffix) && nzchar(suffix)) paste0("_", suffix), "_dagitty.txt")))
  }
  
  invisible(list(amat_true = amat0, amat_used = amat, data = dat))
}

#plot(dag)

# (B) Generate many datasets for power curves: n_grid × reps
generate_grid <- function(n_grid = c(200, 500, 1000),
                          reps = 5,
                          seed = 1L,
                          add = 0, remove = 0, flip = 0,
                          path = "../data",
                          base = .sachs_base("none")) {
  set.seed(seed)
  out <- list()
  for (n in n_grid) {
    for (r in seq_len(reps)) {
      suf <- paste0("n", n, "_rep", r,
                    if (add|remove|flip) paste0("_add",add,"_rem",remove,"_flip",flip) else "")
      res <- generate_dataset(n = n, seed = sample.int(1e9, 1),
                              add = add, remove = remove, flip = flip,
                              path = path, base = base, suffix = suf,
                              save_graph_files = (r == 1))  # save graph once per n
      out[[length(out)+1]] <- list(n = n, rep = r, files_base = paste0(base, "_", suf))
    }
  }
  invisible(out)
}

# --- 6) Backwards-compatible one-shot generator (observational only) ---
generate_observational <- function(n = 1000, seed = 42, path = "../data") {
  generate_dataset(n = n, seed = seed, add = 0, remove = 0, flip = 0,
                   path = path, base = .sachs_base("none"), suffix = NULL)
}


# 1) Baseline: data from the original DAG (observational)
generate_observational(n = 1000, seed = 42, path = "../data")   # -> ../data/cd3cd28.csv

# 2) Perturbed: same n, but flip 2 edges in the DAG before simulating
generate_dataset(
  n = 1000,
  seed = 43,
  flip = 2,                      # change a few edges (keeps graph acyclic)
  suffix = "n1000_flip2",        # -> ../data/cd3cd28_n1000_flip2.csv
  path = "../data",
  save_graph_files = TRUE        # also writes _edges.csv and _dagitty.txt
)

generate_dataset(
  n = 1000,
  seed = 43,
  remove = 2,                 # remove 2 random edges (keeps DAG acyclic)
  suffix = "n1000_remove2",   # -> ../data/cd3cd28_n1000_remove2.csv
  path = "../data",
  save_graph_files = TRUE     # also writes _edges.csv and _dagitty.txt
)


generate_dataset(
  n = 200,
  seed = 43,
  flip = 2,                      # change a few edges (keeps graph acyclic)
  suffix = "n1000_flip2",        # -> ../data/cd3cd28_n1000_flip2.csv
  path = "../data",
  save_graph_files = TRUE        # also writes _edges.csv and _dagitty.txt
)







# short and controlled version

## ======================
#Perturb & Detect (linear model, traceable)
## ======================


plot_dag_adj <- function(amat, file = NULL) {
  g <- adj2dag(amat)
  if (is.null(file)) {
    plot(g)
  } else {
    png(file, width = 1000, height = 800)
    plot(g)
    dev.off()
  }
}

edges_df <- function(amat) {
  el <- which(amat == 1, arr.ind = TRUE)
  if (!nrow(el)) return(data.frame(from=character(), to=character()))
  data.frame(from = rownames(amat)[el[,1]],
             to   = colnames(amat)[el[,2]])
}

diff_edges <- function(before, after) {
  # before/after: data.frames with columns from,to
  suppressWarnings({
    added   <- dplyr::anti_join(after,  before, by = c("from","to"))
    removed <- dplyr::anti_join(before, after,  by = c("from","to"))
  })
  # flips = where reverse appears
  flips <- merge(added, transform(removed, from2 = to, to2 = from),
                 by.x=c("from","to"), by.y=c("from2","to2"))
  if (nrow(flips)) {
    flipped <- data.frame(from = flips$to, to = flips$from)
    added   <- dplyr::anti_join(added,   flipped, by=c("from","to"))
    removed <- dplyr::anti_join(removed, flipped, by=c("from","to"))
  } else flipped <- added[0,]
  list(added = added, removed = removed, flipped = flipped)
}

cis_with_Z <- function(amat) {
  g <- adj2dag(amat)
  cis <- impliedConditionalIndependencies(g)
  cis[vapply(cis, function(x) length(x$Z) > 0, logical(1))]
}
ci_keys <- function(cis) sapply(cis, paste)
compare_ci_sets <- function(amat_true, amat_alt) {
  c_true <- cis_with_Z(amat_true); k_true <- ci_keys(c_true)
  c_alt  <- cis_with_Z(amat_alt);  k_alt  <- ci_keys(c_alt)
  list(only_in_true = setdiff(k_true, k_alt),
       only_in_alt  = setdiff(k_alt,  k_true),
       in_both      = intersect(k_true, k_alt))
}

# simple linear SEM consistent with DAG
simulate_linear_from_dag <- function(amat, n = 1000, seed = 1L,
                                     w_mean = 0.8, w_sd = 0.2, noise_sd = 1) {
  set.seed(seed)
  stopifnot(all(rownames(amat) == colnames(amat)))
  vars <- rownames(amat)
  g <- igraph::graph_from_adjacency_matrix(amat, mode="directed")
  stopifnot(igraph::is_dag(g))
  order <- igraph::as_ids(igraph::topo_sort(g, mode="out"))
  X <- matrix(0, n, length(vars)); colnames(X) <- vars
  for (v in order) {
    parents <- vars[which(amat[, v]==1L)]
    if (!length(parents)) {
      X[, v] <- rnorm(n)  # roots
    } else {
      w <- rnorm(length(parents), w_mean, w_sd)  # one weight per parent (fixed across rows)
      mu <- as.matrix(X[, parents, drop=FALSE]) %*% w
      X[, v] <- as.numeric(mu) + rnorm(n, 0, noise_sd)
    }
  }
  as.data.frame(X)[, nms, drop=FALSE]
}

# --- experiment parameters ---
out_path <- "../data"
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

n_samples <- 1000
seed_sim  <- 101
flip_k    <- 2                     # number of edge flips to introduce
basename <- .sachs_base("none")    # "cd3cd28"
suffix   <- paste0("linear_n", n_samples, "_flip", flip_k)

# --- 1) true DAG (original), plot & save dagitty text ---
amat_true <- get_dag()
plot_dag_adj(amat_true, file.path(out_path, paste0(basename, "_TRUE.png")))
writeLines(as.character(adj2dag(amat_true)),
           file.path(out_path, paste0(basename, "_TRUE_dagitty.txt")))

# --- 2) perturb DAG (flip a few edges), plot & save ---
amat_alt <- perturb_dag(amat_true, flip = flip_k, seed = seed_sim)
plot_dag_adj(amat_alt,  file.path(out_path, paste0(basename, "_PERTURBED.png")))
writeLines(as.character(adj2dag(amat_alt)),
           file.path(out_path, paste0(basename, "_PERTURBED_dagitty.txt")))

# --- 3) log which edges changed ---
e_true <- edges_df(amat_true)
e_alt  <- edges_df(amat_alt)
chg    <- diff_edges(e_true, e_alt)
readr::write_csv(chg$added,   file.path(out_path, paste0(basename, "_", suffix, "_ADDED.csv")))
readr::write_csv(chg$removed, file.path(out_path, paste0(basename, "_", suffix, "_REMOVED.csv")))
readr::write_csv(chg$flipped, file.path(out_path, paste0(basename, "_", suffix, "_FLIPPED.csv")))

# (optional) 4) compare CI sets: what statements changed?
ci_diff <- compare_ci_sets(amat_true, amat_alt)
writeLines(ci_diff$only_in_true, file.path(out_path, paste0(basename, "_", suffix, "_CI_only_in_TRUE.txt")))
writeLines(ci_diff$only_in_alt,  file.path(out_path, paste0(basename, "_", suffix, "_CI_only_in_PERTURBED.txt")))

# --- 5) simulate LINEAR data from the PERTURBED DAG and save ---
dat_lin <- simulate_linear_from_dag(amat_alt, n = n_samples, seed = seed_sim)
csv_path <- write_dataset(dat_lin, path = out_path, base = basename, suffix = suffix)
message("Wrote: ", csv_path)






## Erkenntnis: Anderer DAG kann auch zu gleichen CI's führen! daher nicht immer widerlegt, wenn dag anders wird









