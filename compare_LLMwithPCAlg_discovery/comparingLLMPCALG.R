write_llm_dag_task <- function(path = "llm_dag_task.txt") {
  txt <- paste0(
    "Task: Infer a DAG (Directed Acyclic Graph) from scratch for the following 11 variables measured in single-cell flow cytometry (CD3/CD28 condition; log-transformed).
    
    Output: return ONLY a CSV with columns `from,to`, each row a directed edge `from -> to`. No undirected edges.
    
    The following variables are protein names (use these names exactly):
    Raf, Mek, PLCg, PIP2, PIP3, Erk, Akt, PKA, PKC, p38, JNK
    
    Rules:
      1) Output only CSV (no prose). First line must be `from,to`.
      2) Use only the node names above; no extra nodes; no self-loops.
      3) The graph must be acyclic (a DAG).
      4) Keep it reasonably sparse; choose directions you consider most plausible.
      5) Do not include undirected edges.

  Return format example (structure only — replace with your own edges):
  from,to
  A,B
  C,D"
  )
  writeLines(txt, path)
  message("Wrote LLM DAG task to: ", normalizePath(path))
}


#LLM DAG CSV GPT5
llm_csv_text <- "from,to
PIP2,PIP3
PIP3,Akt
PIP3,PLCg
PLCg,PKC
PKC,Raf
Raf,Mek
Mek,Erk
PKC,JNK
PKC,p38
Akt,Raf
PKA,Raf
"
writeLines(llm_csv_text, "llm_dag.csv")

# ===========================================
# Kook-style falsification: PC DAG vs LLM DAG
# ===========================================

set.seed(426)


data_path <- "~/Desktop/master_thesis/code/MasterThesis/data"  # folder with cd3cd28*.xls
pc_rds    <- "pc_none_cpdag.rds"                                
llm_csv   <- "llm_dag.csv"                                      
alpha     <- 0.05                                               
tests     <- c("gcm", "pcm")                                    


# if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
# BiocManager::install(c("graph","RBGL"))
# install.packages(c("pcalg","readxl","dplyr","dagitty","comets","ranger","igraph"))

library(readxl); library(dplyr)
library(pcalg);  library(graph)
library(dagitty)
library(comets)
library(ranger)
library(igraph)

#  Canonical node order 
nms <- c("Raf","Mek","PLCg","PIP2","PIP3","Erk","Akt","PKA","PKC","p38","JNK")

# Data loader (Kook-style)
read_data_kook <- function(int = c("none","Akt","PIP2","Erk","PKC","PIP3"),
                           path = "../data") {
  int <- match.arg(int)
  file <- switch(int,
                 "none" = "cd3cd28.xls",
                 "Akt"  = "cd3cd28+aktinhib.xls",
                 "PIP2" = "cd3cd28+psitect.xls",
                 "Erk"  = "cd3cd28+u0126.xls",
                 "PKC"  = "cd3cd28+g0076.xls",
                 "PIP3" = "cd3cd28+ly.xls"
  )
  dat <- read_xls(file.path(path, file)) |>
    mutate(across(everything(), log)) |>
    select(where(is.numeric)) |>
    as.data.frame()
  stopifnot(ncol(dat) == 11)
  colnames(dat) <- nms
  dat
}

# Normalize any adjacency to a clean 11×11 integer matrix 
normalize_adj <- function(A, nodes = nms) {
  if (inherits(A, "graphNEL")) A <- as(A, "matrix")
  A <- as.matrix(A)
  if (is.null(rownames(A)) || is.null(colnames(A)))
    stop("Adjacency must have row/col names.")
  if (!setequal(rownames(A), nodes) || !setequal(colnames(A), nodes))
    stop("Adjacency names must match the expected nodes exactly.")
  A <- A[nodes, nodes, drop = FALSE]
  A[is.na(A)] <- 0
  A <- (A > 0) + 0L
  storage.mode(A) <- "integer"
  A
}

# PC CPDAG -> one representative DAG 
pc_cpdag_to_dagamat <- function(pc_fit, nodes = nms) {
  cpdag <- as(pc_fit@graph, "matrix")
  cpdag[cpdag > 0] <- 1L
  out <- pdag2dag(as(cpdag, "graphNEL"))
  if (!out$success) stop("Failed to orient CPDAG to a DAG.")
  A <- as(out$graph, "matrix")
  storage.mode(A) <- "integer"
  if (is.null(rownames(A)) || is.null(colnames(A))) {
    nds <- graph::nodes(out$graph)
    rownames(A) <- colnames(A) <- nds
  }
  normalize_adj(A, nodes)
}

#  Read LLM DAG CSV (from,to); ensure acyclicity
read_llm_dag_csv <- function(path, nodes = nms) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot(all(c("from","to") %in% names(df)))
  A <- matrix(0L, length(nodes), length(nodes), dimnames = list(nodes, nodes))
  for (k in seq_len(nrow(df))) {
    a <- trimws(df$from[k]); b <- trimws(df$to[k])
    if (!(a %in% nodes && b %in% nodes))
      stop("Unknown node in LLM CSV: ", a, " / ", b)
    if (a == b) stop("Self-loop in LLM CSV at row ", k)
    A[a, b] <- 1L
  }
  g <- igraph::graph_from_adjacency_matrix(A, mode = "directed")
  if (!igraph::is_dag(g)) stop("LLM graph is not acyclic; please return a DAG.")
  normalize_adj(A, nodes)
}

#  dagitty conversion & implied CI extraction 
make_dagitty_from_adj <- function(A, nodes = nms) {
  A <- normalize_adj(A, nodes)
  idx <- which(A == 1L, arr.ind = TRUE)
  edges <- if (nrow(idx)) {
    apply(idx, 1, function(rc) sprintf("%s -> %s", nodes[rc[1]], nodes[rc[2]]))
  } else character(0)
  spec <- paste("dag {", paste(nodes, collapse = "; "),
                if (length(edges)) ";" else "",
                paste(edges, collapse = "; "), "}")
  dagitty(spec)
}

cis_from_adj <- function(A, nodes = nms) {
  g <- make_dagitty_from_adj(A, nodes)
  cis <- impliedConditionalIndependencies(g)
  cis[sapply(cis, function(ci) length(ci$Z) > 0)]  # Kook filters to non-empty Z
}

# COMETS testing (GCM/PCM) + Holm correction 
test_ci_set <- function(dat, cis_list, test = c("gcm","pcm")) {
  test <- match.arg(test)
  if (length(cis_list) == 0L) {
    return(data.frame(CI=character(0), p.value=numeric(0),
                      adj.p.value=numeric(0), test=character(0)))
  }
  pv <- vapply(cis_list, function(ci) {
    rhs_y <- if (length(ci$Y)) paste(ci$Y, collapse = "+") else "1"
    rhs_z <- paste(ci$Z, collapse = "+")
    fm <- as.formula(paste0(ci$X, " ~ ", rhs_y, " | ", rhs_z))
    comets(fm, dat, test = test, coin = TRUE)$p.value
  }, numeric(1))
  data.frame(
    CI = sapply(cis_list, paste),
    p.value = pv,
    adj.p.value = p.adjust(pv, "holm"),
    test = test,
    row.names = NULL
  )
}

evaluate_graph_kook <- function(dat, A, label, which_tests = tests) {
  cis <- cis_from_adj(A, nodes = nms)
  outs <- lapply(which_tests, function(tt) transform(test_ci_set(dat, cis, tt), graph = label))
  do.call(rbind, outs)
}

compare_llm_vs_pc_kook <- function(dat, A_pc, A_llm, alpha = 0.05, which_tests = tests) {
  res_pc  <- evaluate_graph_kook(dat, A_pc,  "PC",  which_tests)
  res_llm <- evaluate_graph_kook(dat, A_llm, "LLM", which_tests)
  res <- rbind(res_pc, res_llm)
  res$rejected <- res$adj.p.value < alpha
  summary <- res |>
    dplyr::group_by(graph, test) |>
    dplyr::summarise(n_CIs = dplyr::n(), n_rejected = sum(rejected), .groups = "drop")
  list(results = res, summary = summary)
}

#  (Optional) apply intervention: zero incoming edges 
apply_intervention <- function(A, int) {
  if (identical(int, "none")) return(A)
  A1 <- A; A1[, int] <- 0L; A1
}

# RUN (single condition: "none") this is the first dataset of kook
X     <- read_data_kook(int = "none", path = data_path)
pcfit <- readRDS(pc_rds)
A_pc  <- pc_cpdag_to_dagamat(pcfit, nodes = nms)
A_llm <- read_llm_dag_csv(llm_csv, nodes = nms)

cmp <- compare_llm_vs_pc_kook(X, A_pc, A_llm, alpha = alpha)
print(cmp$summary)

write.csv(cmp$results, "kook_style_results.csv", row.names = FALSE)
write.csv(cmp$summary,  "kook_style_summary.csv", row.names = FALSE)

#  quick DOT exports 
save_dot <- function(A, file) {
  nds <- colnames(A); lines <- c("digraph DAG {", '  node [shape=ellipse];')
  for (i in seq_along(nds)) for (j in seq_along(nds)) if (A[i, j] == 1L)
    lines <- c(lines, sprintf('  "%s" -> "%s";', nds[i], nds[j]))
  writeLines(c(lines, "}"), file)
}

