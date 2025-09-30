# agent_tool_based_dag_discovery

# Agentic DAG Tools for LK/COMETS Falsification (single dataset)
# -----------------------------------------------------------------------------
# This file defines a clean toolbelt of pure R functions and a small controller
# to run a propose -> test -> (optional) refine loop while keeping the testing
# identical to Lukas Kook's paper (COMETS GCM/PCM with Holm correction).
# -----------------------------------------------------------------------------

# Dependencies ----------------------------------------------------------------
needed <- c("jsonlite","dagitty","bnlearn","pcalg","httr2","tidyverse","readxl","comets")
miss <- setdiff(needed, rownames(installed.packages()))
if (length(miss)) install.packages(miss)

suppressPackageStartupMessages({
  library(jsonlite)
  library(dagitty)
  library(bnlearn)
  library(pcalg)
  library(httr2)
  library(tidyverse)
  library(readxl)
  library(comets)
})

# Auto-load .env if available (for OPENAI_API_KEY and model names)
if (requireNamespace("dotenv", quietly = TRUE)) {
  try(dotenv::load_dot_env(".env"), silent = TRUE)
}
library(dotenv)
dotenv::load_dot_env("~/Desktop/master_thesis/code/MasterThesis/one_shot_llm_dag_discovery/.env")

# Config ----------------------------------------------------------------------
alpha_default <- 0.05
TESTS_DEFAULT <- c("gcm","pcm")
MODEL <- Sys.getenv("OPENAI_MODEL", unset = "gpt-4o-mini")        # for propose
EDIT_MODEL <- Sys.getenv("OPENAI_EDIT_MODEL", unset = MODEL)        # for edits

# Utility: stop if no API key present
stop_if_no_key <- function() {
  if (!nzchar(Sys.getenv("OPENAI_API_KEY"))) {
    stop("Please set OPENAI_API_KEY (e.g., via .env or Sys.setenv).", call. = FALSE)
  }
}

# Sachs variable names (order must match data columns)
nms <- c("Raf","Mek","PLCg","PIP2","PIP3","Erk","Akt","PKA","PKC","p38","JNK")

# Baseline DAGs + IO -----------------------------------------------------------
get_dag <- function(which = c("consensus","sachs"), int = c("none","Akt","PIP2","Erk","PKC","PIP3")) {
  which <- match.arg(which); int <- match.arg(int)
  consensus <- matrix(c(
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
  ), ncol=11, byrow=TRUE)
  sachs <- matrix(c(
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
  ), ncol=11, byrow=TRUE)
  amat <- if (which=="consensus") consensus else sachs
  colnames(amat) <- nms; rownames(amat) <- nms
  if (int != "none") amat[, int] <- 0
  amat
}

read_data <- function(int = c("none","Akt","PIP2","Erk","PKC","PIP3"), path = "../data") {
  int <- match.arg(int)
  file <- switch(int,
                 none = "cd3cd28.xls", Akt = "cd3cd28+aktinhib.xls", PIP2 = "cd3cd28+psitect.xls",
                 Erk = "cd3cd28+u0126.xls", PKC = "cd3cd28+g0076.xls", PIP3 = "cd3cd28+ly.xls")
  dat <- readxl::read_xls(file.path(path, file)) |> mutate(across(everything(), log))
  colnames(dat) <- nms
  dat
}

adj2dag <- function(adj_matrix) {
  nodes <- rownames(adj_matrix)
  dag_string <- "dag {"
  for (i in seq_len(nrow(adj_matrix))) for (j in seq_len(ncol(adj_matrix))) {
    if (adj_matrix[i,j] == 1) dag_string <- paste(dag_string, nodes[i], "->", nodes[j], ";")
  }
  dagitty(paste(dag_string, "}"))
}

# ---------------------------- Tools (pure functions) -------------------------

# 1) Propose adjacency matrix via LLM -----------------------------------------
# Input: variables (character vector)
# Output: list(variables, adjacency (p x p numeric matrix), raw json)
tool_propose_matrix <- function(variables, model = MODEL, temperature = 0) {
  stop_if_no_key()
  system_msg <- paste(
    "Return ONLY JSON with keys 'variables' and 'adjacency'.",
    "'variables' must exactly match the provided order.",
    "'adjacency' is a square 0/1 matrix where (i,j)=1 means variables[i] -> variables[j].",
    "Diagonal must be 0. Graph must be acyclic (DAG). No extra text.")
  user_msg <- jsonlite::toJSON(list(
    variables = unname(variables),
    constraints = list(diagonal_zero = TRUE, binary_entries = TRUE, acyclic = TRUE,
                       max_parents_per_node = 3, total_edges_range = c(10,25))
  ), auto_unbox = TRUE)
  body <- list(model = model, temperature = temperature,
               messages = list(list(role="system", content=system_msg),
                               list(role="user", content=as.character(user_msg))))
  req <- request("https://api.openai.com/v1/chat/completions") |>
    req_headers(Authorization = paste("Bearer", Sys.getenv("OPENAI_API_KEY")),
                "Content-Type" = "application/json") |>
    req_body_json(body)
  resp <- req_perform(req)
  resp_json <- resp_body_json(resp, simplifyVector = FALSE)
  if (!is.null(resp_json$error)) stop(paste0("OpenAI API error: ", resp_json$error$message))
  choice1 <- resp_json$choices[[1]]
  content <- choice1$message$content
  if (is.list(content)) {
    content <- paste(vapply(content, function(part) {
      if (!is.null(part$text)) part$text else if (!is.null(part$content)) part$content else ""
    }, character(1)), collapse = "")
  }
  txt <- as.character(content)
  out <- try(jsonlite::fromJSON(txt), silent = TRUE)
  if (inherits(out, "try-error")) stop("Agent did not return valid JSON. Raw output:\n", txt)
  if (!all(c("variables","adjacency") %in% names(out))) stop("Agent JSON missing required keys.")
  if (!identical(unname(out$variables), unname(variables))) stop("Agent variables do not match input order.")
  A <- as.matrix(out$adjacency); mode(A) <- "numeric"
  colnames(A) <- out$variables; rownames(A) <- out$variables
  list(variables = out$variables, adjacency = A, raw = txt)
}

# 2) Validate adjacency matrix -------------------------------------------------
# Returns list(ok=TRUE/FALSE, errors=character(), dag=dagitty or NULL)
tool_validate_matrix <- function(A, variables) {
  errs <- c()
  ok <- TRUE; dag <- NULL
  p <- length(variables)
  if (!is.matrix(A)) { errs <- c(errs, "Adjacency is not a matrix"); ok <- FALSE }
  else if (!all(dim(A) == c(p,p))) { errs <- c(errs, "Wrong shape (not p x p)"); ok <- FALSE }
  if (ok && (!all(rownames(A) == variables) || !all(colnames(A) == variables))) {
    errs <- c(errs, "Row/column names must match variables order"); ok <- FALSE
  }
  if (ok && any(is.na(A))) { errs <- c(errs, "Adjacency contains NA"); ok <- FALSE }
  if (ok && !all(A %in% c(0,1))) { errs <- c(errs, "Adjacency must be binary 0/1"); ok <- FALSE }
  if (ok && !all(diag(A) == 0)) { errs <- c(errs, "Diagonal must be all zeros"); ok <- FALSE }
  if (ok) {
    dag_try <- try(adj2dag(A), silent = TRUE)
    if (inherits(dag_try, "try-error")) { errs <- c(errs, "Failed to construct dagitty DAG"); ok <- FALSE }
    else if (!isAcyclic(dag_try)) { errs <- c(errs, "Graph has cycles (not a DAG)"); ok <- FALSE }
    else dag <- dag_try
  }
  list(ok = ok, errors = errs, dag = dag)
}

# 3) Extract testable CIs (|Z| > 0) -------------------------------------------
# Returns a list of CI objects as given by dagitty (each has $X,$Y,$Z)
tool_extract_cis <- function(dag) {
  cis <- impliedConditionalIndependencies(dag)
  cis[unlist(lapply(cis, function(x) length(x$Z) > 0))]
}

# 4) LK/COMETS tests on one dataset -------------------------------------------
# Returns tidy tibble with columns: test, n_cis_tested, falsified, CI, p.value, adj.p.value
tool_lk_test <- function(data, cis, test = c("gcm","pcm"), alpha = alpha_default) {
  test <- match.arg(test)
  if (length(cis) == 0) {
    return(tibble(test = test, n_cis_tested = 0L, falsified = FALSE,
                  CI = character(), p.value = numeric(), adj.p.value = numeric()))
  }
  pv <- vapply(cis, function(ci) {
    fm <- reformulate(paste0(paste0(ci$Y, collapse = "+"), "|", paste0(ci$Z, collapse = "+")),
                      response = ci$X)
    suppressWarnings(comets(fm, data, test = test, coin = TRUE)$p.value)
  }, numeric(1))
  adjp <- p.adjust(pv, method = "holm")
  tibble(test = test, n_cis_tested = length(pv), falsified = any(adjp < alpha),
         CI = sapply(cis, paste), p.value = as.numeric(pv), adj.p.value = as.numeric(adjp))
}

# 5) Decide across tests -------------------------------------------------------
# Input: list of tibble rows for each test; Output: compact decision list
tool_decide <- function(results_by_test) {
  tib <- bind_rows(results_by_test)
  list(
    falsified_gcm = any(tib$test=="gcm" & tib$adj.p.value < alpha_default),
    falsified_pcm = any(tib$test=="pcm" & tib$adj.p.value < alpha_default),
    n_cis_tested = tib %>% group_by(test) %>% summarise(n=unique(n_cis_tested), .groups='drop'),
    n_violations_by_test = tib %>% group_by(test) %>% summarise(n_viol = sum(adj.p.value < alpha_default), .groups='drop')
  )
}

if (FALSE) {
# 6) Suggest edits via LLM -----------------------------------------------------
# Input: adjacency, variables, violations (tibble with CI and adj.p.value), budget_k
# Output: data.frame(op, from, to) with op in {add, remove, reverse}
tool_suggest_edits <- function(adjacency, variables, violations, budget_k = 3,
                               model = EDIT_MODEL, temperature = 0) {
  stop_if_no_key()
  # Prepare a compact payload: current edges + top-K violated CIs
  idx <- which(adjacency == 1, arr.ind = TRUE)
  edges <- as.data.frame(cbind(from = variables[idx[,1]], to = variables[idx[,2]]), stringsAsFactors = FALSE)
  v_top <- violations %>% arrange(adj.p.value) %>% slice_head(n = budget_k)
  sys <- paste(
    "You propose minimal edits to a DAG to address rejected CIs.",
    "Return ONLY JSON array of edits: [{\"op\":\"add|remove|reverse\", \"from\":\"X\", \"to\":\"Y\"}, ...].",
    "Keep graph acyclic; at most", budget_k, "edits. No extra text.")
  usr <- toJSON(list(
    variables = variables,
    edges = edges,
    violated_cis = v_top %>% transmute(CI = CI, adj_p = adj.p.value, test = test),
    rules = list(
      op = c("add","remove","reverse"),
      semantics = "(i,j)=1 means edge variables[i] -> variables[j]",
      keep_acyclic = TRUE, max_edits = budget_k)
  ), auto_unbox = TRUE)
  body <- list(model = model, temperature = temperature,
               messages = list(list(role="system", content=sys), list(role="user", content=as.character(usr))))
  req <- request("https://api.openai.com/v1/chat/completions") |>
    req_headers(Authorization = paste("Bearer", Sys.getenv("OPENAI_API_KEY")),
                "Content-Type" = "application/json") |>
    req_body_json(body)
  resp <- req_perform(req)
  resp_json <- resp_body_json(resp, simplifyVector = FALSE)
  if (!is.null(resp_json$error)) stop(paste0("OpenAI API error: ", resp_json$error$message))
  content <- resp_json$choices[[1]]$message$content
  if (is.list(content)) {
    content <- paste(vapply(content, function(part) if (!is.null(part$text)) part$text else "", character(1)), collapse = "")
  }
  edits <- try(fromJSON(as.character(content)), silent = TRUE)
  if (inherits(edits, "try-error")) stop("LLM did not return valid JSON edits. Raw output:\n", as.character(content))
  edits <- as_tibble(edits) %>% mutate(op = tolower(op)) %>%
    filter(op %in% c("add","remove","reverse")) %>% distinct()
  # ensure nodes exist
  if (nrow(edits)) {
    stopifnot(all(edits$from %in% variables), all(edits$to %in% variables))
  }
  edits
}
}

# 7) Apply edits safely --------------------------------------------------------
# Returns updated adjacency if valid and acyclic; errors if invalid
tool_apply_edits <- function(adjacency, variables, edits) {
  A <- adjacency
  name_to_idx <- setNames(seq_along(variables), variables)
  apply_one <- function(A, op, from, to) {
    i <- name_to_idx[[from]]; j <- name_to_idx[[to]]
    if (op == "add") {
      A[i,j] <- 1
    } else if (op == "remove") {
      A[i,j] <- 0
    } else if (op == "reverse") {
      A[i,j] <- 0; A[j,i] <- 1
    }
    A
  }
  for (k in seq_len(nrow(edits))) {
    A_new <- apply_one(A, edits$op[k], edits$from[k], edits$to[k])
    # validate acyclicity after each edit
    v <- tool_validate_matrix(A_new, variables)
    if (!v$ok) stop("Edit produced invalid graph: ", paste(v$errors, collapse = "; "))
    A <- A_new
  }
  A
}



# ---------------------------- Controller (agent loop) -------------------------
# Runs propose -> test -> (optional) refine for up to max_iter iterations.
# Returns a list with adjacency_final, history (results per iter), and summaries.
agent_loop <- function(data,
                       variables = colnames(data),
                       max_iter = 0,
                       edit_budget = 3,
                       tests = TESTS_DEFAULT,
                       alpha = alpha_default) {
  stopifnot(identical(unname(variables), unname(colnames(data))))
  # Iteration 0: propose
  message("Proposing initial DAG via LLM...")
  prop <- tool_propose_matrix(variables, model = MODEL, temperature = 0)
  A_cur <- prop$adjacency
  v <- tool_validate_matrix(A_cur, variables)
  if (!v$ok) stop("Initial proposal invalid: ", paste(v$errors, collapse = "; "))
  
  history <- list()
  iter <- 0
  repeat {
    message(sprintf("Testing (iter %d)...", iter))
    dag <- v$dag
    cis <- tool_extract_cis(dag)
    res_list <- lapply(tests, function(tst) tool_lk_test(data, cis, test = tst, alpha = alpha))
    res <- bind_rows(res_list)
    res$iter <- iter
    history[[length(history)+1]] <- res
    
    fals_gcm <- any(res$test=="gcm" & res$adj.p.value < alpha)
    fals_pcm <- any(res$test=="pcm" & res$adj.p.value < alpha)
    if ((!fals_gcm && !fals_pcm) || iter >= max_iter) break
    
    # Prepare violations table for the LLM editor
    viol <- res %>% filter(adj.p.value < alpha) %>% arrange(adj.p.value)
    if (nrow(viol) == 0) break  # should not happen given condition, but safe
    
    message("Suggesting edits via LLM...")
    edits <- tool_suggest_edits(A_cur, variables, viol, budget_k = edit_budget, model = EDIT_MODEL)
    if (!nrow(edits)) break
    message("Applying edits...")
    A_new <- tool_apply_edits(A_cur, variables, edits)
    A_cur <- A_new
    v <- tool_validate_matrix(A_cur, variables)
    if (!v$ok) break
    iter <- iter + 1
  }
  
  results_all <- bind_rows(history)
  summary_tbl <- results_all %>%
    group_by(iter, test) %>%
    summarise(n_cis_tested = unique(n_cis_tested), falsified = any(adj.p.value < alpha),
              n_violations = sum(adj.p.value < alpha), .groups = "drop")
  list(adjacency_final = A_cur, history = results_all, summary = summary_tbl)
}

# ---------------------------- Example usage (single dataset) -----------------
# NOTE: Uncomment to run in your project after setting up .env and data path.
#
DATA_PATH <- "../data"  # adjust as needed
D_obs <- read_data(int = "none", path = DATA_PATH)
stopifnot(identical(colnames(D_obs), nms))
#
# # One-shot (max_iter = 0)
# out0 <- agent_loop(D_obs, variables = nms, max_iter = 0, edit_budget = 3)
# print(out0$summary)
#
# # Agentic (up to 2 refinement rounds)
out2 <- agent_loop(D_obs, variables = nms, max_iter = 2, edit_budget = 3)
print(out2$summary)
