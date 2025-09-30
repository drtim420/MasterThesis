# One-shot LLM DAG -> LK (Kook) CI testing pipeline (single dataset, no blinding)
# ---------------------------------------------------------------------------
# This script:
#  1) Reads ONE observational dataset (Sachs proteins).
#  2) Calls an LLM to propose a DAG as a BINARY ADJACENCY MATRIX given ONLY the variable names.
#  3) Validates the matrix (shape, binary, diag=0, acyclic).
#  4) Tests the agent DAG using the EXACT LK/Kook workflow (COMETS GCM/PCM + Holm),
#     identical to the paper's approach for CI generation and testing.
#  5) Runs the same test for the Consensus DAG as a baseline.
#
# Requirements:
#  - R >= 4.2
#  - Packages: jsonlite, dagitty, bnlearn, pcalg, httr2, comets, tidyverse, readxl
#  - OpenAI API key in env: Sys.setenv(OPENAI_API_KEY = "...")
#  - The Sachs .xls files available (path configured below)
# ---------------------------------------------------------------------------

library(dotenv)

# if you just updated .env and want to refresh in the SAME session:
Sys.unsetenv("OPENAI_API_KEY")      # clear the old value
dotenv::load_dot_env(".env")        # reload from file

# quick sanity check (don’t print full key)
nchar(Sys.getenv("OPENAI_API_KEY")) > 0




# (optional sanity check)
substr(Sys.getenv("OPENAI_API_KEY"), 1, 6)  # don't print full key



# =========================
# 0) Libraries & config
# =========================
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

set.seed(1)
alpha <- 0.05
TESTS <- c("gcm","pcm")
DATA_PATH <- "~/Desktop/master_thesis/code/MasterThesis/one_shot_llm_dag_discovery/data"   # adjust if needed
SAVE_RESULTS <- TRUE
RESULTS_DIR <- "../results"

# =========================
# 1) Sachs names + DAG helpers (from paper, kept identical in spirit)
# =========================
nms <- c(
  "Raf", "Mek", "PLCg", "PIP2", "PIP3", "Erk", "Akt", "PKA", "PKC",
  "p38", "JNK"
)

get_dag <- function(which = c("consensus", "sachs"), int = c("none", "Akt", "PIP2", "Erk", "PKC", "PIP3")) {
  which <- match.arg(which)
  int <- match.arg(int)
  
  ### Consensus graph
  consensus <- matrix(c(
    0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0
  ), ncol = 11, byrow = TRUE)
  
  ### Sachs graph
  sachs <- matrix(c(
    0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0
  ), ncol = 11, byrow = TRUE)
  
  amat <- switch(which,
                 "consensus" = consensus,
                 "sachs" = sachs)
  colnames(amat) <- nms
  rownames(amat) <- nms
  
  ### Remove incoming nodes (if intervened)
  if (int != "none") {
    amat[, int] <- 0
  }
  
  amat
}

read_data <- function(int = c("none", "Akt", "PIP2", "Erk", "PKC", "PIP3"), path = DATA_PATH) {
  int <- match.arg(int)
  file <- switch(int,
                 "none" = "cd3cd28.xls",
                 "Akt"  = "cd3cd28+aktinhib.xls",
                 "PIP2" = "cd3cd28+psitect.xls",
                 "Erk"  = "cd3cd28+u0126.xls",
                 "PKC"  = "cd3cd28+g0076.xls",
                 "PIP3" = "cd3cd28+ly.xls")
  dat <- readxl::read_xls(file.path(path, file)) |>
    mutate_all(log)
  colnames(dat) <- nms
  dat
}

adj2dag <- function(adj_matrix) {
  nodes <- rownames(adj_matrix)
  dag_string <- "dag {"
  for (i in 1:nrow(adj_matrix)) {
    for (j in 1:ncol(adj_matrix)) {
      if (adj_matrix[i, j] == 1) {
        dag_string <- paste(dag_string, nodes[i], "->", nodes[j], ";")
      }
    }
  }
  dagitty(paste(dag_string, "}"))
}

# =========================
# 2) Agent: propose adjacency matrix (LLM)
# =========================
stop_if_no_key <- function() {
  if (!nzchar(Sys.getenv("OPENAI_API_KEY"))) {
    stop("Please set OPENAI_API_KEY via Sys.setenv(OPENAI_API_KEY='...')", call. = FALSE)
  }
}

agent_propose_adj <- function(variables, model = "gpt-4o-mini", temperature = 0) {
  stop_if_no_key()
  # Build prompts (strict JSON output)
  system_msg <- paste(
    "You output ONLY a JSON object with keys 'variables' and 'adjacency'.",
    "'variables' must exactly match the order provided.",
    "'adjacency' is a square binary matrix where entry (i,j)=1 means an edge variables[i] -> variables[j].",
    "The matrix must be a DAG: acyclic, with zero diagonal.")
  user_msg <- jsonlite::toJSON(list(
    instruction = "Return ONLY JSON. No comments.",
    variables = unname(variables),
    constraints = list(diagonal_zero = TRUE, binary_entries = TRUE, acyclic = TRUE)
  ), auto_unbox = TRUE)
  
  body <- list(
    model = model,
    temperature = temperature,
    messages = list(
      list(role = "system", content = system_msg),
      list(role = "user",   content = as.character(user_msg))
    )
  )
  
  resp <- request("https://api.openai.com/v1/chat/completions") |>
    req_headers(Authorization = paste("Bearer", Sys.getenv("OPENAI_API_KEY")),
                "Content-Type" = "application/json") |>
    req_body_json(body) |>
    req_perform()
  
  resp_json <- resp_body_json(resp, simplifyVector = TRUE)
  txt <- resp_json$choices[[1]]$message$content
  
  
  # Parse JSON
  out <- try(jsonlite::fromJSON(txt), silent = TRUE)
  if (inherits(out, "try-error")) stop("Agent did not return valid JSON. Raw output:\n", txt)
  if (!all(c("variables","adjacency") %in% names(out))) stop("Agent JSON missing required keys.")
  if (!identical(unname(out$variables), unname(variables))) stop("Agent variables do not match input order.")
  A <- as.matrix(out$adjacency)
  mode(A) <- "numeric"
  colnames(A) <- out$variables
  rownames(A) <- out$variables
  list(variables = out$variables, adjacency = A, raw = txt)
}

# =========================
# 3) Validation: adjacency matrix must be a DAG
# =========================
validate_adj <- function(A, variables) {
  p <- length(variables)
  if (!is.matrix(A)) stop("Adjacency is not a matrix.")
  if (!all(dim(A) == c(p,p))) stop("Adjacency has wrong shape.")
  if (!all(rownames(A) == variables) || !all(colnames(A) == variables)) {
    stop("Row/column names must exactly match variables (order too).")
  }
  if (any(is.na(A))) stop("Adjacency contains NA.")
  if (!all(A %in% c(0,1))) stop("Adjacency must be binary (0/1).")
  if (!all(diag(A) == 0)) stop("Adjacency diagonal must be all zeros (no self loops).")
  # Acyclicity via dagitty
  dag <- adj2dag(A)
  if (!isAcyclic(dag)) stop("Adjacency encodes cycles; not a DAG.")
  invisible(TRUE)
}

# =========================
# 4) LK/Kook testing (COMETS) on ONE dataset
# =========================
# Extract testable CIs exactly like in the shared paper code: keep |Z| > 0
extract_cis <- function(dag) {
  cis <- impliedConditionalIndependencies(dag)
  # Keep only those with non-empty conditioning set
  tcis <- cis[unlist(lapply(cis, function(x) length(x$Z) > 0))]
  tcis
}

# Run COMETS tests (GCM/PCM) + Holm correction, return tidy data.frame
run_lk_tests <- function(dag, data, tests = TESTS, alpha = 0.05) {
  tcis <- extract_cis(dag)
  if (length(tcis) == 0) {
    return(bind_rows(lapply(tests, function(tst) {
      tibble(test = tst, n_cis_tested = 0L, falsified = FALSE, CI = character(), p.value = numeric(), adj.p.value = numeric())
    })))
  }
  results <- lapply(tests, function(tst) {
    pv <- lapply(tcis, function(ci) {
      fm <- reformulate(
        paste0(paste0(ci$Y, collapse = "+"), "|", paste0(ci$Z, collapse = "+")),
        response = ci$X
      )
      suppressWarnings(comets(fm, data, test = tst, coin = TRUE)$p.value)
    }) |> unlist()
    adjp <- p.adjust(pv, method = "holm")
    tibble(
      test = tst,
      n_cis_tested = length(pv),
      falsified = any(adjp < alpha),
      CI = sapply(tcis, paste),
      p.value = pv,
      adj.p.value = adjp
    )
  })
  bind_rows(results)
}

# Convenience: run both Agent and Consensus on the same dataset
run_once <- function(data) {
  vars <- colnames(data)
  stopifnot(identical(unname(vars), unname(nms)))
  
  # Agent proposal
  message("Requesting agent adjacency matrix...")
  agent_out <- agent_propose_adj(vars, model = "gpt-5-thinking", temperature = 0)
  A_agent <- agent_out$adjacency
  # Validate (retry once if needed)
  ok <- FALSE
  tryCatch({ validate_adj(A_agent, vars); ok <<- TRUE }, error = function(e) message("Validation failed: ", e$message))
  if (!ok) {
    message("Retrying agent once with error hint...")
    # Re-prompt with error hint (lightweight)
    # NOTE: For simplicity we just call again; a more advanced retry would include the error message in the prompt.
    agent_out <- agent_propose_adj(vars, model = "gpt-5-thinking", temperature = 0)
    A_agent <- agent_out$adjacency
    validate_adj(A_agent, vars) # if fails -> error
  }
  
  dag_agent <- adj2dag(A_agent)
  
  # LK tests (Agent)
  message("Running LK tests for Agent DAG...")
  res_agent <- run_lk_tests(dag_agent, data, tests = TESTS, alpha = alpha)
  res_agent$graph <- "Agent"
  
  # Consensus baseline
  message("Running LK tests for Consensus DAG...")
  A_cons <- get_dag("consensus", int = "none")
  dag_cons <- adj2dag(A_cons)
  res_cons <- run_lk_tests(dag_cons, data, tests = TESTS, alpha = alpha)
  res_cons$graph <- "Consensus"
  
  # Bind tidy results
  out <- bind_rows(res_agent, res_cons) |>
    relocate(graph, test, n_cis_tested, falsified, CI, p.value, adj.p.value)
  
  # Save artifacts
  if (SAVE_RESULTS) {
    if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)
    write_csv(out, file.path(RESULTS_DIR, "one_shot_agent_vs_consensus.csv"))
    writeLines(agent_out$raw, file.path(RESULTS_DIR, "agent_raw.json"))
  }
  
  list(
    adjacency_agent = A_agent,
    dag_agent = dag_agent,
    results = out
  )
}

# =========================
# 5) RUN (single observational dataset)
# =========================
# Load the observational dataset (no intervention)
D_obs <- read_data(int = "none", path = DATA_PATH)

# Execute once: Agent vs Consensus on this dataset
run_out <- run_once(D_obs)

# Print a compact summary to console
summary_tbl <- run_out$results |>
  group_by(graph, test) |>
  summarise(n_cis_tested = unique(n_cis_tested), falsified = unique(falsified),
            n_violations = sum(adj.p.value < alpha), .groups = "drop")
print(summary_tbl)

# End of script
