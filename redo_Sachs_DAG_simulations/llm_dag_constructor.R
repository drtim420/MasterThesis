




# OpenAI helper 

.have_key <- function(var = "OPENAI_API_KEY") {
  nzchar(Sys.getenv(var, unset = ""))
}

openai_chat_min <- function(system_msg, user_msg,
                            model = Sys.getenv("OPENAI_MODEL", unset = "gpt-4o-mini"),
                            temperature = 0) {
  if (!.have_key()) return(NULL)
  
  body <- list(
    model = model,
    temperature = temperature,
    messages = list(
      list(role = "system", content = system_msg),
      list(role = "user",   content = user_msg)
    )
  )
  
  req <- httr2::request("https://api.openai.com/v1/chat/completions") |>
    httr2::req_headers(
      Authorization = paste("Bearer", Sys.getenv("OPENAI_API_KEY")),
      "Content-Type" = "application/json"
    ) |>
    httr2::req_body_json(body)
  
  resp <- httr2::req_perform(req)
  j <- httr2::resp_body_json(resp, simplifyVector = FALSE)
  
  if (!is.null(j$error)) {
    stop(paste0("OpenAI API error: ", j$error$message))
  }
  
  content <- j$choices[[1]]$message$content
  if (is.list(content)) {
    # handle "content blocks" if API ever returns them
    content <- paste(vapply(content, function(part) {
      if (!is.null(part$text)) part$text else if (!is.null(part$content)) part$content else ""
    }, character(1)), collapse = "")
  }
  as.character(content)
}



# ---- LLM-based DAG proposal ----------------------------------------------

# Ask an LLM to propose a DAG given only the variable names.
# Output: adjacency matrix amat (same format as get_sachs_amat()).
propose_dag_from_llm <- function(vars,
                                 graph_name = "LLM DAG",
                                 model = Sys.getenv("OPENAI_MODEL", unset = "gpt-4o-mini"),
                                 temperature = 0.2) {
  if (!.have_key()) {
    message("No OPENAI_API_KEY set – falling back to empty DAG.")
    A <- matrix(0L, nrow = length(vars), ncol = length(vars),
                dimnames = list(vars, vars))
    return(A)
  }
  
  sys_msg <- paste(
    "You are an expert in causal discovery and cell signaling.",
    "You will be given a list of variable names from a biological dataset.",
    "Your task is to propose a plausible directed acyclic graph (DAG)",
    "representing causal relations between these variables.",
    "",
    "IMPORTANT OUTPUT FORMAT:",
    "- Output ONLY directed edges, one per line.",
    "- Each line must be: Parent,Child",
    "- Use ONLY variable names from the provided list.",
    "- Do NOT include any headers, explanations, or extra text.",
    "- The edges must define a DAG (no directed cycles)."
  )
  
  usr_msg <- paste0(
    "Variable names:\n",
    paste(vars, collapse = ", "),
    "\n\n",
    "Propose a plausible DAG for these variables in the exact format:\n",
    "Parent,Child\nParent,Child\n...\n"
  )
  
  txt <- openai_chat_min(sys_msg, usr_msg, model = model, temperature = temperature)
  if (is.null(txt) || !nzchar(txt)) {
    message("LLM call failed or returned empty text – falling back to empty DAG.")
    A <- matrix(0L, nrow = length(vars), ncol = length(vars),
                dimnames = list(vars, vars))
    return(A)
  }
  
  # Parse "Parent,Child" lines
  lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  
  A <- matrix(0L, nrow = length(vars), ncol = length(vars),
              dimnames = list(vars, vars))
  
  for (ln in lines) {
    parts <- strsplit(ln, ",", fixed = TRUE)[[1]]
    if (length(parts) != 2) next
    parent <- trimws(parts[1])
    child  <- trimws(parts[2])
    
    if (!(parent %in% vars) || !(child %in% vars)) next
    if (parent == child) next
    
    A[parent, child] <- 1L
  }
  
  # Optional: check acyclicity; if cyclic, warn but still return.
  g <- igraph::graph_from_adjacency_matrix(A, mode = "directed")
  if (!igraph::is_dag(g)) {
    warning("LLM-proposed adjacency is not a DAG (contains cycles). ",
            "You may need to adjust the prompt or post-process the edges.")
  }
  
  message("LLM DAG proposal: ", sum(A), " directed edges created.")
  A
}








# DAG–CI agent with LLM summary 

dag_ci_agent <- function(dat, amat_llm,
                         alpha = 0.05,
                         tests = c("gcm", "pcm"),
                         graph_name = "Hypothesized DAG") {
  
  # existing CI tests
  res <- run_ci_tests(amat_llm, dat, tests = tests, alpha = alpha)
  
  # Split into summary + rejected list
  summary_tbl <- res |>
    dplyr::group_by(test) |>
    dplyr::summarise(
      n_tests    = dplyr::n(),
      n_rejected = sum(rejected),
      min_adj_p  = min(adj.p.value, na.rm = TRUE),
      .groups    = "drop"
    )
  
  rejected_tbl <- res |>
    dplyr::filter(rejected) |>
    dplyr::arrange(adj.p.value) |>
    dplyr::select(test, CI, adj.p.value)
  
  # Turn tables into plain text for the prompt
  smry_txt <- paste(capture.output(print(summary_tbl, n = Inf)), collapse = "\n")
  rej_txt  <- if (nrow(rejected_tbl)) {
    paste(capture.output(print(rejected_tbl, n = min(10, nrow(rejected_tbl)))), collapse = "\n")
  } else {
    "<none>"
  }
  
  # LLM prompt
  sys_msg <- paste(
    "You are a statistician explaining CI-tests for DAGs.",
    "Output 4–6 short bullet points, no prose paragraphs.",
    "Define 'falsified' := CI with Holm-adjusted p < alpha."
  )
  
  usr_msg <- paste0(
    "We tested whether the DAG '", graph_name, "' is consistent with data.\n",
    sprintf("alpha = %.3f\n\n", alpha),
    "Per-test summary (per CI test type):\n",
    smry_txt, "\n\n",
    "Rejected CI statements (Holm-adjusted p-values):\n",
    rej_txt, "\n\n",
    "Explain briefly:\n",
    "1) Does the DAG pass or fail overall?\n",
    "2) How many CIs were tested and rejected per test type?\n",
    "3) Mention the most important 1–3 rejected CIs, if any.\n",
    "4) Remind that 'not falsified' ≠ 'true graph'.\n"
  )
  
  interp <- openai_chat_min(sys_msg, usr_msg)
  
  
  # Return everything
  list(
    raw_results   = res,
    summary       = summary_tbl,
    rejected      = rejected_tbl,
    interpretation = interp
  )
}


# Prepare data
dat  <- read_sachs_observational("../data")

# Get a hypothesized DAG from the LLM using only the column names
vars <- colnames(dat)
amat_llm <- propose_dag_from_llm(vars, graph_name = "LLM Sachs-like DAG")

# Now run exactly the same pipeline as before, just with amat_llm
out_llm <- dag_ci_agent(dat, amat_llm, alpha = 0.05, graph_name = "LLM Sachs DAG")

cat("\n--- LLM-based DAG: CI Interpretation ---\n", out_llm$interpretation, "\n")



########### 
###plot the dag

# LLM DAG object
dag_llm <- adj2dag(amat_llm)

# Same order / layout as the Sachs DAG
paper_order <- c("PIP2","PLCg","Mek","Raf","JNK","p38",
                 "PKC","PKA","Akt","Erk","PIP3")

# Plot LLM DAG with same circular layout (PIP2 at 12 o'clock, mirrored)
plot_sachs_paper_circle(
  dag    = dag_llm,
  order  = paper_order,
  anchor = "PIP2",
  mirror_y = TRUE,
  title  = "LLM–proposed Sachs DAG"
)
