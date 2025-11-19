
# Minimal LLM-assisted agent for DAG–data consistency checks
#
# Given:
#   - a DAG in adjacency matrix
#   - a dataset 
#
# The agent:
#   1) Uses the existing run_ci_tests() function to test all implied
#      conditional independencies (CIs) of the DAG against the data.
#   2) Summarises how many CIs were tested and how many were rejected





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




# DAG–CI agent with LLM summary 

dag_ci_agent <- function(dat, amat,
                         alpha = 0.05,
                         tests = c("gcm", "pcm"),
                         graph_name = "Hypothesized DAG") {
  
  # existing CI tests
  res <- run_ci_tests(amat, dat, tests = tests, alpha = alpha)
  
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



# Prepare data & DAG (from your existing code)
dat  <- read_sachs_observational("../data")
amat <- get_sachs_amat()


# Run minimal agent
out <- dag_ci_agent(dat, amat, alpha = 0.05, graph_name = "Sachs DAG")

cat("\n--- LLM Interpretation ---\n", out$interpretation, "\n")



