# agent_llm_controller

# LLM-Directed Agent Controller using Tool Descriptions (no edit tools)
# ---------------------------------------------------------------------
# The LLM chooses which tool to call next via STRICT-JSON actions:
#   {"call":{"tool":"<name>","args":{...}}}  OR  {"final":{"summary":"..."}}
# Tool outputs are fed back as assistant messages containing:
#   {"tool_result":{"tool":"<name>","result":{...}}}
# No raw data is ever exposed to the LLM.

suppressPackageStartupMessages({
  library(jsonlite)
  library(httr2)
})

# (Optional) auto-load .env (API key & model names)
if (requireNamespace("dotenv", quietly = TRUE)) {
  try(dotenv::load_dot_env(".env"), silent = TRUE)
}

MODEL_PLANNER <- Sys.getenv("OPENAI_PLANNER_MODEL", unset = Sys.getenv("OPENAI_MODEL", unset = "gpt-4o-mini"))

stop_if_no_key <- function() {
  if (!nzchar(Sys.getenv("OPENAI_API_KEY"))) stop("OPENAI_API_KEY not set.")
}

# --------------------------- Tool Registry (for the LLM) ---------------------
# Underlying implementations live in tools_agent_LK.R
.tool_specs <- list(
  list(
    name = "propose_matrix",
    description = "Propose a directed adjacency matrix (DAG) from variable names only.",
    args = list(variables = list(type = "array[string]"))
  ),
  list(
    name = "validate_matrix",
    description = "Validate current adjacency: binary, diag=0, acyclic. Returns ok/errors.",
    args = list()
  ),
  list(
    name = "lk_test",
    description = "Run LK/COMETS GCM and PCM tests using the CURRENT adjacency. Returns per-test falsification and a list of violated CIs (with adjusted p-values).",
    args = list(alpha = list(type = "number", default = 0.05))
  ),
  list(
    name = "finish",
    description = "Finish the task and produce a short summary (mention falsification status, #CIs tested, and any issues).",
    args = list(summary = list(type = "string"))
  )
)

# Helper: render tool card for the system prompt
.tool_card <- function(s) {
  paste0("- ", s$name, ": ", s$description,
         "\n  args: ", toJSON(s$args, auto_unbox = TRUE))
}

# --------------------------- Planner Prompt ----------------------------------
planner_system_prompt <- function() {
  paste(
    "You are an autonomous causal-structure agent. You have a limited set of tools.",
    "Your job: propose a DAG over the given variables, validate it, run LK/COMETS (GCM & PCM), then finish.",
    "Always keep the graph a valid DAG (acyclic, 0/1 adjacency, zero diagonal).",
    "NEVER access raw data; tests are executed by the environment and summarized for you.",
    "Respond ONLY in STRICT JSON as ONE of:",
    '{"call":{"tool":"<tool_name>","args":{...}}}',
    "OR",
    '{"final":{"summary":"<short result>"}}',
    "The environment returns tool outputs as assistant messages whose ENTIRE content is:",
    '{"tool_result":{"tool":"<name>","result":{...}}}. Read that JSON and decide your next action.',
    "",
    "Tools:\n", paste(vapply(.tool_specs, .tool_card, character(1)), collapse = "\n"),
    sep = "\n"
  )
}

# --------------------------- Dispatcher --------------------------------------
# state: list with {data, variables, A (current adjacency), last_results, violations}
.dispatch_tool <- function(tool, args, state) {
  if (tool == "propose_matrix") {
    stop_if_no_key()
    res <- tool_propose_matrix(state$variables)  # from tools_agent_LK.R
    state$A <- res$adjacency
    return(list(state = state,
                result = list(status = "ok",
                              edges = which(state$A == 1, arr.ind = TRUE))))
  }
  
  if (tool == "validate_matrix") {
    v <- tool_validate_matrix(state$A, state$variables)
    return(list(state = state, result = list(ok = v$ok, errors = v$errors)))
  }
  
  if (tool == "lk_test") {
    alpha <- args$alpha %||% 0.05
    v <- tool_validate_matrix(state$A, state$variables)
    if (!v$ok) return(list(state = state, result = list(error = paste(v$errors, collapse = ", "))))
    cis <- tool_extract_cis(v$dag)
    gcm <- tool_lk_test(state$data, cis, test = "gcm", alpha = alpha)
    pcm <- tool_lk_test(state$data, cis, test = "pcm", alpha = alpha)
    res <- dplyr::bind_rows(gcm, pcm)
    state$last_results <- res
    state$violations  <- dplyr::filter(res, adj.p.value < alpha)
    
    # compact summary for planner
    compact <- res |>
      dplyr::group_by(test) |>
      dplyr::summarise(n_cis_tested = unique(n_cis_tested),
                       falsified = any(adj.p.value < alpha),
                       n_violations = sum(adj.p.value < alpha),
                       .groups = "drop")
    viol_top <- state$violations |>
      dplyr::arrange(adj.p.value) |>
      dplyr::mutate(adj.p.value = round(adj.p.value, 6)) |>
      dplyr::select(test, CI, adj.p.value) |>
      head(10)
    
    return(list(state = state, result = list(summary = compact, violations = viol_top)))
  }
  
  if (tool == "finish") {
    return(list(state = state, result = list(done = TRUE, summary = args$summary)))
  }
  
  list(state = state, result = list(error = paste0("Unknown tool: ", tool)))
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# --------------------------- Agent Runner ------------------------------------
run_llm_agent <- function(data, variables, max_steps = 8, alpha = 0.05, verbose = TRUE) {
  stop_if_no_key()
  state <- list(data = data, variables = variables, A = NULL, last_results = NULL, violations = NULL)
  
  # Initialize conversation
  msgs <- list(
    list(role = "system", content = planner_system_prompt()),
    list(role = "user", content = toJSON(list(
      goal = "Propose and validate a DAG over the variables using LK/COMETS; finish when done.",
      variables = variables,
      alpha = alpha
    ), auto_unbox = TRUE))
  )
  
  step <- 1
  while (step <= max_steps) {
    # Call LLM planner
    body <- list(model = MODEL_PLANNER, temperature = 0, messages = msgs)
    
    req <- request("https://api.openai.com/v1/chat/completions") |>
      req_headers(Authorization = paste("Bearer", Sys.getenv("OPENAI_API_KEY")),
                  "Content-Type" = "application/json") |>
      req_body_json(body)
    
    resp <- tryCatch(
      req_perform(req),
      error = function(e) {
        msg <- conditionMessage(e)
        rb <- tryCatch(resp_body_string(e$response), error = function(...) "")
        stop("Planner HTTP error: ", msg, "\n", rb)
      }
    )
    
    resp_json <- resp_body_json(resp, simplifyVector = FALSE)
    if (!is.null(resp_json$error)) stop(paste0("OpenAI API error: ", resp_json$error$message))
    
    content <- resp_json$choices[[1]]$message$content
    if (is.list(content)) content <- paste(vapply(content, function(p) p$text %||% "", character(1)), collapse = "")
    
    # Expect strict JSON with either {call:{...}} or {final:{...}}
    act <- try(fromJSON(as.character(content)), silent = TRUE)
    if (inherits(act, "try-error")) stop("Planner returned non-JSON:\n", as.character(content))
    
    if (!is.null(act$final)) {
      if (verbose) cat("\n[FINAL SUMMARY]\n", act$final$summary, "\n")
      break
    }
    
    if (is.null(act$call$tool)) stop("Planner JSON missing call.tool")
    
    tool <- act$call$tool
    args <- act$call$args %||% list()
    if (verbose) cat(sprintf("\n[STEP %d] Calling tool: %s\n", step, tool))
    
    # Dispatch
    res <- .dispatch_tool(tool, args, state)
    state <- res$state
    tool_result_json <- toJSON(res$result, auto_unbox = TRUE, null = "null")
    
    # Append: echo planner action, then the tool result (as assistant JSON blob)
    msgs <- c(
      msgs,
      list(list(role = "assistant", content = as.character(content))),
      list(list(role = "assistant",
                content = paste0('{"tool_result":{"tool":"', tool, '","result":', tool_result_json, '}}')))
    )
    
    step <- step + 1
  }
  
  # Return final objects for downstream use
  list(adjacency = state$A, last_results = state$last_results, violations = state$violations)
}

# --------------------------- Example usage -----------------------------------
# source("R/tools_agent_LK.R")
# D_obs <- read_data(int = "none", path = "../data")
# out <- run_llm_agent(D_obs, variables = colnames(D_obs), max_steps = 8, alpha = 0.05)
# out$last_results %>% dplyr::group_by(test) %>% dplyr::summarise(
#   n_cis_tested = unique(n_cis_tested),
#   falsified = any(adj.p.value < 0.05),
#   n_violations = sum(adj.p.value < 0.05)
# )
