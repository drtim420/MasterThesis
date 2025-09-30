# agent_llm_controller

# LLM-Directed Agent Controller using Tool Descriptions
# ------------------------------------------------------
# This file turns the previously defined tools into an LLM-controlled agent.
# The LLM chooses which tool to call next based on the goal/state, via a
# constrained JSON action format. No raw data is ever exposed to the LLM.
# ------------------------------------------------------

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
# We expose high-level tools. Each entry declares: name, description, arguments.
# The controller will execute the corresponding R function and return a JSON
# result back to the LLM as a tool message.

# NOTE: The underlying implementations live in tools_agent_LK.R

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
    description = "Run LK/COMETS GCM and PCM tests on the dataset using the CURRENT adjacency. Returns per-test falsification and a list of violated CIs (with adjusted p-values).",
    args = list(alpha = list(type = "number", default = 0.05))
  ),
  list(
    name = "suggest_edits",
    description = "Propose up to K minimal edits (add/remove/reverse) to resolve the CURRENT set of violated CIs.",
    args = list(k = list(type = "integer", default = 3))
  ),
  list(
    name = "apply_edits",
    description = "Apply a list of edits to the CURRENT adjacency and revalidate.",
    args = list(edits = list(type = "array[object]", schema = list(fields = c("op","from","to"))))
  ),
  list(
    name = "finish",
    description = "Finish the task and produce a short summary. Include whether the DAG is falsified by GCM/PCM, number of CIs tested, and remaining issues if any.",
    args = list(summary = list(type = "string"))
  )
)

# Helper: render a compact tool card for the system prompt
.tool_card <- function(s) {
  paste0("- ", s$name, ": ", s$description,
         "\n  args: ", toJSON(s$args, auto_unbox = TRUE))
}

# --------------------------- Planner Prompt ----------------------------------
planner_system_prompt <- function() {
  paste(
    "You are an autonomous causal-structure agent. You have a limited set of tools.",
    "Your job: propose a DAG over the given variables, test it with LK/COMETS (GCM & PCM),",
    "and if falsified, suggest up to K minimal edits to address the rejected CIs,",
    "while always keeping the graph a valid DAG (acyclic, 0/1 adjacency, zero diagonal).",
    "NEVER access raw data; tests are executed by the environment and summarized for you.",
    "Always respond in STRICT JSON with one of these forms:",
    '{"call": {"tool": "<tool_name>", "args": { ... }}}  OR  ',
    '{"final": {"summary": "<short result>"}}',
    "Do not output any extra text.",
    "\nTools:\n", paste(vapply(.tool_specs, .tool_card, character(1)), collapse = "\n"),
    sep = "\n")
}

# --------------------------- Dispatcher --------------------------------------
# state: list with fields {data, variables, A (current adjacency), last_results, violations}

.dispatch_tool <- function(tool, args, state) {
  if (tool == "propose_matrix") {
    stop_if_no_key()
    res <- tool_propose_matrix(state$variables)
    state$A <- res$adjacency
    list(state = state, result = list(status = "ok", edges = which(state$A==1, arr.ind=TRUE)))
  } else if (tool == "validate_matrix") {
    v <- tool_validate_matrix(state$A, state$variables)
    list(state = state, result = list(ok = v$ok, errors = v$errors))
  } else if (tool == "lk_test") {
    alpha <- args$alpha %||% 0.05
    v <- tool_validate_matrix(state$A, state$variables)
    if (!v$ok) return(list(state = state, result = list(error = paste(v$errors, collapse=", "))))
    cis <- tool_extract_cis(v$dag)
    gcm <- tool_lk_test(state$data, cis, test = "gcm", alpha = alpha)
    pcm <- tool_lk_test(state$data, cis, test = "pcm", alpha = alpha)
    res <- dplyr::bind_rows(gcm, pcm)
    state$last_results <- res
    state$violations <- dplyr::filter(res, adj.p.value < alpha)
    # compact return for the LLM
    compact <- res |>
      dplyr::group_by(test) |>
      dplyr::summarise(n_cis_tested = unique(n_cis_tested),
                       falsified = any(adj.p.value < alpha),
                       n_violations = sum(adj.p.value < alpha), .groups = "drop")
    viol_top <- state$violations |>
      dplyr::arrange(adj.p.value) |>
      dplyr::mutate(adj.p.value = round(adj.p.value, 6)) |>
      dplyr::select(test, CI, adj.p.value) |>
      head(10)
    list(state = state, result = list(summary = compact, violations = viol_top))
  } else if (tool == "suggest_edits") {
    k <- args$k %||% 3
    if (is.null(state$violations) || nrow(state$violations) == 0) {
      return(list(state = state, result = list(edits = list(), note = "No violations to fix.")))
    }
    edits <- tool_suggest_edits(state$A, state$variables, state$violations, budget_k = k)
    list(state = state, result = list(edits = edits))
  } else if (tool == "apply_edits") {
    edits <- args$edits
    if (is.null(edits) || nrow(edits) == 0) return(list(state = state, result = list(status = "no_edits")))
    A2 <- tool_apply_edits(state$A, state$variables, edits)
    state$A <- A2
    list(state = state, result = list(status = "ok"))
  } else if (tool == "finish") {
    return(list(state = state, result = list(done = TRUE, summary = args$summary)))
  } else {
    list(state = state, result = list(error = paste0("Unknown tool: ", tool)))
  }
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# --------------------------- Agent Runner ------------------------------------
run_llm_agent <- function(data, variables, max_steps = 8, alpha = 0.05, verbose = TRUE) {
  stop_if_no_key()
  state <- list(data = data, variables = variables, A = NULL, last_results = NULL, violations = NULL)
  
  sys <- planner_system_prompt()
  msgs <- list(list(role = "system", content = sys),
               list(role = "user", content = toJSON(list(goal = "Propose and validate a DAG over the variables using LK/COMETS; refine only if falsified.",
                                                         variables = variables,
                                                         alpha = alpha), auto_unbox = TRUE)))
  
  step <- 1
  while (step <= max_steps) {
    # Call LLM planner
    body <- list(model = MODEL_PLANNER, temperature = 0,
                 messages = msgs)
    req <- request("https://api.openai.com/v1/chat/completions") |>
      req_headers(Authorization = paste("Bearer", Sys.getenv("OPENAI_API_KEY")),
                  "Content-Type" = "application/json") |>
      req_body_json(body)
    resp <- req_perform(req)
    resp_json <- resp_body_json(resp, simplifyVector = FALSE)
    if (!is.null(resp_json$error)) stop(paste0("OpenAI API error: ", resp_json$error$message))
    content <- resp_json$choices[[1]]$message$content
    if (is.list(content)) content <- paste(vapply(content, function(p) p$text %||% "", character(1)), collapse = "")
    
    # Expect strict JSON with either {call:{tool,...}} or {final:{...}}
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
    
    # Add tool result back into conversation
    msgs <- c(msgs,
              list(list(role = "assistant", content = as.character(content))),
              list(list(role = "tool", content = tool_result_json, name = tool)))
    
    # If after a test no violations remain, the planner should choose finish next.
    step <- step + 1
  }
  
  # Return final objects for downstream use
  list(adjacency = state$A, last_results = state$last_results, violations = state$violations)
}

# --------------------------- Example usage -----------------------------------
# source("R/tools_agent_LK.R")
# D_obs <- read_data(int = "none", path = "../data")
# out <- run_llm_agent(D_obs, variables = colnames(D_obs), max_steps = 8, alpha = 0.05)
# out$last_results %>% dplyr::group_by(test) %>% dplyr::summarise(n_cis_tested = unique(n_cis_tested),
#                                                                falsified = any(adj.p.value < 0.05),
#                                                                n_violations = sum(adj.p.value < 0.05))
