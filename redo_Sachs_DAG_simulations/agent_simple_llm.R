# ============================================================
# agent_simple_llm.R
# Minimal agent:
#  - Input: data.frame (log-transformed Sachs vars), hypothesized adjacency (amat)
#  - Check CI consistency (COMETS GCM/PCM + Holm) using your 01_* functions
#  - If falsified: list rejected CIs + suggest likely missing edges
#  - Visual: DAG + red dashed chords for rejected pairs
#  - LLM: ≤6 bullet interpretation (skips if OPENAI_API_KEY absent)
# ============================================================


source("data_simulation_v2.R")

# Load .env if present
if (!requireNamespace("dotenv", quietly = TRUE)) install.packages("dotenv")
dotenv::load_dot_env(".env")  # requires a file named .env in getwd()

ensure_openai_key <- function(var = "OPENAI_API_KEY") {
  val <- Sys.getenv(var, unset = "")
  if (!nzchar(val)) stop(sprintf("%s is not set. Put it in .env or ~/.Renviron as:\n%s=sk-...\n", var, var), call. = FALSE)
  invisible(val)
}

.have_key <- function() nzchar(Sys.getenv("OPENAI_API_KEY"))


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(httr2)
})

# ---- Tiny OpenAI helper (optional) ------------------------------------------
.have_key <- function() nzchar(Sys.getenv("OPENAI_API_KEY"))

openai_chat_min <- function(system_msg, user_msg,
                            model = Sys.getenv("OPENAI_MODEL", unset = "gpt-4o-mini"),
                            temperature = 0) {
  if (!.have_key()) return(NULL)
  body <- list(
    model = model, temperature = temperature,
    messages = list(
      list(role = "system", content = system_msg),
      list(role = "user",   content = user_msg)
    )
  )
  req <- request("https://api.openai.com/v1/chat/completions") |>
    req_headers(
      Authorization = paste("Bearer", Sys.getenv("OPENAI_API_KEY")),
      "Content-Type" = "application/json"
    ) |>
    req_body_json(body)
  resp <- req_perform(req)
  j <- resp_body_json(resp, simplifyVector = FALSE)
  if (!is.null(j$error)) stop(paste0("OpenAI API error: ", j$error$message))
  content <- j$choices[[1]]$message$content
  if (is.list(content)) content <- paste(vapply(content, function(part) {
    if (!is.null(part$text)) part$text else if (!is.null(part$content)) part$content else ""
  }, character(1)), collapse = "")
  as.character(content)
}

# ---- Adapter: glue your 01_* functions into one agent call -------------------
# Requires these from 01_test_vs_sim_sachs.R:
#   adj2dag(), run_ci_tests(), suggest_and_test_orientations(),
#   ensure_circular_coords(), plot_dag_with_tested_cis()

agent_check <- function(data, amat_hyp, out_dir = "../results",
                        alpha = 0.05, tests = c("gcm","pcm"),
                        circle_order = c("PIP2","PLCg","Mek","Raf","JNK","p38","PKC","PKA","Akt","Erk","PIP3"),
                        graph_name = "Hypothesized DAG") {
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # 1) CI tests (LK/Kook: only CIs with non-empty Z)
  res <- run_ci_tests(amat_hyp, data, tests = tests, alpha = alpha)
  rejected <- res |> dplyr::filter(rejected) |> dplyr::arrange(adj.p.value)
  
  # compact per-test summary
  sum_tbl <- res |>
    group_by(test) |>
    summarise(
      n_tests = dplyr::n(),
      n_rejected = sum(adj.p.value < alpha),
      .groups = "drop"
    )
  
  # 2) Suggestions: rank pairs from rejected CIs, test feasible orientations
  sg <- suggest_and_test_orientations(
    viol_tbl = rejected,
    amat0    = amat_hyp,
    dat      = data,
    tests    = tests,
    alpha    = alpha
  )
  sugg <- sg$recommendation |>
    tidyr::separate(pair, into = c("X","Y"), sep = "—", remove = FALSE)
  
  # 3) Plot: show only rejected pairs as red dashed chords
  # Parse CI strings -> ci objects
  ci_from_str <- function(s) {
    parts <- strsplit(s, " _\\|\\|_ ", perl = TRUE)[[1]]
    X <- trimws(parts[1])
    yz <- strsplit(parts[2], "\\|", perl = TRUE)[[1]]
    Y <- trimws(yz[1])
    Zs <- trimws(strsplit(yz[2], ",")[[1]])
    Zs <- Zs[Zs != ""]
    if (identical(Zs, "∅")) Zs <- character(0)
    list(X = X, Y = Y, Z = Zs)
  }
  ci_list_rej <- lapply(rejected$CI, ci_from_str)
  
  dag_hyp <- adj2dag(amat_hyp)
  dag_circ <- ensure_circular_coords(dag_hyp, order = circle_order, anchor = "PIP2", mirror_y = TRUE)
  p <- plot_dag_with_tested_cis(
    dag_circ,
    cis_list = ci_list_rej,
    title = sprintf("%s — Rejected CI chords (Holm α=%.2f): %d",
                    graph_name, alpha, nrow(rejected))
  )
  
  # Save artifacts
  f_rej  <- file.path(out_dir, "agent_rejected_cis.csv")
  f_sugg <- file.path(out_dir, "agent_suggested_edges.csv")
  f_png  <- file.path(out_dir, "agent_rejected_chords.png")
  readr::write_csv(rejected, f_rej)
  readr::write_csv(sugg,    f_sugg)
  ggplot2::ggsave(f_png, p, width = 7, height = 7, dpi = 150)
  
  list(
    eval = list(summary = sum_tbl, rejected = rejected),
    suggestions = sugg,
    plot = p,
    files = list(rejected_csv = f_rej, suggestions_csv = f_sugg, plot_png = f_png),
    graph_name = graph_name,
    alpha = alpha
  )
}

# ---- One-shot agent with LLM interpretation ----------------------------------
run_agent_simple_llm <- function(data, amat_hyp, out_dir = "../results",
                                 alpha = 0.05, tests = c("gcm","pcm"),
                                 graph_name = "Hypothesized DAG") {
  ag <- agent_check(data, amat_hyp, out_dir = out_dir, alpha = alpha, tests = tests, graph_name = graph_name)
  
  smry <- ag$eval$summary |> mutate(across(everything(), as.character))
  rej  <- ag$eval$rejected |> arrange(adj.p.value) |> select(test, CI, adj.p.value)
  sugg <- ag$suggestions
  
  smry_txt <- paste(capture.output(print(smry, n = Inf)), collapse = "\n")
  rej_txt  <- if (nrow(rej))  paste(capture.output(print(head(rej, 25),  n = Inf)), collapse = "\n") else "<none>"
  sug_txt  <- if (nrow(sugg)) paste(capture.output(print(head(sugg, 25), n = Inf)), collapse = "\n") else "<none>"
  
  sys <- paste(
    "You output at most 6 short bullet points.",
    "Audience: stats-savvy but busy.",
    "Define 'falsified' as: any Holm-adjusted p-value < alpha."
  )
  usr <- paste0(
    "Summarize CI-consistency for ", graph_name, ". Use bullets only.\n",
    sprintf("alpha = %.3f\n\n", ag$alpha),
    "Per-test summary:\n", smry_txt, "\n\n",
    "Top rejected CI statements (adj. p):\n", rej_txt, "\n\n",
    "Candidate missing edges (ranked):\n", sug_txt, "\n\n",
    "Instructions: (1) Say pass/fail. (2) #CIs tested and #rejected per test. ",
    "(3) Name 2–4 critical pairs. (4) Mention feasible orientations if any. ",
    "(5) Caveat that 'not falsified' ≠ 'true graph'."
  )
  
  interp <- openai_chat_min(sys, usr)
  if (is.null(interp)) {
    # fallback bullets if no key
    bullets <- c(
      "• LLM not configured (OPENAI_API_KEY missing).",
      paste0("• ", graph_name, ": ", paste(smry$test, " n=", smry$n_tests, ", rej=", smry$n_rejected, collapse = " | ")),
      if (nrow(rej)) paste0("• Top rejection: ", rej$CI[1], " (", rej$test[1], "), adj.p=", signif(rej$adj.p.value[1], 3)) else "• No rejections.",
      "• See CSVs and plot in results folder."
    )
    interp <- paste(bullets, collapse = "\n")
  }
  
  f_md <- file.path(out_dir, "agent_interpretation.md")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  writeLines(interp, f_md)
  
  ag$files$interpretation_md <- f_md
  ag$interpretation <- interp
  ag
}



dat  <- read_sachs_observational("../data")
amat <- get_sachs_amat()


out <- run_agent_simple_llm(dat, amat, out_dir = "../results", alpha = 0.05)
cat("\n--- LLM Interpretation ---\n", out$interpretation, "\n")
print(out$plot)
