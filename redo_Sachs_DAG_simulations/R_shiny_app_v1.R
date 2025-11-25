# app.R
# DAG–CI Agent with LLM-based DAG proposal and Shiny UI

## ---- libraries ----
library(shiny)
library(dplyr)
library(igraph)
library(httr2)
library(visNetwork)
library(comets)
library(dagitty)
library(shinybusy)
## ======================================================
## 1) OpenAI helper
## ======================================================

library(dotenv)

# Load .env from project root
if (file.exists(".env")) {
  dotenv::load_dot_env(".env")
}


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

 

circle_layout <- function(vars, order, radius = 300) {
  angle_step <- 2*pi / length(order)
  layout <- data.frame(
    id = order,
    x  = sapply(seq_along(order), function(i) radius * sin((i-1) * angle_step)),
    y  = sapply(seq_along(order), function(i) radius * cos((i-1) * angle_step))
  )
  layout <- layout[layout$id %in% vars, ]
  layout
}


make_safe_names <- function(x) {
  x2 <- gsub("[^0-9A-Za-z_\\.]", "_", x)   # replace bad chars by "_"
  # ensure starts with a letter
  x2 <- ifelse(grepl("^[A-Za-z]", x2),
               x2,
               paste0("V", x2))
  make.unique(x2)
}


make_circle_nodes <- function(vars,
                              paper_order = NULL,
                              radius = 300) {
  # If a preferred order is given (e.g. Sachs), use that first
  if (!is.null(paper_order)) {
    vars_ordered <- c(
      paper_order[paper_order %in% vars],
      vars[!vars %in% paper_order]
    )
  } else {
    vars_ordered <- vars
  }
  
  n <- length(vars_ordered)
  # angles: 0 at 12 o'clock, then clockwise
  angles <- seq(0, 2*pi - 2*pi/n, length.out = n)
  
  coords <- data.frame(
    id = vars_ordered,
    x  = radius * sin(angles),
    y  = radius * cos(angles),
    stringsAsFactors = FALSE
  )
  
  # Merge back to original vars: if some vars weren't in paper_order,
  # they still get positions.
  nodes <- merge(
    data.frame(id = vars, stringsAsFactors = FALSE),
    coords,
    by = "id",
    all.x = TRUE,
    sort = FALSE
  )
  
  # Fallback: if any x/y are NA (shouldn't normally happen),
  # put them on a small secondary circle.
  na_idx <- is.na(nodes$x)
  if (any(na_idx)) {
    m <- sum(na_idx)
    angles2 <- seq(0, 2*pi - 2*pi/m, length.out = m)
    nodes$x[na_idx] <- (radius * 0.6) * sin(angles2)
    nodes$y[na_idx] <- (radius * 0.6) * cos(angles2)
  }
  
  nodes
}



## ======================================================
## 2) LLM-based DAG proposal
## ======================================================

propose_dag_from_llm <- function(vars,
                                 expert_text = NULL,
                                 graph_name = "LLM DAG",
                                 model = Sys.getenv("OPENAI_MODEL", unset = "gpt-4o-mini"),
                                 temperature = 0.2) {
  if (! .have_key()) {
    message("No OPENAI_API_KEY set – falling back to empty DAG.")
    A <- matrix(0L, nrow = length(vars), ncol = length(vars),
                dimnames = list(vars, vars))
    return(A)
  }
  
  sys_msg <- paste(
    "You are an expert in causal discovery and cell signaling.",
    "You will be given a list of variable names from a biological dataset,",
    "optionally with expert background knowledge.",
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
  
  expert_block <- if (!is.null(expert_text) && nzchar(expert_text)) {
    paste0("\n\nExpert background knowledge (free text):\n", expert_text)
  } else {
    ""
  }
  
  usr_msg <- paste0(
    "Variable names:\n",
    paste(vars, collapse = ", "),
    expert_block,
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
  
  # check acyclicity; if cyclic, warn but still return
  g <- igraph::graph_from_adjacency_matrix(A, mode = "directed")
  if (!igraph::is_dag(g)) {
    warning("LLM-proposed adjacency is not a DAG (contains cycles). ",
            "You may need to adjust the prompt or post-process the edges.")
  }
  
  message("LLM DAG proposal: ", sum(A), " directed edges created.")
  A
}

## ======================================================
## 3) DAG–CI agent (your existing logic)
## ======================================================

# IMPORTANT: this assumes you have a function run_ci_tests(amat, dat, tests, alpha)
# defined somewhere (e.g. in another file or package).

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
  
  ## ---- LLM explanation prompt: for clinicians + scientists ----
  sys_msg <- paste(
    "You are a statistician and medical data scientist.",
    "You explain results of conditional independence (CI) tests used to evaluate causal assumptions in biological or clinical datasets.",
    "Your audience is medical professionals and researchers with basic statistical knowledge but not experts in causal inference.",
    "",
    "WRITE CLEARLY AND CONCISELY:",
    "- Use 5–7 short bullet points.",
    "- Avoid formulas and heavy notation.",
    "- Prefer plain language, explain any technical term briefly.",
    "",
    "DEFINITIONS:",
    "- A DAG (Directed Acyclic Graph) represents assumptions about cause-and-effect between variables (e.g. biomarkers, clinical factors).",
    "- A CI test checks whether the data is compatible with the assumption that two variables are independent after adjusting for other variables.",
    "- 'Falsified' means that evidence from the data contradicts the assumed independence (Holm-adjusted p < alpha).",
    "",
    "FOCUS ON:",
    "1) Whether the DAG is broadly compatible with the data or clearly contradicted.",
    "2) How many independence assumptions were tested and how many were rejected.",
    "3) Which rejected assumptions are most important and what they might mean biologically or clinically.",
    "4) A reminder that 'not falsified' does not prove that the DAG or causal story is correct.",
    "5) Optional: suggest simple next steps (e.g. consider adding a link between specific nodes, revisit biological mechanism, or collect more data)."
  )
  
  usr_msg <- paste0(
    "We tested whether the DAG '", graph_name, "' is consistent with this dataset using CI tests.\n",
    sprintf("Significance level alpha = %.3f.\n\n", alpha),
    "Summary of CI tests per test type (number tested, number falsified, smallest adjusted p-value):\n",
    smry_txt, "\n\n",
    "List of falsified independence assumptions (Holm-adjusted p-values):\n",
    rej_txt, "\n\n",
    "Please provide a short explanation in bullet points for clinicians and researchers.\n",
    "Address:\n",
    "- Does the proposed causal structure overall look acceptable or clearly inconsistent with the data?\n",
    "- Roughly how many assumptions were checked and how many did not hold?\n",
    "- Which 1–3 violated assumptions are most important and how might they change the biological or clinical interpretation?\n",
    "- Remind that a non-falsified DAG is not necessarily the true causal graph, only that it is not contradicted by these tests.\n"
  )
  
  
  interp <- openai_chat_min(sys_msg, usr_msg)
  
  # Return everything
  list(
    raw_results    = res,
    summary        = summary_tbl,
    rejected       = rejected_tbl,
    interpretation = interp
  )
}

## ======================================================
## 4) Helper: parse CI strings → (X, Y)
## ======================================================

parse_ci_pair <- function(ci_str) {
  if (is.null(ci_str) || is.na(ci_str) || !nzchar(ci_str)) {
    return(c(NA_character_, NA_character_))
  }
  
  # 1) Normalise independence symbols first
  #    so that we don't cut inside "_||_"
  main <- trimws(ci_str)
  main <- gsub("_\\|\\|_", "⟂", main)   # "A _||_ B"
  main <- gsub("\\|\\|",   "⟂", main)   # "A || B"
  main <- gsub("⫫",       "⟂", main)   # unicode independence
  main <- gsub("⟂⟂",      "⟂", main)   # just in case
  
  # 2) Now remove the conditioning set:
  #    everything starting with a " |"
  #    e.g. "A ⟂ B | Z,W" -> "A ⟂ B"
  main <- sub("\\s*\\|.*$", "", main)
  main <- trimws(main)
  
  # 3) Split on the unified separator
  parts <- strsplit(main, "⟂", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  
  if (length(parts) < 2) {
    return(c(NA_character_, NA_character_))
  }
  
  c(parts[1], parts[2])
}




## ======================================================
## 5) Build matrices:
##    - testable edges (ICs tested) → red dashed
##    - missing edges (rejected ICs) → red solid
## ======================================================

build_testable_matrix <- function(ci_results, vars) {
  if (is.null(ci_results) || !nrow(ci_results)) {
    return(matrix(0L, nrow = length(vars), ncol = length(vars),
                  dimnames = list(vars, vars)))
  }
  
  A <- matrix(0L, nrow = length(vars), ncol = length(vars),
              dimnames = list(vars, vars))
  
  for (ci in ci_results$CI) {
    pair <- parse_ci_pair(ci)
    X <- pair[1]; Y <- pair[2]
    if (!is.na(X) && !is.na(Y) && X %in% vars && Y %in% vars) {
      A[X, Y] <- 1L
      A[Y, X] <- 1L   # symmetrisch, wir markieren nur "es wurde getestet"
    }
  }
  A
}

infer_missing_edges_from_ci <- function(ci_results, vars) {
  if (is.null(ci_results) || !nrow(ci_results)) {
    return(matrix(0L, nrow = length(vars), ncol = length(vars),
                  dimnames = list(vars, vars)))
  }
  
  rejs <- ci_results[ci_results$rejected, , drop = FALSE]
  if (!nrow(rejs)) {
    return(matrix(0L, nrow = length(vars), ncol = length(vars),
                  dimnames = list(vars, vars)))
  }
  
  A_new <- matrix(0L, nrow = length(vars), ncol = length(vars),
                  dimnames = list(vars, vars))
  
  for (ci in rejs$CI) {
    pair <- parse_ci_pair(ci)
    X <- pair[1]; Y <- pair[2]
    if (!is.na(X) && !is.na(Y) && X %in% vars && Y %in% vars) {
      A_new[X, Y] <- 1L
      A_new[Y, X] <- 1L  # skeleton style; orientation not specified
    }
  }
  A_new
}

## ======================================================
## 6) Plot DAG with annotations
## ======================================================

plot_dag_with_annotations <- function(amat_base,
                                      amat_testable = NULL,
                                      amat_new = NULL,
                                      main = "DAG",
                                      layout = NULL) {
  vars <- rownames(amat_base)
  if (is.null(vars)) {
    stop("amat_base must have dimnames (variable names).")
  }
  
  # base directed graph
  g <- igraph::graph_from_adjacency_matrix(amat_base, mode = "directed")
  
  # Default aesthetics for DAG edges
  E(g)$color <- "black"
  E(g)$lty   <- 1
  E(g)$lwd   <- 1
  
  ## 1) Add new “missing” edges in solid red (from rejected CIs)
  if (!is.null(amat_new)) {
    for (i in seq_len(nrow(amat_new))) {
      for (j in seq_len(ncol(amat_new))) {
        if (amat_new[i, j] == 1L &&
            amat_base[i, j] == 0L && amat_base[j, i] == 0L) {
          v1 <- rownames(amat_new)[i]
          v2 <- colnames(amat_new)[j]
          if (v1 %in% vars && v2 %in% vars) {
            g <- igraph::add_edges(
              g, c(v1, v2),
              color = "red", lty = 1, lwd = 3
            )
          }
        }
      }
    }
  }
  
  ## 2) Layout (shared with Shiny via rv_layout)
  if (!is.null(layout)) {
    lay <- layout[igraph::V(g)$name, , drop = FALSE]
  } else {
    lay <- igraph::layout_in_circle(g)
  }
  
  # we control scaling ourselves so we can overlay segments()
  xlim <- range(lay[, 1]) * 1.2
  ylim <- range(lay[, 2]) * 1.2
  
  plot(
    g,
    main            = main,
    vertex.label.cex = 1.1,
    vertex.size      = 25,
    edge.arrow.size  = 0.4,
    layout           = lay,
    rescale          = FALSE,
    xlim             = xlim,
    ylim             = ylim
  )
  
  ## 3) Overlay CI pairs as red dotted lines (undirected)
  if (!is.null(amat_testable)) {
    ci_vars <- rownames(amat_testable)
    if (is.null(ci_vars)) ci_vars <- colnames(amat_testable)
    
    for (i in seq_len(nrow(amat_testable))) {
      for (j in seq_len(ncol(amat_testable))) {
        if (i < j && amat_testable[i, j] == 1L) {
          v1 <- rownames(amat_testable)[i]
          v2 <- colnames(amat_testable)[j]
          
          if (v1 %in% vars && v2 %in% vars) {
            idx1 <- match(v1, vars)
            idx2 <- match(v2, vars)
            # draw an undirected dotted red line between the two nodes
            segments(
              x0 = lay[idx1, 1],
              y0 = lay[idx1, 2],
              x1 = lay[idx2, 1],
              y1 = lay[idx2, 2],
              col = "red",
              lty = 3,     # dotted
              lwd = 2
            )
          }
        }
      }
    }
  }
}




## ======================================================
## 7) Shiny app
## ======================================================

ui <- fluidPage(
  titlePanel("DAG–CI Agent"),
  tags$head(
    tags$style(HTML("
    #llm_interpretation {
      white-space: pre-wrap;   /* allow wrapped lines */
      overflow: visible !important;
      height: auto !important;
      max-height: none !important;
      font-size: 16px;
      line-height: 1.4;
    }
  "))
  ),
  add_busy_spinner(
    spin = "fading-circle",   # type of spinner
    position = "top-right",   # or "top-left", "top-center"
    margins = c(10, 10),
    timeout = 100             # ms delay before it shows
  ),
  tags$head(
    tags$style(HTML("
    #llm_interpretation {
      white-space: pre-wrap;
      overflow: visible !important;
      height: auto !important;
      max-height: none !important;
      font-size: 16px;
      line-height: 1.4;
      border: 1px solid #ccc;
      padding: 12px;
      background-color: #fafafa;
      border-radius: 6px;
    }
  "))
  ),
  
  sidebarLayout(
    sidebarPanel(
      h4("1. Daten hochladen"),
      fileInput("datafile", "Daten hochladen",
                accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls",
                           "text/csv", "application/vnd.ms-excel",
                           "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")),
      hr(),
      uiOutput("var_info"),
      
      hr(),
      h4("2. DAG-Hypothese wählen"),
      radioButtons("dag_source", "DAG-Quelle",
                   choices = c("LLM-Vorschlag"          = "llm",
                               "Adjazenzmatrix hochladen" = "upload",
                               "Manuell zeichnen"        = "manual"),
                   inline = FALSE),
      
      conditionalPanel(
        condition = "input.dag_source == 'llm'",
        textAreaInput("expert_text", "Expert Knowledge (optional)",
                      placeholder = "z.B. 'Akt beeinflusst Erk und PIP3 wirkt auf Akt'",
                      rows = 4),
        actionButton("btn_llm_dag", "LLM-DAG vorschlagen")
      ),
      
      conditionalPanel(
        condition = "input.dag_source == 'upload'",
        fileInput("amat_file", "Adjazenzmatrix-CSV (Zeilen/Spalten = Variablen)"),
        actionButton("btn_upload_dag", "DAG laden")
      ),
      
      hr(),
      h4("3. CI-Tests"),
      numericInput("alpha", "Signifikanzniveau α",
                   value = 0.05, min = 0, max = 0.5, step = 0.01),
      checkboxGroupInput("tests", "CI-Tests",
                         choices = c("GCM" = "gcm", "PCM" = "pcm"),
                         selected = c("gcm", "pcm")),
      actionButton("btn_run_ci", "CI-Tests starten")
    ),
    
    
    
    mainPanel(
      
      conditionalPanel(
        condition = "input.dag_source == 'manual'",
        h4("Manueller DAG-Editor"),
        p("Ziehe Knoten, füge Kanten über das '+'-Symbol oben links hinzu,",
          "oder lösche Kanten über das 'Mülleimer'-Symbol."),
        visNetworkOutput("dag_editor", height = "400px"),
        br(),
        actionButton("btn_manual_apply", "Kanten aus Editor übernehmen und als DAG verwenden")
      ),
      hr(),
      
      
      h4("DAG-Plot"),
      plotOutput("dag_plot", height = "500px"),
      
      hr(),
      h4("Konsistenz-Check"),
      verbatimTextOutput("consistency_msg"),
      
      hr(),
      h4("Zusammenfassung CI-Tests"),
      tableOutput("tbl_summary"),
      
     
      
    
      
      hr(),
      
      
      h4("Abgelehnte ICs (Top 10)"),
      tableOutput("tbl_rejected"),
      
      hr(),
      h4("DAG mit verworfenen CIs"),
      plotOutput("dag_ci_plot", height = "500px"),   
      
      hr(),
      h4("LLM-Erklärung"),
      verbatimTextOutput("llm_interpretation")
      
      
    )
  )
)


adjacency_to_edges <- function(A) {
  vars_row <- rownames(A)
  vars_col <- colnames(A)
  if (is.null(vars_row)) vars_row <- vars_col
  if (is.null(vars_col)) vars_col <- vars_row
  
  idx <- which(A != 0, arr.ind = TRUE)
  if (nrow(idx) == 0) {
    return(data.frame(
      from = character(0),
      to   = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  data.frame(
    from = vars_row[idx[, "row"]],
    to   = vars_col[idx[, "col"]],
    stringsAsFactors = FALSE
  )
}





server <- function(input, output, session) {
  
  # Manual DAG edges from visNetwork
  rv_manual_edges <- reactiveVal(NULL)
  rv_layout       <- reactiveVal(NULL)   
  
  # Adjacency matrix used for CI tests / plotting
  rv_amat <- reactiveVal(NULL)
  
  # CI results
  rv_ci       <- reactiveVal(NULL)
  rv_testable <- reactiveVal(NULL)
  rv_missing  <- reactiveVal(NULL)
  
  ## ---- DATA with auto header/sep detection ----
  dat <- reactive({
    req(input$datafile)
    ext  <- tools::file_ext(input$datafile$name)
    path <- input$datafile$datapath
    
    # 1) Read data into df
    if (ext %in% c("xlsx", "xls")) {
      df <- readxl::read_excel(path)
      df <- as.data.frame(df)
      colnames(df) <- make.names(colnames(df), unique = TRUE)
    } else {
      # auto-detect separator from first line
      line1 <- readLines(path, n = 1)
      sep_candidates <- c("," , ";" , "\t")
      counts <- vapply(
        sep_candidates,
        function(s) {
          m <- gregexpr(s, line1, fixed = TRUE)[[1]]
          if (identical(m, -1L)) 0L else length(m)
        },
        integer(1)
      )
      if (all(counts == 0)) {
        sep <- ","
      } else {
        sep <- sep_candidates[which.max(counts)]
      }
      
      # peek at first few rows to guess header
      peek <- read.table(path,
                         sep = sep,
                         nrows = 5,
                         header = FALSE,
                         stringsAsFactors = FALSE,
                         check.names = FALSE)
      
      row1 <- peek[1, , drop = TRUE]
      row2 <- if (nrow(peek) >= 2) peek[2, , drop = TRUE] else NULL
      
      is_letter <- function(x) grepl("[A-Za-z]", x)
      is_number <- function(x) suppressWarnings(!is.na(as.numeric(x)))
      
      has_letters_row1 <- any(sapply(row1, is_letter))
      all_unique_row1  <- length(unique(as.character(row1))) == ncol(peek)
      has_numbers_row2 <- !is.null(row2) && any(sapply(row2, is_number))
      
      header_guess <- has_letters_row1 && all_unique_row1 && has_numbers_row2
      
      df <- read.csv(path,
                     header = header_guess,
                     sep    = sep,
                     check.names = TRUE)
    }
    
    # 2) Compute shared circular layout for ALL file types
    vars <- colnames(df)
    n <- length(vars)
    
    angles <- seq(from = pi/2,           # 12 o'clock
                  by   = -2*pi/n,       # clockwise
                  length.out = n)
    coords <- cbind(x = cos(angles), y = sin(angles))
    rownames(coords) <- vars
    rv_layout(coords)
    
    df
  })
  
  
  ## ---- VAR INFO ----
  output$var_info <- renderUI({
    req(dat())
    tagList(
      h5("Variablen im Datensatz (verwendete Namen):"),
      verbatimTextOutput("var_names")
    )
  })
  
  output$var_names <- renderText({
    req(dat())
    paste(colnames(dat()), collapse = ", ")
  })
  
  output$dag_editor <- renderVisNetwork({
    req(dat(), rv_layout())
    vars   <- colnames(dat())
    coords <- rv_layout()
    
    # scale coords so the circle is big enough in the canvas
    coords_scaled <- coords * 300  # radius ~300
    
    nodes <- data.frame(
      id      = vars,
      label   = vars,
      x       = coords_scaled[vars, "x"],
      y       = coords_scaled[vars, "y"],
      fixed   = TRUE,     # <- do not let visNetwork move them
      physics = FALSE,    # <- no physics for these nodes
      stringsAsFactors = FALSE
    )
    
    edges <- rv_manual_edges()
    if (is.null(edges) || nrow(edges) == 0) {
      edges <- data.frame(
        from = character(0),
        to   = character(0),
        stringsAsFactors = FALSE
      )
    } else {
      edges <- as.data.frame(edges)
      edges <- edges[, c("from", "to"), drop = FALSE]
    }
    
    visNetwork(nodes, edges, height = "500px") %>%
      visEdges(arrows = "to") %>%
      visOptions(manipulation = TRUE,
                 nodesIdSelection = TRUE)
  })
  
  
  
  # Keep track of edges the user draws / deletes in the editor
  observe({
    edges <- input$dag_editor_edges
    if (!is.null(edges)) {
      rv_manual_edges(as.data.frame(edges))
    }
  })
  
  # When user clicks "Kanten aus Editor übernehmen..."
  observeEvent(input$btn_manual_apply, {
    req(dat())
    vars  <- colnames(dat())
    edges <- rv_manual_edges()
    
    A <- matrix(0L, nrow = length(vars), ncol = length(vars),
                dimnames = list(vars, vars))
    
    if (!is.null(edges) && nrow(edges) > 0) {
      for (i in seq_len(nrow(edges))) {
        from <- as.character(edges$from[i])
        to   <- as.character(edges$to[i])
        if (from %in% vars && to %in% vars && from != to) {
          A[from, to] <- 1L
        }
      }
    }
    
    rv_amat(A)
    showNotification("Manuelle DAG-Kanten übernommen.", type = "message")
  })
  
  ## ---- LLM DAG proposal (and load into editor) ----
  observeEvent(input$btn_llm_dag, {
    req(dat())
    vars <- colnames(dat())
    
    A <- propose_dag_from_llm(vars,
                              expert_text = input$expert_text,
                              graph_name = "LLM DAG")
    
    # store DAG for CI tests
    rv_amat(A)
    
    # edges for manual editor
    edges_llm <- adjacency_to_edges(A)
    rv_manual_edges(edges_llm)
    
    # switch UI to manual editor so user can edit
    updateRadioButtons(session, "dag_source", selected = "manual")
    
    showNotification("LLM-DAG erstellt und in den Editor geladen.", type = "message")
  })
  
  ## ---- Upload adjacency matrix ----
  observeEvent(input$btn_upload_dag, {
    req(input$amat_file, dat())
    amat_raw <- as.matrix(read.csv(input$amat_file$datapath,
                                   row.names = 1,
                                   check.names = FALSE))
    
    if (nrow(amat_raw) != ncol(amat_raw)) {
      showNotification("Adjazenzmatrix ist nicht quadratisch.", type = "error")
      return()
    }
    
    data_vars <- colnames(dat())
    A_vars    <- rownames(amat_raw)
    missing_in_data <- setdiff(A_vars, data_vars)
    if (length(missing_in_data) > 0) {
      showNotification(
        paste("Diese Knoten sind nicht im Datensatz vorhanden:",
              paste(missing_in_data, collapse = ", ")),
        type = "error"
      )
      return()
    }
    
    rv_amat(amat_raw)
    
    # also load into manual editor
    rv_manual_edges(adjacency_to_edges(amat_raw))
    updateRadioButtons(session, "dag_source", selected = "manual")
    
    showNotification("Adjazenzmatrix geladen und in den Editor übernommen.", type = "message")
  })
  
  ## ---- CI-tests ----
  observeEvent(input$btn_run_ci, {
    req(dat(), rv_amat())
    A <- rv_amat()
    
    data_vars <- colnames(dat())
    A_vars    <- rownames(A)
    if (is.null(A_vars)) A_vars <- colnames(A)
    
    missing_in_data <- setdiff(A_vars, data_vars)
    if (length(missing_in_data) > 0) {
      showNotification(
        paste("Fehler: folgende Knoten sind nicht als Daten-Variablen vorhanden:",
              paste(missing_in_data, collapse = ", ")),
        type = "error"
      )
      return()
    }
    
    out <- dag_ci_agent(dat(),
                        amat_llm = A,
                        alpha = input$alpha,
                        tests = input$tests,
                        graph_name = "Hypothesized DAG")
    rv_ci(out)
    
    vars <- rownames(A)
    if (is.null(vars)) vars <- colnames(A)
    
    A_testable <- build_testable_matrix(out$raw_results, vars)
    rv_testable(A_testable)
    
    A_missing <- infer_missing_edges_from_ci(out$raw_results, vars)
    rv_missing(A_missing)
  })
  
  ## ---- Consistency message ----
  output$consistency_msg <- renderPrint({
    ci <- rv_ci()
    req(ci)
    n_rej <- nrow(ci$rejected)
    if (n_rej == 0) {
      cat("Ergebnis: DAG ist NICHT falsifiziert für α =",
          input$alpha, "\n",
          "(keine IC-Statements wurden verworfen).\n")
    } else {
      cat("Ergebnis: DAG ist NICHT konsistent mit den Daten.\n",
          n_rej, "IC-Statements wurden bei α =", input$alpha, "verworfen.\n")
    }
  })
  
  ## ---- Tables ----
  output$tbl_summary <- renderTable({
    ci <- rv_ci()
    req(ci)
    ci$summary
  })
  
  output$tbl_rejected <- renderTable({
    ci <- rv_ci()
    req(ci)
    head(ci$rejected, 10)
  })
  
  output$llm_interpretation <- renderText({
    ci <- rv_ci()
    req(ci)
    if (is.null(ci$interpretation)) "" else ci$interpretation
  })
  
  ## ---- DAG-Plot: nur Hypothesen-DAG ----
  output$dag_plot <- renderPlot({
    req(dat(), rv_layout())
    coords <- rv_layout()
    
    # Basis-DAG: falls noch keiner gesetzt ist, leerer Graph mit allen Variablen
    A_base <- rv_amat()
    if (is.null(A_base)) {
      vars <- colnames(dat())
      A_base <- matrix(0L, nrow = length(vars), ncol = length(vars),
                       dimnames = list(vars, vars))
    }
    
    vars <- rownames(A_base)
    if (is.null(vars)) vars <- colnames(A_base)
    
    layout_mat <- coords[vars, , drop = FALSE]
    
    # keine CI-Information markieren – nur der DAG selbst
    plot_dag_with_annotations(
      amat_base     = A_base,
      amat_testable = NULL,
      amat_new      = NULL,
      main          = "Hypothesized DAG",
      layout        = layout_mat
    )
  })
  
  
  
  
  
  
  ## ---- DAG-Plot mit verworfenen CIs (rot gestrichelt) ----
  output$dag_ci_plot <- renderPlot({
    req(dat(), rv_layout(), rv_ci())
    
    A_base <- rv_amat()
    if (is.null(A_base)) return(invisible(NULL))
    
    coords <- rv_layout()
    
    # DAG-Knotennamen
    vars <- rownames(A_base)
    if (is.null(vars)) vars <- colnames(A_base)
    layout_mat <- coords[vars, , drop = FALSE]
    
    # Hilfsfunktion: Namen kanonisieren (gleicher Raum für DAG & CI-Strings)
    canon <- function(x) {
      x <- tolower(x)
      x <- sub("_(ll|ul)$", "", x)     # z.B. PIP2_ll -> pip2
      x <- gsub("[^a-z0-9]", "", x)    # entferne Sonderzeichen, Punkte, Slash etc.
      x
    }
    
    canon_vars  <- canon(vars)
    name_map    <- setNames(vars, canon_vars)  # canon -> original Name
    
    # 1) Basis-DAG zeichnen
    g <- igraph::graph_from_adjacency_matrix(A_base, mode = "directed")
    
    xlim <- range(layout_mat[, 1]) * 1.2
    ylim <- range(layout_mat[, 2]) * 1.2
    
    plot(
      g,
      main             = "DAG mit verworfenen CIs (rot gestrichelt)",
      vertex.label.cex = 1.1,
      vertex.size      = 25,
      edge.arrow.size  = 0.4,
      layout           = layout_mat,
      rescale          = FALSE,
      xlim             = xlim,
      ylim             = ylim
    )
    
    # 2) Verworfene CIs als rote gestrichelte Linien überlagern
    out  <- rv_ci()
    rejs <- out$rejected
    if (is.null(rejs) || nrow(rejs) == 0) return(invisible(NULL))
    
    # optional: debug in Konsole sehen
    message("Anzahl verworfene CIs: ", nrow(rejs))
    
    for (k in seq_len(nrow(rejs))) {
      ci_str <- as.character(rejs$CI[k])
      pair   <- parse_ci_pair(ci_str)
      X_raw  <- pair[1]
      Y_raw  <- pair[2]
      
      if (is.na(X_raw) || is.na(Y_raw)) next
      
      # auf kanonische Form bringen
      Xc <- canon(X_raw)
      Yc <- canon(Y_raw)
      
      if (!(Xc %in% names(name_map) && Yc %in% names(name_map))) {
        next
      }
      
      X <- name_map[[Xc]]
      Y <- name_map[[Yc]]
      
      i <- match(X, vars)
      j <- match(Y, vars)
      if (is.na(i) || is.na(j)) next
      
      segments(
        x0 = layout_mat[i, 1],
        y0 = layout_mat[i, 2],
        x1 = layout_mat[j, 1],
        y1 = layout_mat[j, 2],
        col = "red",
        lty = 3,   # dotted
        lwd = 3
      )
    }
  })
  
  
  
  
}



shinyApp(ui, server)
