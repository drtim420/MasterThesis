# run_kook_on_data_folder.R
# Runs Kook/COMETs on ALL datasets found in ../data (csv/xlsx/xls).
# Data/DAG generation stays in data_simulation_with_DAG.R.

set.seed(1)
alpha <- 0.05
save  <- TRUE
data_path <- "../data"
tests <- c("gcm","pcm")

library(comets)
library(tidyverse)
library(readxl)
library(dagitty)
library(tidyr)

# bring node list + helpers (NO generation here)
source("data_simulation_with_DAG.R")   # provides nms, get_dag(), adj2dag(), etc.
use_saved_dag <- FALSE   # <- FALSE = test ALL datasets against get_dag() (the original)

# ---- discover all datasets (csv/xlsx/xls) ----
files_csv  <- list.files(data_path, pattern = "\\.csv$",  full.names = TRUE, ignore.case = TRUE)
files_xlsx <- list.files(data_path, pattern = "\\.xlsx$", full.names = TRUE, ignore.case = TRUE)
files_xls  <- list.files(data_path, pattern = "\\.xls$",  full.names = TRUE, ignore.case = TRUE)
files <- unique(c(files_csv, files_xlsx, files_xls))
if (length(files) == 0) stop("No datasets found in ", data_path)

# ---- loader: reads data + finds matching DAG (basename_dagitty.txt) ----
load_dataset_with_dag <- function(file_path) {
  ext <- tools::file_ext(file_path) |> tolower()
  dat <- switch(ext,
                "csv"  = readr::read_csv(file_path, show_col_types = FALSE),
                "xlsx" = readxl::read_xlsx(file_path),
                "xls"  = readxl::read_xls(file_path),
                stop("Unsupported extension: ", ext)
  )
  
  # Keep only variables we know; reorder to nms if possible
  if (!all(nms %in% colnames(dat))) {
    warning("Skipping (missing required columns): ", basename(file_path))
    return(NULL)
  }
  dat <- dat[, nms, drop = FALSE]
  colnames(dat) <- nms
  
  base <- tools::file_path_sans_ext(basename(file_path))
  dir  <- dirname(file_path)
  dag_file <- file.path(dir, paste0(base, "_dagitty.txt"))
  
  # choose hypothesis DAG
  if (use_saved_dag) {
    if (file.exists(dag_file)) {
      dag <- dagitty::dagitty(paste(readLines(dag_file), collapse = " "))
      amat <- matrix(0, length(nms), length(nms), dimnames = list(nms, nms))
      eds  <- dagitty::edges(dag)
      if (nrow(eds)) for (k in seq_len(nrow(eds))) amat[eds[k,1], eds[k,2]] <- 1
    } else {
      amat <- get_dag()  # fallback to original if no saved DAG
    }
  } else {
    amat <- get_dag()    # <-- ORIGINAL DAG as hypothesis for ALL datasets
  }
  list(data = dat, amat = amat, dataset = basename(file_path))
  
}

# ---- run COMETs on one dataset ----
run_ci_tests <- function(amat, dat, tests = c("gcm","pcm")) {
  dag  <- adj2dag(amat)
  cis  <- impliedConditionalIndependencies(dag)
  tcis <- cis[vapply(cis, function(x) length(x$Z) > 0, logical(1))]
  
  run_one <- function(tst) {
    pv <- vapply(tcis, function(ci) {
      fm <- reformulate(
        paste0(paste0(ci$Y, collapse = "+"), "|", paste0(ci$Z, collapse = "+")),
        response = ci$X
      )
      comets(fm, dat, test = tst, coin = TRUE)$p.value
    }, numeric(1))
    tibble(test = tst, CI = sapply(tcis, paste), p.value = pv,
           adj.p.value = p.adjust(pv, "holm"))
  }
  bind_rows(lapply(tests, run_one))
}

# ---- batch runner over ALL files in ../data ----
batched <- list()
for (f in files) {
  obj <- try(load_dataset_with_dag(f), silent = TRUE)
  if (inherits(obj, "try-error") || is.null(obj)) next
  res <- run_ci_tests(obj$amat, obj$data, tests = tests) %>%
    mutate(dataset = obj$dataset)
  batched[[length(batched)+1]] <- res
}
if (length(batched) == 0) stop("No valid datasets with required columns found.")

res <- bind_rows(batched)

# ---- summaries & reporting ----
summary_tbl <- res %>%
  mutate(viol = adj.p.value < alpha) %>%
  group_by(dataset, test) %>%
  summarise(n_CI = n(), n_viol = sum(viol), prop_viol = mean(viol), .groups = "drop")

cat("Datasets tested:", length(unique(res$dataset)), "\n")
print(summary_tbl %>% arrange(dataset, test))

viol <- res %>% filter(adj.p.value < alpha) %>% arrange(adj.p.value)
if (nrow(viol) == 0) {
  message("No CI violations at Holm-adjusted alpha = ", alpha, ".")
} else {
  print(viol, n = Inf)
}

# optional: top-10 smallest raw p per test
res %>% arrange(p.value) %>% group_by(test) %>% slice_head(n = 10) %>% ungroup() %>% print(n = Inf)

# ---- save outputs ----
if (save) {
  dir.create("../results", showWarnings = FALSE)
  readr::write_csv(res, "../results/all-tests_multi.csv")
  readr::write_csv(summary_tbl, "../results/summary_multi.csv")
  if (nrow(viol) > 0) readr::write_csv(viol, "../results/violations_multi.csv")
}
