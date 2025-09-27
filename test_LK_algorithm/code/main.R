# testing the LK method on synthetic data and given DAG's

# 1) Build DAG + generate observational synthetic data (CSV)
source("data_simulation_with_DAG.R")     # defines nms, get_dag(), adj2dag(), simulate_from_dag(), generate_observational()
generate_observational(n = 1000, seed = 42)  # writes ../data/cd3cd28.csv

set.seed(1)
save <- TRUE

library(comets)
library(tidyverse)
library(readxl)
library(dagitty)

# 2) Load the synthetic observational dataset once
dat <- readr::read_csv("../data/cd3cd28.csv", show_col_types = FALSE)
colnames(dat) <- nms

# 3) Helper to get the DAG (matrix) — your get_dag() has no args, so ignore graph/int
get_dag_any <- function(graph_name = "New", int = "none") {
  get_dag()
}

# Use only the synthetic graph and only the observational setting
graphs <- c("New")
ints   <- c("none")

# 4) Run the COMET tests (GCM/PCM) on CIs implied by the synthetic DAG
res <- lapply(ints, \(int) {
  cat("Running intervention:", int, "\n")
  do.call("rbind", lapply(graphs, \(wg) {
    cat("Running graph:", wg, "\n")
    do.call("rbind", lapply(c("gcm","pcm"), \(tst) {
      cat("Running test:", tst, "\n")
      
      # DAG + CIs
      amat <- get_dag_any(wg, int)
      dag  <- adj2dag(amat)
      cis  <- impliedConditionalIndependencies(dag)
      tcis <- cis[vapply(cis, function(x) length(x$Z) > 0, logical(1))]
      
      # Test all CIs with COMETs
      pv <- vapply(tcis, \(ci) {
        fm <- reformulate(
          paste0(paste0(ci$Y, collapse = "+"), "|", paste0(ci$Z, collapse = "+")),
          response = ci$X
        )
        comets(fm, dat, test = tst, coin = TRUE)$p.value
      }, numeric(1))
      
      data.frame(
        graph = wg, intervention = int, CI = sapply(tcis, paste),
        p.value = pv, adj.p.value = p.adjust(pv, "holm"), test = tst
      )
    }))
  }))
}) |> do.call("rbind", args = _)




# === After 'res' is created ===
alpha <- 0.05

# 1) List CI violations (Holm-adjusted)
viol <- res %>%
  dplyr::filter(adj.p.value < alpha) %>%
  dplyr::arrange(adj.p.value)

if (nrow(viol) == 0) {
  message("No CI violations at Holm-adjusted alpha = ", alpha, ".")
} else {
  print(viol, n = Inf)
  readr::write_csv(viol, "../results/violations.csv")
}

# 2) Quick counts by test
res %>%
  dplyr::mutate(viol = adj.p.value < alpha) %>%
  dplyr::count(test, viol) %>%
  print()

# 3) (Optional) See which CIs fail across tests (wide view)
library(tidyr)
viol_wide <- res %>%
  dplyr::select(CI, test, adj.p.value) %>%
  tidyr::pivot_wider(names_from = test, values_from = adj.p.value)
viol_wide %>% dplyr::filter(dplyr::if_all(-CI, ~ . < alpha)) %>% print(n = Inf)



# 5) Save results
if (save) {
  dir.create("../results", showWarnings = FALSE)
  readr::write_csv(res, "../results/all-tests.csv")
}

# Optional: quick peek
print(head(res))
