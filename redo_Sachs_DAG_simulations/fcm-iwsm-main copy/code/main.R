### Analysis of the observational distribution
### Falsifying consensus and Sachs graph
### LK 2025

set.seed(1)
save <- TRUE

### Dependencies
library("comets")
library("tidyverse")
library("readxl")
library("dagitty")
source("functions.R")

res <- lapply(
  c("none", "Akt", "PIP2", "Erk", "PKC", "PIP3"),
  \(int) {
    cat("Running intervention:", int, "\n")
    lapply(c("Consensus", "Sachs"), \(wg) {
      cat("Running graph:", wg, "\n")

      lapply(c("gcm", "pcm"), \(tst) {
        cat("Running test:", tst, "\n")

        ### Read data
        dat <- read_data(int = int)
        amat <- get_dag(tolower(wg), int = int)
        dag <- adj2dag(amat)
        cis <- impliedConditionalIndependencies(dag)

        tcis <- cis[unlist(lapply(cis, \(x) length(x$Z) > 0))]
        cat("Found", length(tcis), "CI constraints\n")

        pv <- lapply(tcis, \(ci) {
          fm <- reformulate(
            paste0(
              paste0(ci$Y, collapse = "+"), "|",
              paste0(ci$Z, collapse = "+")
            ),
            response = ci$X
          )
          comets(fm, dat, test = tst, coin = TRUE)$p.value
        }) |> unlist()
        adjp <- p.adjust(pv, "holm")

        data.frame(
          graph = wg,
          intervention = int,
          CI = sapply(tcis, paste),
          p.value = pv,
          adj.p.value = adjp,
          test = tst
        )
      }) |> do.call("rbind", args = _)
    }) |> do.call("rbind", args = _)
  }
) |> do.call("rbind", args = _)

### Only JNK _||_ p38 | PKA, PKC is rejected on the interventional datasets
res |> filter(adj.p.value < 0.05, graph == "Sachs", intervention != "none")

### Filter the CIs for the table
cis <- res |>
  filter(graph == "Consensus", intervention == "none", adj.p.value < 0.05) |>
  pull("CI")

### Make table
res |>
  filter(CI %in% cis) |>
  filter((graph == "Consensus" & intervention == "none") |
    (graph != "Consensus" & intervention != "none")) |>
  select(intervention, graph, test, CI, adj.p.value) |>
  mutate(adj.p.value = format.pval(adj.p.value,
    eps = 0.0001,
    scientific = FALSE, digits = 1
  )) |>
  knitr::kable(
    escape = FALSE, row.names = FALSE,
    booktabs = TRUE, format = "latex", align = "lrrrr"
  )

### Save results
if (save) {
  if (!dir.exists("../results")) {
    dir.create("../results")
  }
  write_csv(res, "../results/all-tests.csv")
}
