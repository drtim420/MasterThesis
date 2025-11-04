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
  ), ncol = 11)

  ### SACHS graph
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
  ), ncol = 11)

  amat <- switch(which,
    "consensus" = consensus,
    "sachs" = sachs
  )
  colnames(amat) <- nms
  rownames(amat) <- nms

  ### Remove incoming nodes (if intervened)
  if (int != "none") {
    amat[, int] <- 0
  }

  amat
}

read_data <- function(int = c("none", "Akt", "PIP2", "Erk", "PKC", "PIP3"), path = "../data") {
  file <- switch(int,
    "none" = "cd3cd28.xls",
    "Akt" = "cd3cd28+aktinhib.xls",
    "PIP2" = "cd3cd28+psitect.xls",
    "Erk" = "cd3cd28+u0126.xls",
    "PKC" = "cd3cd28+g0076.xls",
    "PIP3" = "cd3cd28+ly.xls"
  )
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

tcis2txt <- function(x) {
  sapply(x, paste0)
}
