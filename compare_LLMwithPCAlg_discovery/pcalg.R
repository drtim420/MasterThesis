# --- PC on Kook-style Sachs data (pure pcalg) -------------------------------
# - Mirrors his read_data(): picks one of the cd3cd28*.xls files, log-transforms,
#   and renames columns to Raf, Mek, PLCg, PIP2, PIP3, Erk, Akt, PKA, PKC, p38, JNK
# - Runs PC (gaussCItest) and saves CPDAG + adjacency + plot
# ---------------------------------------------------------------------------

set.seed(426)

# ------------ Config (edit these two) ------------
int  <- "none"   # one of: "none","Akt","PIP2","Erk","PKC","PIP3"
data_path <- "~/Desktop/master_thesis/code/MasterThesis/data"  # folder containing the xls files
alpha <- 0.01           # PC significance level (try 0.01–0.05)

# ------------ Dependencies ------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# Bioconductor deps required by pcalg
BiocManager::install(c("graph", "RBGL", "Rgraphviz"), ask = FALSE, update = FALSE)
# CRAN deps
install.packages(setdiff(c("pcalg","Matrix","readxl","dplyr"), rownames(installed.packages())),
                 repos = "https://cloud.r-project.org")

library(readxl)
library(dplyr)
library(pcalg)
library(graph)
quietly_load_Rgraphviz <- requireNamespace("Rgraphviz", quietly = TRUE)

# ------------ Mirror Kook's loader ------------
nms <- c("Raf","Mek","PLCg","PIP2","PIP3","Erk","Akt","PKA","PKC","p38","JNK")

read_data_kook <- function(int = c("none","Akt","PIP2","Erk","PKC","PIP3"),
                           path = "../data") {
  int <- match.arg(int)
  file <- switch(int,
                 "none" = "cd3cd28.xls",
                 "Akt"  = "cd3cd28+aktinhib.xls",
                 "PIP2" = "cd3cd28+psitect.xls",
                 "Erk"  = "cd3cd28+u0126.xls",
                 "PKC"  = "cd3cd28+g0076.xls",
                 "PIP3" = "cd3cd28+ly.xls"
  )
  fp <- file.path(path, file)
  if (!file.exists(fp)) stop("Data file not found: ", fp)
  
  dat <- read_xls(fp) |> mutate(across(everything(), log))
  if (ncol(dat) != length(nms)) {
    warning("Column count != 11; keeping numeric columns and renaming.")
    dat <- dat |> select(where(is.numeric))
  }
  if (ncol(dat) != length(nms)) {
    stop("Expected 11 numeric columns after loading; got ", ncol(dat), ".")
  }
  colnames(dat) <- nms
  as.data.frame(dat)
}

X <- read_data_kook(int = int, path = data_path)
cat(sprintf("Loaded %s (n=%d, p=%d)\n", int, nrow(X), ncol(X)))

# ------------ PC algorithm (Gaussian CI test) ------------
suffStat <- list(C = cor(X), n = nrow(X))
pc_fit <- pc(suffStat    = suffStat,
             indepTest   = gaussCItest,
             alpha       = alpha,
             labels      = colnames(X),
             skel.method = "stable",
             u2pd        = "relaxed",
             maj.rule    = TRUE,
             solve.confl = TRUE,
             verbose     = TRUE)

# ------------ Save results ------------
saveRDS(pc_fit, file = sprintf("pc_%s_cpdag.rds", int))
amat <- as(pc_fit@graph, "matrix")
write.csv(amat, file = sprintf("pc_%s_adj.csv", int), row.names = TRUE)

if (quietly_load_Rgraphviz) {
  png(sprintf("pc_%s_graph.png", int), width = 1400, height = 900, res = 150)
  plot(pc_fit, main = sprintf("PC on %s (alpha=%.3f)", int, alpha))
  dev.off()
}

cat("Done. Files written:\n",
    sprintf("  - pc_%s_cpdag.rds\n", int),
    sprintf("  - pc_%s_adj.csv\n", int),
    if (quietly_load_Rgraphviz) sprintf("  - pc_%s_graph.png\n", int) else "")


# now we view it


library(Rgraphviz)
fit <- pc_fit
plot(fit, main = "PC (CD3/CD28, alpha=0.01)")

