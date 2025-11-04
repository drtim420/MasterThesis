## ============================================
## Standalone DAG data diagnostics: linear vs nonlinear (Sachs)
## ============================================

set.seed(123)


data_path    <- "../data"       
int_choice   <- "none"          
graph_choice <- "sachs"        
out_dir      <- "../results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


need_pkg <- function(p) if (!requireNamespace(p, quietly = TRUE))
  stop(sprintf("Please install '%s' first: install.packages('%s')", p, p))
need_pkg("readxl")
need_pkg("mgcv")


nms <- c("Raf","Mek","PLCg","PIP2","PIP3","Erk","Akt","PKA","PKC","p38","JNK")


get_dag <- function(which = c("consensus","sachs"),
                    int = c("none","Akt","PIP2","Erk","PKC","PIP3")) {
  which <- match.arg(which)
  int   <- match.arg(int)
  
  consensus <- matrix(c(
    0,0,0,0,0,0,0,1,1,0,0,
    1,0,0,0,0,0,0,1,1,0,0,
    0,0,0,0,1,0,0,0,0,0,0,
    0,0,1,0,1,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,
    0,1,0,0,0,0,0,1,0,0,0,
    0,0,0,0,1,0,0,1,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,
    0,0,1,1,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,1,1,0,0,
    0,0,0,0,0,0,0,1,1,0,0
  ), ncol = 11, byrow = TRUE)
  
  sachs <- matrix(c(
    0,0,0,0,0,0,0,1,1,0,0,
    1,0,0,0,0,0,0,1,1,0,0,
    0,0,0,0,0,0,0,0,0,0,0,
    0,0,1,0,1,0,0,0,0,0,0,
    0,0,1,0,0,0,0,0,0,0,0,
    0,1,0,0,0,0,0,1,0,0,0,
    0,0,0,0,0,1,0,1,0,0,0,
    0,0,0,0,0,0,0,0,1,0,0,
    0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,1,1,0,0,
    0,0,0,0,0,0,0,1,1,0,0
  ), ncol = 11, byrow = TRUE)
  
  amat <- switch(which, "consensus" = consensus, "sachs" = sachs)
  colnames(amat) <- nms
  rownames(amat) <- nms
  
  
  if (int != "none") amat[, int] <- 0L
  amat
}

## Read one dataset (and log-transform) 
read_sachs_data <- function(int = c("none","Akt","PIP2","Erk","PKC","PIP3"),
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
  full <- file.path(path, file)
  if (!file.exists(full)) {
    stop("File not found: ", full,
         "\nSet 'data_path' correctly or place the Sachs .xls files there.")
  }
  dat <- readxl::read_xls(full)
  # Log-transform all numeric columns (Sachs convention)
  for (j in seq_along(dat)) if (is.numeric(dat[[j]])) dat[[j]] <- log(dat[[j]])
  # Keep/rename columns to nms order if possible
  if (ncol(dat) >= length(nms)) {
    dat <- dat[, seq_len(length(nms))]
  }
  colnames(dat) <- nms
  dat
}

## Parents finder 
get_parents <- function(amat, v) {
  vars <- rownames(amat)
  vars[which(amat[, v] == 1L)]
}


lm_formula  <- function(v, parents) as.formula(paste(v, "~", paste(parents, collapse = " + ")))
gam_formula <- function(v, parents) {
  sm <- paste(sprintf("s(%s)", parents), collapse = " + ")
  as.formula(paste(v, "~", sm))
}

##  5-fold CV RMSE 
cv_rmse <- function(fit_fun, dat, form, K = 5) {
  n <- nrow(dat)
  id <- sample(rep(1:K, length.out = n))
  se <- numeric(K)
  resp <- all.vars(form)[1]
  for (k in 1:K) {
    tr <- dat[id != k, , drop = FALSE]
    te <- dat[id == k, , drop = FALSE]
    fit <- fit_fun(form, tr)
    pr  <- try(predict(fit, newdata = te), silent = TRUE)
    if (inherits(pr, "try-error")) { se[k] <- NA; next }
    se[k] <- mean((te[[resp]] - pr)^2, na.rm = TRUE)
  }
  sqrt(mean(se, na.rm = TRUE))
}

fit_lm  <- function(form, dat) stats::lm(form, data = dat)
fit_gam <- function(form, dat) mgcv::gam(form, data = dat, method = "REML")

## Load data + DAG (with surgery) 
cat("\n--- Loading data & DAG ---\n")
dat  <- read_sachs_data(int = int_choice, path = data_path)
amat <- get_dag(which = graph_choice, int = int_choice)
stopifnot(all(colnames(dat) == rownames(amat)))

cat("Dataset:", int_choice, "| Graph:", graph_choice,
    "| Rows:", nrow(dat), "Cols:", ncol(dat), "\n")

## Quick EDA
cat("\n--- EDA: summary statistics (per variable) ---\n")
summ <- t(vapply(dat, function(x) {
  c(mean = mean(x), sd = sd(x), min = min(x),
    p25 = quantile(x, 0.25), median = median(x),
    p75 = quantile(x, 0.75), max = max(x))
}, numeric(7)))
print(round(summ, 3))

## Correlation heatmap (PNG) 
corr <- stats::cor(dat, use = "pairwise.complete.obs")
png(file.path(out_dir, sprintf("corr_%s_%s.png", int_choice, graph_choice)),
    width = 1100, height = 1000)
par(mar = c(8, 8, 3, 1))
image(1:ncol(corr), 1:ncol(corr), t(corr[nrow(corr):1, ]),
      axes = FALSE, main = sprintf("Correlation heatmap: %s / %s",
                                   toupper(int_choice), toupper(graph_choice)))
axis(1, at = 1:ncol(corr), labels = colnames(corr), las = 2, cex.axis = 0.7)
axis(2, at = 1:ncol(corr), labels = rev(rownames(corr)), las = 2, cex.axis = 0.7)
dev.off()
cat("Saved correlation heatmap to:", file.path(out_dir, sprintf("corr_%s_%s.png", int_choice, graph_choice)), "\n")

## Main: per-node linear vs nonlinear comparison
cat("\n--- Linear vs Nonlinear (LM vs GAM) ---\n")
vars <- colnames(dat)
rows <- list()

for (v in vars) {
  parents <- get_parents(amat, v)
  if (length(parents) == 0) next  # root nodes: skip (no predictors)
  
  f_lm  <- lm_formula(v, parents)
  f_gam <- gam_formula(v, parents)
  
  m_lm  <- try(fit_lm(f_lm, dat),  silent = TRUE)
  m_gam <- try(fit_gam(f_gam, dat), silent = TRUE)
  
  if (inherits(m_lm, "try-error") || inherits(m_gam, "try-error")) {
    rows[[length(rows)+1]] <- data.frame(
      node = v, k_parents = length(parents),
      lm_ok = !inherits(m_lm, "try-error"),
      gam_ok = !inherits(m_gam, "try-error"),
      AIC_lm = NA, AIC_gam = NA,
      RMSE_lm = NA, RMSE_gam = NA,
      better = NA_character_,
      stringsAsFactors = FALSE
    )
    next
  }
  
  AIC_lm  <- AIC(m_lm)
  AIC_gam <- AIC(m_gam)
  RMSE_lm  <- cv_rmse(fit_lm,  dat, f_lm,  K = 5)
  RMSE_gam <- cv_rmse(fit_gam, dat, f_gam, K = 5)
  
  better <- if (is.finite(RMSE_gam) && is.finite(RMSE_lm)) {
    if (RMSE_gam + 1e-8 < RMSE_lm) "GAM (nonlinear)" else "LM (linear)"
  } else if (is.finite(AIC_gam) && is.finite(AIC_lm)) {
    if (AIC_gam + 1e-8 < AIC_lm) "GAM (nonlinear)" else "LM (linear)"
  } else NA_character_
  
  rows[[length(rows)+1]] <- data.frame(
    node = v, k_parents = length(parents),
    lm_ok = TRUE, gam_ok = TRUE,
    AIC_lm = AIC_lm, AIC_gam = AIC_gam,
    RMSE_lm = RMSE_lm, RMSE_gam = RMSE_gam,
    better = better,
    stringsAsFactors = FALSE
  )
}

summary_tbl <- if (length(rows)) do.call(rbind, rows) else
  data.frame(node=character(), k_parents=integer(), lm_ok=logical(),
             gam_ok=logical(), AIC_lm=double(), AIC_gam=double(),
             RMSE_lm=double(), RMSE_gam=double(), better=character())

print(summary_tbl[order(summary_tbl$better, -(summary_tbl$RMSE_lm - summary_tbl$RMSE_gam)), ], row.names = FALSE)

# Save CSV
csv_out <- file.path(out_dir, sprintf("linearity_summary_%s_%s.csv", int_choice, graph_choice))
utils::write.csv(summary_tbl, csv_out, row.names = FALSE)
cat("Saved summary CSV to:", csv_out, "\n")

## plot GAM smooths for top nonlinear nodes
top_nonlinear <- subset(summary_tbl, better == "GAM (nonlinear)" & is.finite(RMSE_gam) & is.finite(RMSE_lm))
if (nrow(top_nonlinear) == 0) {
  cat("\nNo strong nonlinear wins to plot.\n")
} else {
  # rank by gain
  top_nonlinear$gain <- top_nonlinear$RMSE_lm - top_nonlinear$RMSE_gam
  top3 <- head(top_nonlinear[order(-top_nonlinear$gain), ], 3)
  for (i in seq_len(nrow(top3))) {
    v <- top3$node[i]
    parents <- get_parents(amat, v)
    f_gam   <- gam_formula(v, parents)
    m_gam   <- mgcv::gam(f_gam, data = dat, method = "REML")
    
    png(file.path(out_dir, sprintf("gam_smooths_%s_%s_%s.png", v, int_choice, graph_choice)),
        width = 1100, height = 900)
    op <- par(no.readonly = TRUE)
    par(mfrow = c(ceiling(length(parents)/2), 2))
    plot(m_gam, shade = TRUE, pages = 1,
         main = sprintf("%s ~ smooth(parents) — %s / %s", v, toupper(int_choice), toupper(graph_choice)))
    par(op)
    dev.off()
    cat("Saved GAM smooths for", v, "\n")
  }
}

cat("\nDone. Outputs in:", normalizePath(out_dir), "\n")
