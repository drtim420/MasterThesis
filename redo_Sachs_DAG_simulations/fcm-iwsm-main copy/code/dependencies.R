pkgs <- c(
  "comets",
  "tidyverse",
  "readxl",
  "dagitty",
  "ggdag"
)

install.packages("remotes", repos = "https://cloud.r-project.org")
remotes::install_cran(pkgs, repos = "https://cloud.r-project.org")
