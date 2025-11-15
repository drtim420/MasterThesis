### Analysis of the observational distribution
### Falsifying consensus and Sachs graph
### LK 2025

set.seed(1)
graphs <- c("Consensus", "Sachs") # Which graph
# Observational or intervened graph
interventions <- c("none", "Akt", "PIP2", "Erk", "PKC", "PIP3")
stylize <- FALSE
stylized <- c("Akt", "Erk") # c("JNK", "p38")
col <- ifelse(stylize, "cornflowerblue", "orange")
save <- TRUE

### Dependencies
library("comets")
library("tidyverse")
library("readxl")
library("dagitty")
library("ggdag")
source("functions.R")

nop <- lapply(graphs, \(wg) {
  lapply(interventions, \(int) {
    ### Read data
    dat <- read_data(int = int)
    amat <- get_dag(tolower(wg), int = int)
    dag <- adj2dag(amat)

    # Number of sides
    n <- 11

    # Radius of the polygon
    r <- 1 # You can change this value if needed

    # Compute coordinates
    theta <- seq(0, 2 * pi, length.out = n + 1)[-(n + 1)]
    # Exclude duplicate last point
    x <- r * cos(theta)
    y <- r * sin(theta)

    # Combine into a data frame
    names(x) <- nms
    names(y) <- nms
    coords <- list(x = x, y = y)

    coordinates(dag) <- coords

    if (stylize) {
      tidy_dag <- tidy_dagitty(dag) |>
        mutate(color = ifelse(name %in% stylized, "int", "obs"))
    } else {
      tidy_dag <- tidy_dagitty(dag) |>
        mutate(color = ifelse(name == int, "int", "obs"))
    }

    ggplot(tidy_dag, aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_dag_edges() +
      geom_dag_point() +
      geom_dag_text(aes(color = color), show.legend = FALSE) +
      theme_dag() +
      scale_color_manual(values = c("int" = col, "obs" = "white")) +
      labs(subtitle = paste(wg, "graph"))

    if (save) {
      if (!dir.exists("../figures")) {
        dir.create("../figures")
      }
      ggsave(
        file.path("../figures", paste0(
          tolower(wg), "-", int,
          ifelse(stylize, paste0("-stylized-", paste0(stylized, collapse = "-")), ""),
          ".pdf"
        )),
        height = 4, width = 4
      )
    }
  })
})
