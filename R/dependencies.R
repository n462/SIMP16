#' Install Dependencies
#'
#' This function installs all necessary CRAN, Bioconductor, and GitHub packages.
#' @export
install_dependencies <- function() {
    cran_packages <- c("ape", "dplyr", "ggplot2", "plotly", "tidyr", "vegan",
                       "ggpubr", "gplots", "lme4", "gt", "DT", "tidyverse",
                       "knitr", "usethis", "gmp", "nloptr", "Rmpfr",
                       "lmerTest", "Cairo")

    bioc_packages <- c("phangorn", "phyloseq", "mia", "RCM", "ANCOMBC",
                       "ShortRead", "dada2", "microbiome", "microbiomeMarker",
                       "DECIPHER")

    github_packages <- c("EESI/themetagenomics",
                         "pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

    # Install CRAN packages
    install.packages(cran_packages)

    # Install Bioconductor packages
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install(bioc_packages, update = TRUE, ask = FALSE)

    # Install GitHub packages
    if (!requireNamespace("remotes", quietly = TRUE))
        install.packages("remotes")
    remotes::install_github(github_packages)
}

NULL
