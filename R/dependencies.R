#' Load Required Libraries
#'
#' This function loads all the libraries required for this package.
#'
#' @export
load_required_libraries <- function() {
    libraries <- c(
        "ape", "dplyr", "ggplot2", "plotly", "tidyr", "vegan", "mia",
        "devtools", "ShortRead", "dada2", "ggpubr", "gplots", "lme4",
        "pairwiseAdonis", "microbiome", "gt", "DT", "knitr"
    )

    for (lib in libraries) {
        if (!requireNamespace(lib, quietly = TRUE)) {
            install.packages(lib, dependencies = TRUE)
        }
        library(lib, character.only = TRUE)
    }
}

# Handle Bioconductor packages separately
load_bioconductor_libraries <- function() {
    bioc_libraries <- c("phangorn", "phyloseq", "RCM", "ANCOMBC")
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    for (lib in bioc_libraries) {
        if (!requireNamespace(lib, quietly = TRUE)) {
            BiocManager::install(lib)
        }
        library(lib, character.only = TRUE)
    }
}
