#' Install Dependencies
#'
#' This function installs all necessary CRAN, Bioconductor, and GitHub packages.
#' @import ape
#' @import dplyr
#' @import ggplot2
#' @import plotly
#' @import tidyr
#' @import vegan
#' @import ggpubr
#' @import gplots
#' @import lme4
#' @import gt
#' @import DT
#' @import tidyverse
#' @import knitr
#' @import usethis
#' @import gmp
#' @import nloptr
#' @import Rmpfr
#' @import lmerTest
#' @import Cairo
#' @import phangorn
#' @import phyloseq
#' @import mia
#' @import RCM
#' @import ANCOMBC
#' @import ShortRead
#' @import dada2
#' @import microbiome
#' @import microbiomeMarker
#' @import DECIPHER
#' @importFrom EESI themetagenomics
#' @importFrom pmartinezarbizu pairwiseAdonis
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
