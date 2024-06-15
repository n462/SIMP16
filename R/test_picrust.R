#' Title
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param destination
#' @param asvtax
#' @param asvpi
#'
#' @return
#' @export
#'
#' @examples
picrust <- function(destination, asvtax, asvpi) {
    # Ensure necessary packages are installed
    install_dependencies()

    # Load the necessary libraries
    library(ape)
    library(dplyr)
    library(ggplot2)
    library(plotly)
    library(tidyr)
    library(vegan)
    library(ggpubr)
    library(gplots)
    library(lme4)
    library(gt)
    library(DT)
    library(tidyverse)
    library(knitr)
    library(usethis)
    library(gmp)
    library(nloptr)
    library(Rmpfr)
    library(lmerTest)
    library(Cairo)
    library(phangorn)
    library(phyloseq)
    library(mia)
    library(RCM)
    library(ANCOMBC)
    library(ShortRead)
    library(dada2)
    library(microbiome)
    library(microbiomeMarker)
    library(DECIPHER)
    library(remotes)
    library(mothur)
    # Download the PICRUSt reference files
    download_ref(destination = destination, reference = 'gg_ko')

    # Convert the ASV and taxonomy files to a suitable format
    df <- asvtax
    df <- df[, -1]
    new_header <- df[1, ]
    df <- df[-1, ]
    names(df) <- new_header

    dfasv <- asvpi
    names(dfasv) <- sub("^ASV", "", names(dfasv))
    col_names <- names(dfasv)
    names(dfasv)[1] <- ""
    rownames(dfasv) <- dfasv[, 1]
    dfasv <- dfasv[, -1]

    # Perform PICRUSt analysis
    t <- t4f(otu_table = dfasv, rows_are_taxa = FALSE, tax_table = df,
             reference_path = destination, type = 'uproc', short = TRUE,
             cn_normalize = TRUE, sample_normalize = TRUE, drop = TRUE)

    # Extract KEGG descriptions
    fxn_table_225_final <- t$fxn_meta
    kegg <- fxn_table_225_final$KEGG_Description

    # Write the KEGG descriptions to a file
    write.table(kegg, 'Picrust.txt', sep = '\t')
}

# Example usage:
#destination <- "/media/aya/One Touch/GRAD1/references"
#asvtax <- read.csv("taxtest.csv")  # Assuming you have a file named "taxonomy.csv" containing taxonomy data
#asvpi <- read.csv("ASVtest.csv")    # Assuming you have a file named "asv_data.csv" containing ASV data

# Call the function
#picrust(destination, asvtax, asvpi)
