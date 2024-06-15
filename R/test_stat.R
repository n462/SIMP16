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
#' @param main_factor
#' @param confounder
#' @param method
#' @param perm
#'
#' @return
#' @export
#'
#' @examples
statdiversity<- function(main_factor, confounder, method, perm) {
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
    output_dir = getwd()
    # =====================================
    # RCM Analysis
    # =====================================
    phylo <- readRDS("phyloseq.rds")
    asv <- RCM(phylo)
    plot <- plot(asv, samColour = main_factor)
    ggsave(file.path(output_dir, "RCM.png"), plot = plot, width = 10, height = 8)

    asv_confounder <- RCM(phylo, confounders = confounder, prevCutOff = 0.2)
    plot <- plot(asv_confounder, samColour = main_factor)
    ggsave(file.path(output_dir, "RCM_confounder.png"), plot = plot, width = 10, height = 8)

    # =====================================
    # Alpha Diversity
    # =====================================

    # Chao1 Alpha Diversity
    alpha_data <- estimate_richness(phylo, measures = "Chao1")
    alpha_data_df <- as.data.frame(alpha_data)
    alpha_data_df[[main_factor]] <- as.factor(sample_data(phylo)[[main_factor]])
    p <- ggplot(alpha_data_df, aes(x = .data[[main_factor]], y = Chao1, fill = .data[[main_factor]])) +
        geom_boxplot(outlier.shape = NA) +
        ylab("Chao1 Alpha Diversity") +
        xlab(main_factor) +
        theme_classic()
    ggsave(file.path(output_dir, "Chao1.png"), plot = p, width = 10, height = 8)

    # InvSimpson Alpha Diversity
    alpha_data <- estimate_richness(phylo, measures = "InvSimpson")
    alpha_data_df <- as.data.frame(alpha_data)
    alpha_data_df[[main_factor]] <- as.factor(sample_data(phylo)[[main_factor]])
    p <- ggplot(alpha_data_df, aes(x = .data[[main_factor]], y = InvSimpson, fill = .data[[main_factor]])) +
        geom_boxplot(outlier.shape = NA) +
        ylab("InvSimpson Alpha Diversity") +
        xlab(main_factor) +
        theme_classic()
    ggsave(file.path(output_dir, "InvSimpson.png"), plot = p, width = 10, height = 8)

    # Shannon Alpha Diversity
    alpha_data <- estimate_richness(phylo, measures = "Shannon")
    alpha_data_df <- as.data.frame(alpha_data)
    alpha_data_df[[main_factor]] <- as.factor(sample_data(phylo)[[main_factor]])
    p <- ggplot(alpha_data_df, aes(x = .data[[main_factor]], y = Shannon, fill = .data[[main_factor]])) +
        geom_boxplot(outlier.shape = NA) +
        ylab("Shannon Alpha Diversity") +
        xlab(main_factor) +
        theme_classic()
    ggsave(file.path(output_dir, "Shannon.png"), plot = p, width = 10, height = 8)

    # =====================================
    # PERMANOVA
    # =====================================
    pseq.rel <- microbiome::transform(phylo, "clr")
    otu <- abundances(pseq.rel)
    meta <- meta(pseq.rel)

    clr_dist_matrix <- phyloseq::distance(pseq.rel, method = method)
    h1 <- with(meta, how(nperm = perm, within = Within(type = "free")))

    formula <- reformulate(main_factor, response = "clr_dist_matrix")
    formulaConf <- reformulate(paste(main_factor, confounder, sep="+"), response = "clr_dist_matrix")

    perm <- pairwise.adonis2(formula, data = meta) #, permutations = h1
    perm$parent_call <- NULL
    attr(perm, "class") <- NULL
    file_path <- file.path(output_dir, "PERMANOVA.txt")
    con <- file(file_path, open = "wt")
    for (comparison in names(perm)) {
        writeLines(comparison, con)
        write.table(perm[[comparison]], con, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
        writeLines("\n", con)
    }
    close(con)

    permconf <- pairwise.adonis2(formulaConf, data = meta)
    permconf$parent_call <- NULL
    attr(permconf, "class") <- NULL
    file_path <- file.path(output_dir, "PERMANOVA_Confounder.txt")
    con <- file(file_path, open = "wt")
    for (comparison in names(permconf)) {
        writeLines(comparison, con)
        write.table(permconf[[comparison]], con, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
        writeLines("\n", con)
    }
    close(con)
}



# Test the statistical_analysis function
#statdiversity(main_factor = "", confounder = "", method=c('euclidean','jaccard'), perm=c(999,1000))
