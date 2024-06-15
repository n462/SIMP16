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
#' @param mapfile
#' @param colID
#'
#' @return
#' @export
#'
#' @examples
get_phylo <- function(mapfile, colID) {
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

    # Read the ASV data and sample names
    asv <- read.csv("asv_data.csv")
    group_data <- read.csv("Samples_names.csv", header = TRUE)

    # Rename columns
    new_column_names <- paste0("ASV", seq_along(asv))
    colnames(asv) <- new_column_names

    # Create additional columns for shared format
    label_column <- rep(0.03, nrow(asv))
    numAsvs_column <- rep(ncol(asv), nrow(asv))
    asvshared <- cbind(label = label_column, Group = group_data[, 1], numAsvs = numAsvs_column, asv)
    asvpi <- cbind(s=group_data[,1], asv)
    write.csv(asvpi, file = "ASVtest.csv")

    # Write data frame to CSV with tab-separated values
    write.table(asvshared, file = "ASVshared.shared" , sep = "\t", row.names = FALSE, quote = FALSE)

    # Read the taxonomy data
    new_data <- read.csv("taxonomy.csv")

    # Rename the first column to 'ASV'
    names(new_data)[1] <- "ASV"

    # Rename the rows of the first column
    new_data$ASV <- paste0("ASV", 1:nrow(new_data))
    asvtax <- new_data
    write.csv(asvtax, file = "taxtest.csv")

    # Sum each row and create a new column with the sums
    data <- colSums(asvshared[4:ncol(asvshared)])

    # Prepare a new data frame with only the sum column
    size_data <- data.frame(size = data)

    # If you want to combine the sum column with the original data and write it to a new file:
    # Prepare a new data frame with the sum column added as the second column
    contaxonomy <- cbind(ASV = new_data[, 1], size = size_data$size, new_data[, 2:ncol(new_data)])

    # Merge columns from the second to the last, using ";" as the delimiter
    merged_column <- apply(contaxonomy[, -c(1, 2)], 1, paste, collapse = ";")

    # Create a new data frame with the merged column under the header "Taxonomy"
    contax <- data.frame(Taxonomy = merged_column)

    # Write the combined data to a new taxonomy CSV file
    write.table(cbind(contaxonomy[, 1:2], contax), file = "contaxonomy.taxonomy", sep = "\t", row.names = FALSE, quote = FALSE)

    taxoncom <- cbind(contaxonomy[, 1:2], contax)

    write.csv(taxoncom, file = "taxoncom.csv")

    # Import data into mothur format
    asvshared = "ASVshared.shared"
    contaxonomy = "contaxonomy.taxonomy"
    mothur_data <- import_mothur(
        mothur_shared_file = asvshared,
        mothur_constaxonomy_file = contaxonomy
    )

    # Read map file
    map <- read.csv(mapfile, sep = ",")

    # Transform dataframe into phyloseq class object
    map <- sample_data(map)

    # Assign rownames to be Sample IDs
    rownames(map) <- map[[colID]]

    # Merge mothur data with sample metadata
    mothur_merge <- merge_phyloseq(mothur_data, map)
    file_path <- file.path(getwd(), "phyloseq.rds")
    saveRDS(mothur_merge, file_path)
}

# Example usage:
#get_phylo("G://tets_DATA//metadata.csv", colID="ID")
