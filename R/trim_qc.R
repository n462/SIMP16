# Required packages
required_packages <- c(
    "ape", "dplyr", "ggplot2", "phangorn", "phyloseq",
    "plotly", "tidyr", "vegan", "RCM", "ANCOMBC",
    "devtools", "ShortRead", "dada2"
)

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

# Function to install and load required packages
install_and_load <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        if (pkg %in% c("phangorn", "phyloseq", "RCM", "ANCOMBC")) {
            BiocManager::install(pkg)
        } else if (pkg == "dada2") {
            devtools::install_github("benjjneb/dada2", ref = "v1.16")
        } else {
            install.packages(pkg)
        }
    }
    library(pkg, character.only = TRUE)
}

# Install and load each package
lapply(required_packages, install_and_load)

# Define the function
trim_qc <- function(fastq_files, output_dir, paired = TRUE) {
    setwd(output_dir)
    path <- file.path(fastq_files)

    fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
    if (paired == TRUE)
        fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))

    sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
    write.csv(sample.names, file = "Samples_names.csv", row.names = FALSE)

}

# Example usage
#fastq_files <- "G:/tets_DATA"
#output_dir <-  "G:/tets_DATA"
#trim_qc(fastq_files, output_dir)
