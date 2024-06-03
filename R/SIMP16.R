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

# Function to process FASTQ files
process_fastq_files <- function(fastq_files, output_dir) {
    path <- file.path(fastq_files)

    fnFs <- sort(list.files(path, pattern = "_1.fastq", full.names = TRUE))
    fnRs <- sort(list.files(path, pattern = "_2.fastq", full.names = TRUE))

    sample_names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
    folder_path <- output_dir
    file_path <- file.path(folder_path, "Samples_names.csv")

    write.csv(sample_names, file = file_path, row.names = FALSE)

    readfiles_length <- countLines(fnRs) / 4
    sum_readfiles_length <- sum(readfiles_length)

    if (sum_readfiles_length > 9e6) {
        stop(paste("There is a maximum number of reads (", sum_readfiles_length, ") to be processed by this pipeline."))
    }
}

# Example usage
fastq_files <- "path_to_fastq_files" # Update with the correct path
output_dir <- "path_to_output_dir" # Update with the correct path

tryCatch({
    process_fastq_files(fastq_files, output_dir)
}, error = function(e) {
    cat("An error occurred: ", e$message, "\n")
})
