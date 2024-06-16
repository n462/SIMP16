# Load reticulate
library(reticulate)

# Function to source the Python script from the package directory
source_trim_script <- function() {
    package_dir <- system.file("python", "Trim.py", package = "SIMP16")
    source_python(package_dir)
}

# Example R function that calls a Python function from Trim.py
run_trim <- function(forward_stats_file, reverse_stats_file, paired) {
    source_trim_script()
    trunc_positions <- py$trunc(forward_stats_file, reverse_stats_file, paired)
    return(trunc_positions)
}

#' ASV Processing and Taxonomy Assignment
#'
#' @description
#' This function processes 16S rRNA sequences from .fastq files, performs ASV inference using DADA2, filters and trims reads,
#' removes chimeras, assigns taxonomy using the SILVA database, and saves the ASV table and taxonomy assignments.
#'
#' @param fastq_files A string specifying the input directory containing the 16S rRNA .fastq files.
#' @param output_dir  A string specifying the output directory where the results will be saved.
#' @param paired      A logical value indicating whether the reads are paired-end (TRUE) or single-end (FALSE).
#' @param silva       A string specifying the directory containing the silva.gz file used for assigning taxonomy.
#' @param env         A logical value indicating the operating system environment: Linux (TRUE), Windows (FALSE).
#' @param trunc       A numeric vector specifying the trimming positions for the forward and reverse reads. If NULL, trimming positions will be automatically detected by default.
#'
#' @return ASV and Taxonomy tables in the current directory
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' asv_seq(fastq_files = "/path/to/.fastq", output_dir = "/path/to/output", paired=TRUE,
#'         silva="/path/to/silva.gz", env=TRUE, trunc = NULL)
#' }
#'
asv_seq <- function(fastq_files, output_dir, paired, silva, env, trunc) {
    setwd(output_dir)

    # Sort fastq files
    path <- file.path(fastq_files)
    fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
    fnRs <- if (paired) sort(list.files(path, pattern="_2.fastq", full.names = TRUE)) else NULL

    # Extract sample names
    sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
    write.csv(sample.names, file = "Samples_names.csv", row.names = FALSE)

    # DADA2 Quality Profile
    plotF <- plotQualityProfile(fnFs)
    plotR <- if (paired) plotQualityProfile(fnRs) else NULL
    dataF <- plotF$data
    dataR <- if (paired) plotR$data else NULL

    if (is.null(trunc)) {
        # Write quality stats
        dfF <- dataF %>%
            group_by(Cycle) %>%
            top_n(3, Count) %>%
            summarise(
                Score = median(Score),
                Count = median(Count)
            )
        write.table(t(dfF), file = "qual_statsF.tsv", sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)

        if (paired) {
            dfR <- dataR %>%
                group_by(Cycle) %>%
                top_n(3, Count) %>%
                summarise(
                    Score = median(Score),
                    Count = median(Count)
                )
            write.table(t(dfR), file = "qual_statsR.tsv", sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
        }

        # Use the Python script to determine truncation positions
        trunc <- run_trim("qual_statsF.tsv", if (paired) "qual_statsR.tsv" else NULL, paired)
        trunc <- as.numeric(unlist(trunc))
    }

    # Filter and trim
    filtFs <- file.path(output_dir, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
    filtRs <- if (paired) file.path(output_dir, "filtered", paste0(sample.names, "_R_filt.fastq.gz")) else NULL

    if (paired) {
        out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = trunc,
                             maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE,
                             compress = TRUE, multithread = env)
    } else {
        out <- filterAndTrim(fnFs, filtFs, truncLen = trunc,
                             maxN = 0, maxEE = 2, truncQ = 2, rm.phix = TRUE,
                             compress = TRUE, multithread = env)
    }

    # Sample Inference
    errF <- learnErrors(filtFs, multithread = env)
    errR <- if (paired) learnErrors(filtRs, multithread = env) else NULL
    dadaFs <- dada(filtFs, err = errF, multithread = env)
    dadaRs <- if (paired) dada(filtRs, err = errR, multithread = env) else NULL

    # Merging
    mergers <- if (paired) mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE) else NULL

    # Construct ASV table
    seqtab <- if (paired) makeSequenceTable(mergers) else makeSequenceTable(dadaFs)
    write.csv(seqtab, file = "seqtab.csv", row.names = FALSE)

    # Remove chimeras
    seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = env, verbose = TRUE)
    write.csv(seqtab.nochim, file = "asv_data.csv", row.names = FALSE)

    # Assign taxonomy
    taxa <- assignTaxonomy(seqtab.nochim, silva, multithread = env)
    write.csv(taxa, file = "taxonomy.csv")
    saveRDS(taxa, file = "ASVs_taxonomy.rds")

    return("ASV and taxonomy tables are successfully saved in the current directory...")
}
