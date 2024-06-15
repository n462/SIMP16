#' Title
#'
#' @param fastq_files
#' @param output_dir
#' @param paired
#' @param silva
#' @param env
#' @param trunc
#'
#' @return
#' @export
#'
#' @examples
asv_seq <- function(fastq_files, output_dir, paired, silva,env, trunc) {
    setwd(output_dir)

    # Sort fastq files
    path <- file.path(fastq_files)
    fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
    if (paired)
        fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))

    # Extract sample names
    sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
    write.csv(sample.names, file = "Samples_names.csv", row.names = FALSE)

    # DADA2
    plotF <- plotQualityProfile(fnFs)
    plotR <- plotQualityProfile(fnRs)
    dataF <- plotF$data
    dataR <- plotR$data

    # Write quality stats
    df <- dataF %>%
        group_by(Cycle) %>%
        top_n(3, Count) %>%
        summarise(
            Score = median(Score),
            Count = median(Count)
        )
    write.table(t(df), file = "qual_statsF.tsv", sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)

    df <- dataR %>%
        group_by(Cycle) %>%
        top_n(3, Count) %>%
        summarise(
            Score = median(Score),
            Count = median(Count)
        )
    write.table(t(df), file = "qual_statsR.tsv", sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)

    # Filter and trim
    filtFs <- file.path(output_dir, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
    if (paired)
        filtRs <- file.path(output_dir, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

    names(filtFs) <- sample.names

    if (paired) {
        out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=trunc,
                             maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                             compress=TRUE, multithread=FALSE)
    } else {
        out <- filterAndTrim(fnFs, filtFs,truncLen=trunc,
                             maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
                             compress=TRUE, multithread=env)
    }

    # Sample Inference
    errF <- learnErrors(filtFs, multithread=FALSE)
    errR <- learnErrors(filtRs, multithread=FALSE)
    dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
    if (paired)
        dadaRs <- dada(filtRs, err=errR, multithread=env)

    # Merging
    if (paired)
        mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

    # Construct ASV table
    if (paired) {
        seqtab <- makeSequenceTable(mergers)
    } else {
        seqtab <- makeSequenceTable(dadaFs)
    }
    write.csv(seqtab, file = "seqtab.csv", row.names = FALSE)

    # Remove chimeras
    seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
    write.csv(seqtab.nochim, file = "asv_data.csv", row.names = FALSE)


    # Assign taxonomy
    taxa <- assignTaxonomy(seqtab.nochim, silva, multithread=env)
    write.csv(taxa, file = "taxonomy.csv")
    saveRDS(taxa, file = "ASVs_taxonomy.rds")

    # Return relevant data
    return(list(seqtab = seqtab, seqtab_nochim = seqtab.nochim, taxonomy = taxa))
}

# # Test the function
# fastq_files <- "G:/tets_DATA"
# output_dir <-  "G:/tets_DATA"
#result <- asv_seq(fastq_files, output_dir, paired=TRUE, silva="",env=, trunc = NULL) # env = Windows, Ubuntu FALSE, TRUE
