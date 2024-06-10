perform_dada2 <- function(
        fnFs, fnRs, output_dir, sample.names, paired = TRUE, errF = NULL, errR = NULL
) {
    # Plot quality profiles
    plotF <- plotQualityProfile(fnFs)
    plotR <- plotQualityProfile(fnRs)

    dataF <- plotF$data
    dataR <- plotR$data

    df <- dataF %>%
        group_by(Cycle) %>%
        top_n(3, Count) %>%
        summarise(
            Score = median(Score),
            Count = median(Count)
        )

    # Write the data frame to the CSV file
    write.table(t(df), file = "qual_statsF.tsv", sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)

    df <- dataR %>%
        group_by(Cycle) %>%
        top_n(3, Count) %>%
        summarise(
            Score = median(Score),
            Count = median(Count)
        )

    write.table(t(df), file = "qual_statsR.tsv", sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)

    library(reticulate)
    reticulate::repl_python()

    # Create new directory for filtered files
    dir.create(file.path(output_dir, "filtered"), showWarnings = FALSE, recursive = TRUE)

    filtFs <- file.path(output_dir, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
    if (paired == TRUE)
        filtRs <- file.path(output_dir, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
    names(filtFs) <- sample.names

    if (paired == TRUE){
        out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(251, 250),
                             maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE,
                             compress = TRUE, multithread = TRUE)
    } else {
        out <- filterAndTrim(fnFs, filtFs, truncLen = 130,
                             maxN = 0, maxEE = 2, truncQ = 2, rm.phix = TRUE,
                             compress = TRUE, multithread = TRUE)
    }

    dadaFs <- dada(filtFs, err = errF, multithread = FALSE)
    if (paired == TRUE)
        dadaRs <- dada(filtRs, err = errR, multithread = FALSE)

    if (paired == TRUE)
        mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

    if (paired == TRUE){
        seqtab <- makeSequenceTable(mergers)
    } else {
        seqtab <- makeSequenceTable(dadaFs)
    }

    file_path <- file.path(output_dir, "seqtab.csv")
    write.csv(seqtab, file = file_path, row.names = FALSE)
    dim(seqtab)

    table(nchar(getSequences(seqtab)))

    seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE, verbose = TRUE)

    file_path <- file.path(output_dir, "asv_data.csv")
    write.csv(seqtab.nochim, file = file_path, row.names = FALSE)
    dim(seqtab.nochim)

    sum(seqtab.nochim) / sum(seqtab)

    if (paired == TRUE){
        getN <- function(x) sum(getUniques(x))
        track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
        colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
        rownames(track) <- sample.names
        file_path <- file.path(output_dir, "Summary.csv")
        write.csv(track, file = file_path, row.names = FALSE)
    }

    taxa <- assignTaxonomy(seqtab.nochim, "/media/aya/One Touch/Study 2/silva_nr_v132_train_set.fa.gz", multithread = TRUE)
    taxa.print <- taxa
    rownames(taxa.print) <- NULL

    file_path <- file.path(output_dir, "taxonomy.csv")
    write.csv(taxa, file = file_path)

    file_path <- file.path(output_dir, "ASVs_taxonomy.rds")
    saveRDS(taxa, file_path)
}
