convert_to_shared_format <- function(asv_file, samples_file, output_file) {
    # Read the ASV data and sample names
    asv <- read.csv(asv_file)
    group_data <- read.csv(samples_file, header = TRUE)

    # Rename columns
    new_column_names <- paste0("ASV", seq_along(asv))
    colnames(asv) <- new_column_names

    # Create additional columns for shared format
    label_column <- rep(0.03, nrow(asv))
    numAsvs_column <- rep(ncol(asv), nrow(asv))
    asvshared <- cbind(label = label_column, Group = group_data[, 1], numAsvs = numAsvs_column, asv)

    # Write data frame to CSV with tab-separated values
    write.table(asvshared, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
}

# Example usage:
# convert_to_shared_format("asv_data.csv", "Samples_names.csv", "ASVshared.shared")
