process_taxonomy_and_metadata <- function(taxonomy_file, metadata_file, output_taxonomy_file, output_shared_file, mapfile) {
    # Read the taxonomy data
    new_data <- read.csv(taxonomy_file)

    # Rename the first column to 'ASV'
    names(new_data)[1] <- "ASV"

    # Rename the rows of the first column
    new_data$ASV <- paste0("ASV", 1:nrow(new_data))
    asvtax <- new_data

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
    write.table(cbind(contaxonomy[, 1:2], contax), file = output_taxonomy_file, sep = "\t", row.names = FALSE, quote = FALSE)

    taxoncom <- cbind(contaxonomy[, 1:2], contax)

    # Read metadata file
    metadata <- read.csv(metadata_file)

    # Save metadata file
    write.csv(metadata, file = output_shared_file)
}

# Example usage:
# process_taxonomy_and_metadata("taxonomy.csv", "/media/aya/One Touch/Study 5/Metadata.csv", "contaxonomy.taxonomy", "ASVshared.shared", "mapfile.csv")
