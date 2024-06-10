convert_to_phyloseq <- function(sharedfile, constaxfile, mapfile) {
    # Import data into mothur format
    mothur_data <- import_mothur(
        mothur_shared_file = sharedfile,
        mothur_constaxonomy_file = constaxfile
    )

    # Read map file
    map <- read.csv

    # transform dataframe into phyloseq class object.
    map <- sample_data(map)

    # Assign rownames to be Sample ID's
    # sample id########
    rownames(map) <- map$host_subject_id

    # Merge mothur data with sample metadata
    mothur_merge <- merge_phyloseq(mothur_data, map)

    return(mothur_merge)
}

# Example usage:
# phylo <- convert_to_phyloseq("ASVshared.shared", "contaxonomy.taxonomy", "mapfile.csv")
