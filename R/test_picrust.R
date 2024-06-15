#' Title
#'
#' @param destination
#' @param asvtax
#' @param asvpi
#'
#' @return
#' @export
#'
#' @examples
picrust <- function(destination, asvtax, asvpi) {
    # Download the PICRUSt reference files
    download_ref(destination = destination, reference = 'gg_ko')

    # Convert the ASV and taxonomy files to a suitable format
    df <- asvtax
    df <- df[, -1]
    new_header <- df[1, ]
    df <- df[-1, ]
    names(df) <- new_header

    dfasv <- asvpi
    names(dfasv) <- sub("^ASV", "", names(dfasv))
    col_names <- names(dfasv)
    names(dfasv)[1] <- ""
    rownames(dfasv) <- dfasv[, 1]
    dfasv <- dfasv[, -1]

    # Perform PICRUSt analysis
    t <- t4f(otu_table = dfasv, rows_are_taxa = FALSE, tax_table = df,
             reference_path = destination, type = 'uproc', short = TRUE,
             cn_normalize = TRUE, sample_normalize = TRUE, drop = TRUE)

    # Extract KEGG descriptions
    fxn_table_225_final <- t$fxn_meta
    kegg <- fxn_table_225_final$KEGG_Description

    # Write the KEGG descriptions to a file
    write.table(kegg, 'Picrust.txt', sep = '\t')
}

# Example usage:
#destination <- "/media/aya/One Touch/GRAD1/references"
#asvtax <- read.csv("taxtest.csv")  # Assuming you have a file named "taxonomy.csv" containing taxonomy data
#asvpi <- read.csv("ASVtest.csv")    # Assuming you have a file named "asv_data.csv" containing ASV data

# Call the function
#picrust(destination, asvtax, asvpi)
