perform_lefse <- function(phylo, main_factor) {
    # Get the taxonomy table from your phyloseq object
    tax_table_phy <- tax_table(phylo)
    # Ensure that taxonomic ranks are correctly named
    colnames(tax_table_phy) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
    # Create a new phyloseq object with the updated taxonomy table
    phy_updated <- phyloseq(otu_table(phylo), sample_data(phylo), tax_table_phy)
    # Replace NA values in taxonomy table with placeholders
    tax_table_phy[is.na(tax_table_phy)] <- "Unknown"
    # Run LEfSe analysis
    lef_out <- run_lefse(phy_updated, group = main_factor, norm = "CPM",
                         kw_cutoff = 0.05, lda_cutoff = 0.01)

    ## Plot enrichment analysis plot ##
    p <- plot_ef_bar(lef_out)
    ggsave("LDA.png", plot = p, width = 10, height = 8)

    unique_groups <- length(unique(sample_data(phy_updated)[[main_factor]]))
    color_palette <- rainbow(unique_groups)

    ## Plot cladogram ##
    # Make sure that the colour is the same as group number, and save the image with high dimensions
    p <- plot_cladogram(lef_out, color = color_palette, clade_label_level = 4, only_marker = TRUE )
    suppressWarnings(ggsave("Cladogram.png", plot = p, width = 20, height = 8))

    # Store everything into a list
    results_list <- list(
        lef_out = lef_out,
        enrichment_plot = p,
        color_palette = color_palette
    )

    return(results_list)
}

# Example usage:
# lefse_results <- perform_lefse(phylo = your_phyloseq_object, main_factor = "Category")
