perform_alpha_diversity_analysis <- function(phylo, main_factor, output_dir) {
    # Assuming 'Category' is the metadata column with four groups
    metadata <- sample_data(phylo)
    category <- metadata$Category

    # Extract the OTU table from the phyloseq object
    otu_table <- as(otu_table(phylo), "matrix")

    # Perform Kruskal-Wallis test for each OTU
    kruskal_results <- apply(otu_table, 1, function(otu_row) {
        kruskal.test(otu_row ~ category)
    })

    # Extract p-values from the Kruskal-Wallis test results
    p_values <- sapply(kruskal_results, function(x) x$p.value)

    # Adjust p-values for multiple testing if necessary
    adjusted_p_values <- p.adjust(p_values, method = "bonferroni")

    # Identify significant OTUs
    significant_otus <- which(adjusted_p_values < 0.05)

    # Display significant OTUs and their adjusted p-values
    significant_results <- data.frame(OTU = names(significant_otus), Adjusted_p_value = adjusted_p_values[significant_otus])
    print(significant_results)

    # Calculate Chao1 alpha diversity
    alpha_data <- estimate_richness(phylo, measures = "Chao1")

    # Convert the data to a format suitable for plotting
    alpha_data_df <- as.data.frame(alpha_data)
    alpha_data_df[[main_factor]] <- as.factor(sample_data(phylo)[[main_factor]])

    # Create boxplot
    p <- ggplot(alpha_data_df, aes(x = .data[[main_factor]], y = Chao1, fill = .data[[main_factor]])) +
        geom_boxplot(outlier.shape = NA) +
        ylab("Chao1 Alpha Diversity") +
        xlab(main_factor) +
        theme_classic()
    ggsave(file.path(output_dir, "Chao1.png"), plot = p, width = 10, height = 8)

    # Calculate InvSimpson alpha diversity
    alpha_data <- estimate_richness(phylo, measures = "InvSimpson")

    # Convert the data to a format suitable for plotting
    alpha_data_df <- as.data.frame(alpha_data)
    alpha_data_df[[main_factor]] <- as.factor(sample_data(phylo)[[main_factor]])

    # Create boxplot
    p <- ggplot(alpha_data_df, aes(x = .data[[main_factor]], y = InvSimpson, fill = .data[[main_factor]])) +
        geom_boxplot(outlier.shape = NA) +
        ylab("InvSimpson Alpha Diversity") +
        xlab(main_factor) +
        theme_classic()
    ggsave(file.path(output_dir, "InvSimpson.png"), plot = p, width = 10, height = 8)

    # Calculate Shannon alpha diversity
    alpha_data <- estimate_richness(phylo, measures = "Shannon")

    # Convert the data to a format suitable for plotting
    alpha_data_df <- as.data.frame(alpha_data)
    alpha_data_df[[main_factor]] <- as.factor(sample_data(phylo)[[main_factor]])

    # Create boxplot
    p <- ggplot(alpha_data_df, aes(x = .data[[main_factor]], y = Shannon, fill = .data[[main_factor]])) +
        geom_boxplot(outlier.shape = NA) +
        ylab("Shannon Alpha Diversity") +
        xlab(main_factor) +
        theme_classic()
    ggsave(file.path(output_dir, "Shannon.png"), plot = p, width = 10, height = 8)
}

# Example usage:
# Specify the phyloseq object, main factor, and output directory
perform_alpha_diversity_analysis(phylo = your_phyloseq_object, main_factor = "Category", output_dir = "./output")
