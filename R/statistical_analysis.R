statistical_analysis <- function(phylo, main_factor, confounder = "", output_dir) {
    library(phyloseq)
    library(microbiome)
    library(ggplot2)
    library(dplyr)
    library(vegan)
    library(RCM)

    setwd(output_dir)

    # =====================================
    # RCM Analysis
    # =====================================
    asv <- RCM(phylo)
    plot <- plot(asv, samColour = main_factor)
    ggsave(file.path(output_dir, "RCM.png"), plot = plot, width = 10, height = 8)

    asv_confounder <- RCM(phylo, confounders = confounder, prevCutOff = 0.2)
    plot <- plot(asv_confounder, samColour = main_factor)
    ggsave(file.path(output_dir, "RCM_confounder.png"), plot = plot, width = 10, height = 8)

    # =====================================
    # Alpha Diversity
    # =====================================
    metadata <- sample_data(phylo)
    category <- metadata$Category
    otu_table <- as(otu_table(phylo), "matrix")

    kruskal_results <- apply(otu_table, 1, function(otu_row) {
        kruskal.test(otu_row ~ category)
    })

    p_values <- sapply(kruskal_results, function(x) x$p.value)
    adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
    significant_otus <- which(adjusted_p_values < 0.05)
    significant_results <- data.frame(OTU = names(significant_otus), Adjusted_p_value = adjusted_p_values[significant_otus])
    print(significant_results)

    # Chao1 Alpha Diversity
    alpha_data <- estimate_richness(phylo, measures = "Chao1")
    alpha_data_df <- as.data.frame(alpha_data)
    alpha_data_df[[main_factor]] <- as.factor(sample_data(phylo)[[main_factor]])
    p <- ggplot(alpha_data_df, aes(x = .data[[main_factor]], y = Chao1, fill = .data[[main_factor]])) +
        geom_boxplot(outlier.shape = NA) +
        ylab("Chao1 Alpha Diversity") +
        xlab(main_factor) +
        theme_classic()
    ggsave(file.path(output_dir, "Chao1.png"), plot = p, width = 10, height = 8)

    # InvSimpson Alpha Diversity
    alpha_data <- estimate_richness(phylo, measures = "InvSimpson")
    alpha_data_df <- as.data.frame(alpha_data)
    alpha_data_df[[main_factor]] <- as.factor(sample_data(phylo)[[main_factor]])
    p <- ggplot(alpha_data_df, aes(x = .data[[main_factor]], y = InvSimpson, fill = .data[[main_factor]])) +
        geom_boxplot(outlier.shape = NA) +
        ylab("InvSimpson Alpha Diversity") +
        xlab(main_factor) +
        theme_classic()
    ggsave(file.path(output_dir, "InvSimpson.png"), plot = p, width = 10, height = 8)

    # Shannon Alpha Diversity
    alpha_data <- estimate_richness(phylo, measures = "Shannon")
    alpha_data_df <- as.data.frame(alpha_data)
    alpha_data_df[[main_factor]] <- as.factor(sample_data(phylo)[[main_factor]])
    p <- ggplot(alpha_data_df, aes(x = .data[[main_factor]], y = Shannon, fill = .data[[main_factor]])) +
        geom_boxplot(outlier.shape = NA) +
        ylab("Shannon Alpha Diversity") +
        xlab(main_factor) +
        theme_classic()
    ggsave(file.path(output_dir, "Shannon.png"), plot = p, width = 10, height = 8)

    # =====================================
    # PERMANOVA
    # =====================================
    pseq.rel <- microbiome::transform(phylo, "clr")
    otu <- abundances(pseq.rel)
    meta <- meta(pseq.rel)

    clr_dist_matrix <- phyloseq::distance(pseq.rel, method = "euclidean")
    h1 <- with(meta, how(nperm = 999, within = Within(type = "free")))

    formula <- reformulate(main_factor, response = "clr_dist_matrix")
    formulaConf <- reformulate(paste(main_factor, confounder, sep="+"), response = "clr_dist_matrix")

    perm <- pairwise.adonis2(formula, data = meta) #, permutations = h1
    perm$parent_call <- NULL
    attr(perm, "class") <- NULL
    file_path <- file.path(output_dir, "PERMANOVA.txt")
    con <- file(file_path, open = "wt")
    for (comparison in names(perm)) {
        writeLines(comparison, con)
        write.table(perm[[comparison]], con, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
        writeLines("\n", con)
    }
    close(con)

    permconf <- pairwise.adonis2(formulaConf, data = meta)
    permconf$parent_call <- NULL
    attr(permconf, "class") <- NULL
    file_path <- file.path(output_dir, "PERMANOVA_Confounder.txt")
    con <- file(file_path, open = "wt")
    for (comparison in names(permconf)) {
        writeLines(comparison, con)
        write.table(permconf[[comparison]], con, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
        writeLines("\n", con)
    }
    close(con)
}
