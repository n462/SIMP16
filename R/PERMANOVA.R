perform_permanova <- function(phylo, main_factor, confounder = NULL) {
    # clr transformation
    pseq.rel <- microbiome::transform(phylo, "clr")

    # Get abundances and metadata
    otu <- abundances(pseq.rel)
    meta <- meta(pseq.rel)

    # Method: euclidean, jaccard, ...
    # Perm = 999, 1000, ...

    # clr distances
    clr_dist_matrix <- phyloseq::distance(pseq.rel, method = "euclidean")
    # clr permanova
    h1 <- with(meta, how(nperm = 999,
                         within = Within(type = "free")))

    formula <- reformulate(main_factor, response = "clr_dist_matrix")
    if (!is.null(confounder)) {
        formulaConf <- reformulate(paste(main_factor, confounder, sep = "+"), response = "clr_dist_matrix")
        permconf <- pairwise.adonis2(formulaConf, data = meta)
        # Remove the parent_call element
        permconf$parent_call <- NULL

        # Remove the class attribute
        attr(permconf, "class") <- NULL

        file_path_conf <- "PERMANOVA_Confounder.txt"
        con_conf <- file(file_path_conf, open = "wt")

        # Write each comparison block to the file
        for (comparison in names(permconf)) {
            writeLines(comparison, con_conf)
            write.table(permconf[[comparison]], con_conf, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
            writeLines("\n", con_conf)
        }

        close(con_conf)
    }

    perm <- pairwise.adonis2(formula, data = meta)

    # Remove the parent_call element
    perm$parent_call <- NULL

    # Remove the class attribute
    attr(perm, "class") <- NULL

    file_path <- "PERMANOVA.txt"
    con <- file(file_path, open = "wt")

    # Write each comparison block to the file
    for (comparison in names(perm)) {
        writeLines(comparison, con)
        write.table(perm[[comparison]], con, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
        writeLines("\n", con)
    }

    close(con)
}

# Example usage:
# perform_permanova(phylo = your_phyloseq_object, main_factor = "Category", confounder = "FinalComment")
