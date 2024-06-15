biomarkers <- function(main_factor, confounder = "") {

    phylo <- readRDS("phyloseq.rds")

    if (confounder == "") {
        f <- main_factor
    } else {
        f <- main_factor + confounder
    }

    phylo@sam_data[[main_factor]] <- factor(phylo@sam_data[[main_factor]], levels = c("Uninfected", "STH", "P. vivax", "Both"))

    out <- ancombc(phyloseq = phylo, formula = f,
                   p_adj_method = "BH", prv_cut = 0.10, lib_cut = 1000,
                   group = main_factor, struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
                   max_iter = 100, conserve = TRUE, alpha = 0.1, global = TRUE)

    # ALL files in res
    res <- out$res
    res_global <- out$res_global

    ## ANCOMBC primary result
    # Coefficients
    tab_coef <- res$lfc
    write.csv(tab_coef, file = "Coeff.csv", row.names = FALSE)

    # SEs
    tab_se <- res$se
    write.csv(tab_se, file = "SEs.csv", row.names = FALSE)

    # Test statistics
    tab_w <- res$W
    write.csv(tab_w, file = "t.statistics.csv", row.names = FALSE)

    # P-values
    tab_p <- res$p_val
    write.csv(tab_p, file = "p-values.csv", row.names = FALSE)

    # Adjusted p-values
    tab_q <- res$q_val
    write.csv(tab_q, file = "q-values.csv", row.names = FALSE)

    # Differentially abundant taxa
    tab_diff <- res$diff_abn
    write.csv(tab_diff, file = "Differentially abundant taxa.csv", row.names = FALSE)

    # Bias-adjusted abundances
    samp_frac <- out$samp_frac
    # Replace NA with 0
    samp_frac[is.na(samp_frac)] <- 0

    # Add pseudo-count (1) to avoid taking the log of 0
    log_obs_abn <- log(microbiome::abundances(phylo) + 1)
    # Adjust the log observed abundances
    log_obs_abn_adj <- t(t(log_obs_abn) - samp_frac)
    write.csv(round(log_obs_abn_adj, 2), file = "Log observed abundances.csv", row.names = TRUE)


    # =====================================
    # LEFSE
    # =====================================


    # Get the taxonomy table from your phyloseq object
    tax_table_phy <- tax_table(phylo)
    # Ensure that taxonomic ranks are correctly named
    colnames(tax_table_phy) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
    # Create a new phyloseq object with the updated taxonomy table
    phy_updated <- phyloseq(otu_table(phy), sample_data(phy), tax_table_phy)
    # Replace NA values in taxonomy table with placeholders
    tax_table_phy[is.na(tax_table_phy)] <- "Unknown"
    # Run LEfSe analysis
    lef_out <- run_lefse(phy_updated, group = main_factor, norm = "CPM",
                         kw_cutoff = 0.05, lda_cutoff = 0.01)#0.01,2)
    ## Plot enrichment analysis plot ##
    p <- plot_ef_bar(lef_out)

    ggsave("LDA.png", plot = p, width = 10, height = 8)

    unique_groups <- length(unique(sample_data(phy_updated)[[main_factor]]))
    color_palette <- rainbow(unique_groups)


    ## Plot cladogram ## make sure that the colour is the same as group number , and save the image with high dimensions
    p <- plot_cladogram(lef_out, color = color_palette, clade_label_level =4 , only_marker = TRUE )

    suppressWarnings(ggsave("Cladogram.png", plot = p, width = 20, height = 8))



}

# Example usage:
# biomarkers(main_factor = "", confounder = "")
