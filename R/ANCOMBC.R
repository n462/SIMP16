perform_ancombc <- function(phylo, main_factor, confounder = "") {
    if (confounder == "") {
        f <- main_factor
    } else {
        f <- paste(main_factor, confounder, sep = "+")
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
}

# Example usage:
# perform_ancombc(phylo = your_phyloseq_object, main_factor = "Category", confounder = "FinalComment")
