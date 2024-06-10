perform_RCM_analysis <- function(phylo, main_factor, confounder = NULL, output_dir = "./") {
    # Convert main_factor to a factor if necessary
    if (!is.factor(phylo@sam_data[[main_factor]])) {
        phylo@sam_data[[main_factor]] <- factor(phylo@sam_data[[main_factor]])
    }

    # Choose the order of the target columns if needed

    # Run RCM
    if (is.null(confounder)) {
        asv <- RCM(phylo)
    } else {
        asv <- RCM(phylo, confounders = confounder, prevCutOff = 0.2)
    }

    # Save RCM object
    saveRDS(asv, file.path(output_dir, "RCM_object.rds"))

    # Plot results
    plot_main <- plot(asv, samColour = main_factor)
    plot_confounder <- plot(asv_confounder, samColour = main_factor)

    # Save plots
    ggsave(file.path(output_dir, "RCM_plot.png"), plot = plot_main, width = 10, height = 8)
    ggsave(file.path(output_dir, "RCM_plot_confounder.png"), plot = plot_confounder, width = 10, height = 8)
}

# Example usage:
# perform_RCM_analysis(phylo, main_factor = "comparison_group", confounder = NULL, output_dir = "/path/to/save/directory/")
