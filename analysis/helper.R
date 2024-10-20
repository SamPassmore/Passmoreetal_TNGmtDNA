# helper script for Supplemantary materials

summarizer = function(posterior){
  posterior$generation_time = as.numeric(gsub("g_", replacement = "", posterior$generation))
  # posterior$clade_probability = posterior$numerator / posterior$numerator
  posterior %>%
    filter(!is.na(clade_probability)) %>% # remove NA clades. Usually occur because of singleton clades
    group_by(generation_time) %>%
    summarise(mean.prob = mean(clade_probability, na.rm = TRUE),
              ci_hdi = bayestestR::ci(clade_probability), # 95% credible interval
              n_matches = sum(numerator_sum)) %>% 
    tidyr::unpack(ci_hdi)
}

## Tree + Line Plots
tree_line_plot = function(tree_data, posterior_summary, permuted_summary, variable){
  variable <- ensym(variable)
  p_tree = ggtree(tree_data) + 
    geom_tippoint(aes(col = factor(!!variable)), shape = "\u2588", size = 1) +
    scale_y_continuous(breaks = c(10000, 20000, 30000, 40000),
                       labels = c("10k", "20k", "30k", "40k")) +
    layout_dendrogram() + 
    theme(legend.position = "bottom", 
          axis.ticks.y = element_line(),
          axis.line.y = element_line(),
          axis.text.y = element_text())
  
  # match scales in each plot
  x_range = layer_scales(p_tree)$x$range$range
  
  p_prob = ggplot(data = posterior_summary,
                  aes(ymin = CI_low, ymax = CI_high, x = generation_time, y = mean.prob)) +
    geom_ribbon(fill = "lightblue", alpha = 0.5) +
    geom_line(col = "blue") +
    geom_ribbon(data = permuted_summary,
                aes(ymin = CI_low, ymax = CI_high, x = generation_time, y = mean.prob),
                fill = "lightgreen", alpha = 0.5) +
    geom_line(data = permuted_summary,
              aes(x = generation_time, y = mean.prob),
              col = "forestgreen") +
    ylim(c(0, 1)) +
    xlim(x_range) + ## max tree height (or max p_tree?)
    coord_flip() +
    xlab("") + 
    ylab("") + 
    # scale_x_continuous(expand = c(0, 0)) + 
    # scale_y_continuous(expand = c(0, 0)) + 
    theme_classic() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", linewidth=0.5))
  
  plot = p_tree + p_prob + 
    plot_layout(widths = c(1, 0.25))
  
  return(plot)
}

# Distance plot
Dist_Plot = function(posterior_summary, distance_summary, facet_names){
  ggplot(data = distance_summary,
         aes(
           ymin = CI_low,
           ymax = CI_high,
           x = generation_time,
           y = mean.prob
         )) +
    geom_ribbon(fill = "lightgreen", alpha = 0.5) +
    geom_line(data = distance_summary,
              aes(x = generation_time, y = mean.prob),
              col = "forestgreen") +
    facet_grid(~ distance, labeller = as_labeller(facet_names)) +
    geom_ribbon(
      data = posterior_summary,
      aes(
        ymin = CI_low,
        ymax = CI_high,
        x = generation_time,
        y = mean.prob
      ),
      fill = "lightblue",
      alpha = 0.5
    ) + geom_line(data = posterior_summary,
                  aes(x = generation_time,
                      y = mean.prob),
                  col = "blue") +
    scale_x_reverse(breaks = c(0, 10000, 20000, 30000, 40000, 50000, 60000, 70000),
                    labels = c("0", "10k", "20k", "30k", "40k", "50k", "60k", "70k")) + 
    xlab("Time (years)") + 
    ylab("Probability") + 
    theme_classic()
}

# Posterior average
posterior_aggregator = function(language_probability){
  language_probability = data.table(language_probability)
  
  numerator = language_probability[, .(numerator_sum = sum(numerator)), by = list(posterior, generation)]
  denominator = language_probability[, .(unique_denominator = unique(denominator)), 
                                     by = list(posterior, generation, clade)]
  denominator = denominator[, .(denominator_sum = sum(unique_denominator)), 
                            by = list(posterior, generation)]
  
  language_probability_agg = numerator[denominator, on = list(posterior, generation),]
  language_probability_agg[, clade_probability := numerator_sum / denominator_sum]
  
  return(language_probability_agg)
}
