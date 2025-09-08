# helper script for Supplemantary materials

summarizer = function(posterior, permuted_results = FALSE){
  posterior$generation = as.numeric(gsub("g_", replacement = "", posterior$generation))
  
  if(permuted_results){
    out = posterior %>%
      filter(!is.na(clade_probability)) %>% # remove NA clades. Usually occur because of singleton clades
      # group_by(generation, posterior) %>%
      # summarise(wntree_cladeprob = mean(clade_probability, na.rm = TRUE)) %>% 
      # ungroup() %>% 
      group_by(generation) %>% 
      # summarise(mean_clade_probability = mean(wntree_cladeprob, na.rm = TRUE),
      #           ci_hdi = bayestestR::ci(wntree_cladeprob)) %>% 
      summarise(mean_clade_probability = mean(clade_probability, na.rm = TRUE),
                ci_hdi = bayestestR::ci(clade_probability, method = "HDI", ci = 0.99)) %>% 
      tidyr::unpack(ci_hdi)  
  } else {
    out = posterior %>%
      filter(!is.na(clade_probability)) %>% # remove NA clades. Usually occur because of singleton clades
      group_by(generation) %>%
      summarise(mean_clade_probability = mean(clade_probability, na.rm = TRUE),
                ci_hdi = bayestestR::ci(clade_probability, method = "HDI", ci = 0.99), # 95% credible interval
                n_matches = sum(numerator_sum)) %>% 
      tidyr::unpack(ci_hdi)  
  }
  out
}

## Tree + Line Plots
tree_line_plot = function(tree_data, posterior_summary, permuted_summary, variable){
  variable <- ensym(variable)
  p_tree = ggtree(tree_data, size = 0.05) + 
    geom_tippoint(aes(col = factor(!!variable)), shape = "I", size = 4) +
    scale_fill_discrete(name = NULL) + 
    scale_color_manual(values = c("#FFCD00", "#C8102E")) +
    geom_vline(xintercept = -10000, linetype = "dashed") +
    annotate("text", x = -8000, y = 50, label = "10kya: Evidence of Agriculture", size = 3) + 
    layout_dendrogram() + 
    theme(legend.position = "bottom", 
          axis.ticks.y = element_line(),
          axis.line.y = element_line(),
          axis.text.y = element_text(),
          panel.grid.major.x = element_blank())
  
  # match scales in each plot
  x_range = layer_scales(p_tree)$x$range$range
  
  p_prob = ggplot(data = posterior_summary,
                  aes(ymin = CI_low, ymax = CI_high, x = generation, y = mean_clade_probability)) +
    geom_ribbon(fill = "#F8D1D7", alpha = 0.5) +
    geom_line(col = "#C8102E") +
    geom_ribbon(data = permuted_summary,
                aes(ymin = CI_low, ymax = CI_high, x = generation, y = mean_clade_probability),
                fill = "lightblue", alpha = 0.5) +
    geom_line(data = permuted_summary,
              aes(x = generation, y = mean_clade_probability),
              col = "#00274D") +
    ylim(c(0.3, 1.0)) +
    xlim(x_range) + ## max tree height (or max p_tree?)
    coord_flip() +
    xlab("") + 
    ylab("") + 
    scale_fill_discrete(name = NULL) + 
    theme_classic(base_size = 10) + 
    geom_vline(xintercept = 10000, linetype = "dashed") +
    theme(panel.grid.major.y = element_line(), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", linewidth=0.5),
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  plot = p_tree + p_prob + 
    plot_layout(guides = "collect", widths = c(1, 0.25)) &
    theme(legend.position = "bottom", legend.title = element_blank())  
  
  return(plot)
}

# Distance plot
Dist_Plot = function(posterior_summary, distance_summary, facet_names){
  ggplot(data = distance_summary,
         aes(
           ymin = CI_low,
           ymax = CI_high,
           x = generation,
           y = mean_clade_probability
         )) +
    geom_ribbon(fill = "lightgreen", alpha = 0.5) +
    geom_line(data = distance_summary,
              aes(x = generation, y = mean_clade_probability),
              col = "forestgreen") +
    facet_wrap(~ distance, labeller = as_labeller(facet_names), ncol = 2) +
    geom_ribbon(
      data = posterior_summary,
      aes(
        ymin = CI_low,
        ymax = CI_high,
        x = generation,
        y = mean_clade_probability
      ),
      fill = "lightblue",
      alpha = 0.5
    ) + geom_line(data = posterior_summary,
                  aes(x = generation,
                      y = mean_clade_probability),
                  col = "blue") +
    scale_x_reverse(breaks = c(0, 10000, 20000, 30000, 40000, 50000, 60000, 70000),
                    labels = c("0", "10k", "20k", "30k", "40k", "50k", "60k", "70k")) + 
    xlab("Time (years)") + 
    ylab("Probability") + 
    theme_classic()  + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
    
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

mainplot_prep = function(){
  posterior = posterior[posterior$generation < mcct_depth,]
  permutations = permutations[permutations$generation < mcct_depth,]

  posterior_summary = summarizer(posterior)
  permutation_summary_p = summarizer(permutations)
  
  list(posterior = posterior_summary, permutation = permutation_summary_p)
}
