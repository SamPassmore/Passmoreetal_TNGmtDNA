## TRMCA code

library(tmrca)
library(ape)
library(dplyr)
library(ggplot2)
library(ggridges)
library(ggtree)
library(phangorn)

## 
trees = read.nexus("analysis/m_relaxed_262E8_empirical/m_962trees.tree")

genetic_samples = read.csv('processed_data/genetic_samples.csv')

# reduce to m haploytype
genetic_samples = genetic_samples[genetic_samples$haplogroup_MN == "M",]
# add a TNG variable
genetic_samples$transnew_guinea = ifelse(genetic_samples$Family_ID == "nucl1709", "TNG", "Not TNG")
# add papuan variable 
genetic_samples$papuan = ifelse(genetic_samples$Family_ID %in% c("aust1307", "indo1319", ""), "Not Papuan", "Papuan")

# samples of interest variable
genetic_samples$of_interest = ifelse(genetic_samples$papuan == "Papuan" & genetic_samples$transnew_guinea == "TNG",
                                     "Papuan & TNG", 
                                     ifelse(genetic_samples$papuan == "Papuan" & genetic_samples$transnew_guinea != "TNG",
                                            "Papuan & Not TNG", NA))

genetic_samples$of_interest[genetic_samples$Family_ID == "aust1307"] = "Austronesian"
genetic_samples$of_interest[genetic_samples$Family_ID == "indo1319"] = "Indo-European"
genetic_samples$of_interest[genetic_samples$Family_ID == ""] = NA

states = genetic_samples$of_interest
names(states) = genetic_samples$ID

states = states[!is.na(states)]

## Prune the tree 
trees = lapply(trees, function(t) keep.tip(phy = t, tip = names(states)))
class(trees) = "multiPhylo"

## Perform tmrca
fit = tmrca::tmrca(trees = trees, states = states)

# Clade stability for only women 
women_samples = genetic_samples %>% 
  filter(gender == "F") %>% 
  dplyr::select(ID, of_interest) %>% 
  na.omit() 
  
states_w = women_samples$of_interest
names(states_w) = women_samples$ID

tree_w = lapply(trees, function(t) keep.tip(phy = t, tip = names(states_w)))
class(tree_w) = "multiPhylo"

fit_w = tmrca::tmrca(trees = tree_w, states = states_w)

# Clade stability for only men
men_samples = genetic_samples %>% 
  filter(gender == "M") %>% 
  dplyr::select(ID, of_interest) %>% 
  na.omit() 

states_m = men_samples$of_interest
names(states_m) = men_samples$ID

tree_m= lapply(trees, function(t) keep.tip(phy = t, tip = names(states_m)))
class(tree_m) = "multiPhylo"

fit_m = tmrca::tmrca(trees = tree_m, states = states_m)


#### Plot Results ####
# How many clades exist in more than 50% of all trees
n_trees = length(trees)
      
mode <- function(codes){
  which.max(tabulate(codes))
}

stable_clades = fit %>% 
  group_by(tips) %>% 
  summarise(clade_count = length(TREE),
            clade_proportion = clade_count / n_trees,
            mean_Date = mean(Date),
            mode_node = mode(factor(Node))) %>%
  filter(clade_proportion > 0.5)


stable_clades_w = fit_w %>% 
  group_by(tips) %>% 
  summarise(clade_count = length(TREE),
            clade_proportion = clade_count / n_trees,
            mean_Date = mean(Date),
            mode_node = mode(factor(Node))) %>%
  filter(clade_proportion > 0.5)

stable_clades_m = fit_m %>% 
  group_by(tips) %>% 
  summarise(clade_count = length(TREE),
            clade_proportion = clade_count / n_trees,
            mean_Date = mean(Date),
            mode_node = mode(factor(Node))) %>%
  filter(clade_proportion > 0.5)

density_plots = function(fit, n_trees){
  stable_clades = fit %>% 
    group_by(tips) %>% 
    summarise(clade_count = length(TREE),
              clade_proportion = clade_count / n_trees,
              mean_Date = mean(Date),
              mode_node = mode(factor(Node))) %>%
    filter(clade_proportion > 0.5)
  
  stable_clades_df = fit[fit$tips %in% stable_clades$tips,]
  
  ggplot(stable_clades_df, aes(y = reorder(tips, Date, mean), x = Date, fill = state)) +
    geom_density_ridges2(
      rel_min_height = 0.01,
      scale = 2,
      alpha = 0.6,
      jittered_points = TRUE,
      position = position_points_jitter(width = 0.05, height = 0),
      point_shape = '|',
      point_size = 0.75,
      point_alpha = 1,
    )
}


p1 = density_plots(fit, n_trees)
p_w = density_plots(fit_w, length(tree_w))
p_m = density_plots(fit_m, length(tree_m))

ggsave(plot = p1, filename = "figures/tmrca_m.png")
ggsave(plot = p_w, filename = "figures/tmrca_m.png")
ggsave(plot = p_m, filename = "figures/tmrca_m.png")

### Tree plot ####

mcct_tree = read.nexus("analysis/m_relaxed_262E8_empirical/BEAST_2_m-papua_m_MCCT_subset")

tree_plot = function(mcct_tree, states, data, stable_clades){
  mcct_tree = keep.tip(phy = mcct_tree, tip = names(states))
  states_df = data.frame(states = states, label = names(states))
  states_df = left_join(states_df, data[,c("ID", "gender")], by = c("label" = "ID"))
  
  nodes = sapply(stable_clades$tips, function(x) phytools::findMRCA(mcct_tree, tips = unlist(strsplit(x, "; "))))
  nodes_df = data.frame(nodes = nodes, type = "red")
  
  tree_df = tidytree::full_join(mcct_tree, states_df, by = "label")
  
  ggtree(tree_df, layout="circular", ladderize = TRUE) + 
    geom_tippoint(aes(colour = states, shape = gender)) +
    geom_hilight(data=nodes_df, aes(node=nodes), fill = "red",
                 type = "roundrect") +
    labs(title = "M Haplotypes") 
}

p2 = tree_plot(mcct_tree, states, genetic_samples, stable_clades)
p_w = tree_plot(mcct_tree, states_w, genetic_samples, stable_clades_w)
p_m = tree_plot(mcct_tree, states_m, genetic_samples, stable_clades_m)

ggsave(plot = p2, filename = "figures/tmrca_m_phylo.pdf", width = 300, height = 300, units = "mm")
ggsave(plot = p_w, filename = "figures/tmrca_m_phylo.pdf", width = 300, height = 300, units = "mm")
ggsave(plot = p_m, filename = "figures/tmrca_m_phylo.pdf", width = 300, height = 300, units = "mm")
