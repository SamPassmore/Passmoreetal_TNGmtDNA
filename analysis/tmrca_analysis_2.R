suppressPackageStartupMessages({
  library(tmrca)
  library(ape)
  library(dplyr)
  library(ggplot2)
  library(ggridges)
  library(ggtree)
  library(phangorn)
  library(treeio)
})

## 
trees = read.nexus("analysis/darum_strict/darum_strict_posterior.trees")
trees = trees[1:200]

genetic_samples = read.csv('processed_data/genetic_samples.csv')

# change genders
genetic_samples$gender = toupper(genetic_samples$gender)

# reduce to m haploytype
genetic_samples = genetic_samples[genetic_samples$haplogroup_MN == "M",]

# add a TNG variable
genetic_samples$transnew_guinea = dplyr::if_else(genetic_samples$Family_ID == "nucl1709", "TNG", "Not TNG", missing = NA)
# add papuan variable 
genetic_samples$papuan = dplyr::if_else(genetic_samples$Family_ID %in% c("aust1307", "indo1319", ""), "Not Papuan", "Papuan", missing = NA)

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

# remove problematic value
states = states[names(states) != '575_M4"67']

x = assertthat::assert_that(any(names(states) != '575_M4"67'))

## Prune the tree 
trees = lapply(trees, function(t) keep.tip(phy = t, tip = names(states)))
class(trees) = "multiPhylo"

## Perform tmrca
fit = tmrca::tmrca(trees = trees, states = states)

n_trees = length(trees)

## stable clades
stable_clades = fit %>% 
  group_by(tips) %>% 
  summarise(clade_count = length(TREE),
            clade_proportion = clade_count / n_trees,
            mean_Date = mean(Date),
            mode_node = mode(factor(Node)),
            state = first(state)) %>%
  filter(clade_proportion > 0.5)

stable_clades$n_tips = sapply(stable_clades$tips, function(x) length(unlist(strsplit(x, "; "))))
# hist(stable_clades$n_tips)

# sum(stringr::str_detect(stable_clades$tips, "NV10"))

#### Plots for Ray

mcct_tree = read.beast("analysis/darum_strict/daru_m_hkyg_mcct.trees")
plot_tree = mcct_tree@phylo

states_df = data.frame(states = states, label = names(states))
states_df = left_join(states_df, genetic_samples[,c("ID", "gender")], by = c("label" = "ID"))

plot_tree = keep.tip(plot_tree, states_df$label)

nodes = sapply(stable_clades$tips, function(x) phytools::findMRCA(plot_tree, tips = unlist(strsplit(x, "; "))))
nodes_df = data.frame(nodes = nodes, type = stable_clades$state)


tree_df = tidytree::left_join(plot_tree, states_df, by = "label") 

## TRMCA plot
p1 = ggtree(tree_df, layout="roundrect", ladderize = TRUE, mrsd = "2023-01-01") + 
  geom_tippoint(aes(colour = states, shape = gender)) +
  geom_hilight(data=nodes_df, aes(node = nodes, fill = type),
               type = "roundrect",
               to.bottom = FALSE) +
  labs(title = "M Haplotypes") + 
  theme_classic(base_size = 20) + 
  theme_tree2()
p1
ggsave(plot = p1, filename = "figures/tmrca_m_allsamples.png", height = 290, units = "mm")

## Node stability 
# height_idx = mcct_tree@data
# height_idx$length_0.95_HPD[height_idx$length < 40000] = NA
# mcct_tree@data = height_idx

p2 = ggtree(mcct_tree, layout="roundrect", ladderize = TRUE, mrsd = "2023-01-01") + 
  labs(title = "M Haplotypes") + 
  geom_range(range='length_0.95_HPD', color='red', alpha=.6, size=2) + 
  theme_tree2()

p2

## Daru Samples
tree_df = tidytree::left_join(mcct_tree, genetic_samples, by = c("label" = "ID"))

p1.1 = ggtree(tree_df, layout="roundrect", ladderize = TRUE, mrsd = "2023-01-01") + 
  geom_tippoint(aes(colour = daru, shape = gender)) +
  scale_color_manual(values = c("red", "transparent")) + 
  # geom_hilight(data=nodes_df, aes(node = nodes, fill = type),
  #              type = "roundrect",
  #              to.bottom = FALSE) +
  labs(title = "N Haplotypes") + 
  theme_tree2() + 
  theme(text=element_text(size=21), legend.position = "left") 
p1.1
ggsave(plot = p1.1, filename = "figures/tmrca_m_allsamples_highlightDARU.png", height = 290, units = "mm")


#### For Daru subsample
daru_id = genetic_samples$ID[genetic_samples$daru == "Daru"]
rem_id = which(daru_id == '575_M4"67')
daru_id = daru_id[-rem_id]

daru_tree = ape::keep.tip(phy = mcct_tree@phylo, tip = daru_id)
states = genetic_samples$of_interest[genetic_samples$daru == "Daru"]
states = states[-rem_id]
names(states) = daru_id

## All samples speak non-tng languages so states are unnecessary
daru_md = genetic_samples[genetic_samples$daru == "Daru",]
daru_md$label = daru_md$ID
daru_treedf = tidytree::left_join(daru_tree, daru_md, by = c("label"))

p3 = ggtree(daru_treedf, layout = "roundrect", ladderize = TRUE,  mrsd = "2023-01-01") + 
  geom_tippoint(aes(colour = Family_ID, shape = gender), size = 4) +
  geom_tiplab(aes(label = Name), size = 6, offset = 300) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 120, 6, 6)) + 
  theme(text=element_text(size=21), legend.position = "left")
p3
ggsave(plot = p3, filename = "figures/tmrca_m_daru.png", height = 290, units = "mm")
