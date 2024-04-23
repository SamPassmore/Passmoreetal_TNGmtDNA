suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
  library(ggtree)
  library(dplyr)
  library(stringr)
  library(assertthat)
})

# tree_m = read.nexus("analysis/m_relaxed_262E8_empirical/BEAST_2_m-papua_m_MCCT")
tree_m = read.nexus("analysis/daru_m_strictprior/papuanm_daru_strictprior-papua_m_daru.mcct.trees")
tree_m = phytools::read.newick("analysis/daru_m_strictprior/papuanm_daru_strictprior-papua_m_daru.mcct.trees")
tree_m = ladderize(tree_m)

tree_n = read.nexus("analysis/daru_n_stricprior/")
tree_n = ladderize(tree_n)

genetic_samples = read.csv("processed_data/genetic_samples.csv")
genetic_samples$ID = stringr::str_replace_all(genetic_samples$ID, pattern = '[\\*\\+\\"]', "_")

## Add language variable
# add a new guinea variable
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

x = assert_that(
  all(genetic_samples$Family_ID[genetic_samples$assigned_language == "Pawaia [pwa]"] == "pawa1255")
)

x = assert_that(
  all(genetic_samples$Family_ID[genetic_samples$assigned_language == "Pele-Ata"] == "pele1245")
)


my_phylogeny_plot = function(samples, tree, haplotype, prune = FALSE, title = ""){
  # Only look at M haplotypes
  samples = samples %>% 
    filter(haplogroup_MN == haplotype)
  
  ## I have removed some languages in the pipeline.
  ## These will filter out eventually, but I force them out for now
  samples = samples[!is.na(samples$ID),]
  # prune the tree
  

  if(prune){
    tree = keep.tip(tree, samples$ID[!is.na(samples$of_interest)])  
  } else {
    tree = keep.tip(tree, samples$ID)  
  }
    
  rownames(samples) = samples$ID
  samples = samples[tree$tip.label,]
  
  if(!all(rownames(samples) == tree$tip.label)){
    stop("There is a matching problem")
  }
  

  
  samples$label = tree$tip.label
  tree_df = tidytree::full_join(tree, samples, by = "label")
  
  ggtree(tree_df, layout="roundrect", ladderize = TRUE) + 
    geom_tippoint(aes(colour = of_interest)) + 
    theme(legend.position = "bottom") + 
    labs(title = title) 
}


my_phylogeny_plot(genetic_samples, tree_m, "M", title = "All M haplotype mtDNA samples")
my_phylogeny_plot(genetic_samples, tree_m, "M", prune = TRUE, "M haplotype mtDNA samples with language data")
  

my_phylogeny_plot(genetic_samples, tree_n, "N", title = "All N haplotype mtDNA samples")
my_phylogeny_plot(genetic_samples, tree_n, "N", prune = TRUE, "N haplotype mtDNA samples with language data")

## Language breakdown
genetic_samples %>% 
  group_by(haplogroup_MN, of_interest) %>% 
  summarise(n = n()) %>% 
  mutate(prop = n / sum(n)) 

genetic_samples %>% 
  group_by(haplogroup_MN, assigned_language) %>% 
  summarise(n = n(), of_interest = first(of_interest)) %>% 
  mutate(prop = n / sum(n)) %>% 
  filter(is.na(of_interest)) %>%
  View()
  

pdf("figures/papuandaru_Mhaplotype_iqtree.pdf")
my_phylogeny_plot(genetic_samples, tree_m, "M", prune = TRUE, "M haplotype mtDNA samples with language data") 
dev.off()

pdf("figures/papuandaru_Mhaplotype_iqtree_labels.pdf", height = 60)
my_phylogeny_plot(genetic_samples, tree_m, "M", prune = TRUE, "M haplotype mtDNA samples with language data") + 
  geom_tiplab(aes(label = assigned_language_code)) 
dev.off()

