## Haplogroup table
suppressPackageStartupMessages({
  library(ape)
  library(forcats)
  library(dplyr)
  library(kableExtra)
  library(knitr)
  library(janitor)
})

options(knitr.kable.NA = '-')

# Load Metadata
TNG_metadata = read.csv("processed_data/TNGData.csv")

## Add in duplicate metadata
TNG_dups = read.csv("processed_data/genetic_duplicate_data.csv")
TNG_dups = TNG_dups[!TNG_dups$PhyID %in% TNG_metadata$PhyID,] # subset to the duplicates only

TNG_metadata = rbind(TNG_metadata, TNG_dups)

TNG_metadata = TNG_metadata[TNG_metadata$interview_language != "Tok Pisin [tpi]",] 
# Baining speakers are unlikely to have come from New Guinea Mainland. This is a metadata mistake. 
TNG_metadata = TNG_metadata[TNG_metadata$Branch_Family != "Baining",]

# remove data not used in the phylogeny
mcct_tree = read.nexus("/Volumes/PassmoreSSD/PapuanPast/20250106_LanguageSample_CBES_Relaxed_Constrained/20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.mcct.duptips.trees") 

TNG_metadata = TNG_metadata[TNG_metadata$PhyID %in% mcct_tree$tip.label,]

# Summarise Haplogroups
TNG_metadata$HG = fct_collapse(TNG_metadata$haplogroup,
                               `M/E` = c( "M27a2a", "M27b1", "M28a+204",
                                          "M29a", "E1a1b4", "E1a2+(16261)"),
                               N13 = c("N13"),
                               P = c("P"),
                               P1 = c("P1", "P1+152"),
                               P1d = c("P1d", "P1d1", "P1d1a", "P1d2", "P1d2a"),
                               P2 = c("P2"),
                               `P3/4` = c("P3b", "P3b1", "P4a", "P4a1"),
                               Q1 = c("Q1", "Q1+@16223"), 
                               `Q1a/d` = c("Q1a", "Q1a1", "Q1a1a", "Q1d"),
                               `Q2/3` = c("Q2a", "Q2a4", "Q3", "Q3a+61  62", 
                                          "Q3a+61 62", "Q3a1", "Q3b"))

# Proportions and differences of TNG
# round(prop.table(table(TNG_metadata$TNG, TNG_metadata$HG), margin = 1), 3)
# sort(abs(round(prop.table(table(TNG_metadata$TNG, TNG_metadata$HG), margin = 1), 3)[1,] - round(prop.table(table(TNG_metadata$TNG, TNG_metadata$HG), margin = 1), 3)[2,]))

haplo_table = TNG_metadata %>% 
  # dplyr::filter(!Branch_Family %in% c("Eleman", "Mailuan", "Ok-Oksapmin", "Sepik",
  #                                     "South Bougainville", "Wiru", "Madang", "Ndu",
  #                                     "Indo-European")) %>% 
  dplyr::filter(!Branch_Family %in% c("Indo-European")) %>% 
  dplyr::filter(!is.na(Branch_Family)) %>% 
  group_by(Branch_Family, HG) %>% 
  mutate(n=n()) %>% 
  group_by(HG) %>% 
  distinct(HG, Branch_Family, n) %>%
  mutate(Per=n/sum(n)) %>% #, np=paste0(n," (",round(Per*100,0)," %)")) %>%
  select(-Per) %>% 
  tidyr::spread(HG, n) 

non_tng = c("Arafundi", "Austronesian", "Bosavi", "Eastern Trans Fly", "Eleman", "Lower Sepik",
            "Ndu", "Pahoturi", "Sepik", "South Bougainville", "Yam"
)
tng = c("Angan", "Chimbu-Wahgi", "Enga-Kewa-Huli", "Kainantu-Goroka", "Kiwaian", "Madang",
        "Mailuan", "Ok-Oksapmin", "Suki-Gogodala", "Wiru")


css_row_height <- "
    <style>
    table {
      border-collapse: collapse; /* Optional: collapse borders for tighter look */
    }
    
    table tr {
      height: 20px; /* Adjust as needed */
    }
    </style>
    "

ht = haplo_table %>% 
  arrange(factor(Branch_Family, levels = c(tng, non_tng))) %>% 
  adorn_totals("row") %>%  # Add total row
  adorn_totals("col") %>%  # Add total column
  rename(" "= Branch_Family) %>% 
  kable(., escape = FALSE, padding=0L, linesep = "") %>%
  kableExtra::group_rows("Trans-New Guinea", 1, 10) %>%
  kableExtra::group_rows("Not Trans-New Guinea", 11, nrow(haplo_table)) %>% 
  kable_styling(font_size = 9, full_width = FALSE) %>% 
  paste0(css_row_height, sep = "\n") 


rio::export(ht, file = "figures/haplogroup_table.html", format = "html")

## Chi Sq Tests
tng = TNG_metadata$TNG
p1d = TNG_metadata$HG == "P1d"
tng_p1d = table(tng, p1d)
dimnames(tng_p1d) = list(c("TNG", "NotTNG"), c("P1d", "NotP1d"))
chisq.test(tng_p1d)

tng = TNG_metadata$TNG
q1 = TNG_metadata$HG == "Q1"
tng_q1 = table(tng, q1)
dimnames(tng_q1) = list(c("TNG", "NotTNG"), c("Q1", "Notq1"))
chisq.test(tng_q1)
