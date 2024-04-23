# create language metadata
library(dplyr)

genetic_samples = read.csv('data/genetic_samples.csv')

speaker_counts = genetic_samples %>% 
  group_by(assigned_language) %>% 
  summarise(glottocode = first(Glottocode),
            speakers = n(),
            language_group = first(LANG_GROUP),
            language_phylum = first(LANG_PHYLU),
            latitude = first(Point_Y),
            longitude = first(Point_X),
            in_kinbank = first(in_kinbank),
            in_pronouns = first(in_pronouns),
            in_grambank = first(in_grambank))
            
write.csv(speaker_counts, 'data/language_metadata.csv')
