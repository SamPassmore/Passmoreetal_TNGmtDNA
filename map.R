library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)

languages = read.csv('data/language_metadata.csv')

world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot(data = world) +
  geom_sf() + 
  xlim(c(140, 156)) + 
  ylim(c(-12, -2)) + 
  geom_jitter(data = languages, 
              aes(
                x = longitude, 
                y = latitude,
                size = speakers,
                col = )) + 
  theme_minimal()

genetic_samples = read.csv("processed_data/genetic_samples.csv")
genetic_samples$transnew_guinea = ifelse(genetic_samples$Family_ID == "nucl1709", "TNG", "Not TNG")

ggplot(data = world) +
  geom_sf() + 
  xlim(c(140, 156)) + 
  ylim(c(-12, -2)) + 
  geom_jitter(data = genetic_samples, 
              aes(
                x = longitude_ofsampling, 
                y = latitude_ofsampling,
                col = transnew_guinea)) + 
  theme_minimal()

ggplot(data = world) +
  geom_sf() + 
  xlim(c(140, 156)) + 
  ylim(c(-12, -2)) + 
  geom_jitter(data = genetic_samples, 
              aes(
                x = longitude_ofsampling, 
                y = latitude_ofsampling,
                col = haplogroup_MN)) + 
  theme_minimal()
