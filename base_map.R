## Base map of PNG

library(terra)
library(tidyterra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(ggspatial)

indonesia = ne_states(country = "Indonesia", returnclass = "sf")
# only include West Papua
west_papua = indonesia[indonesia$name %in% c("Papua", "Papua Barat"),]

# Map of PNG
png = ne_states(country = "Papua New Guinea", returnclass = "sf")

# Merge
new_guinea_sf = rbind(west_papua, png)

base_map = ggplot(st_union(new_guinea_sf)) + geom_sf() + theme_classic()

save(base_map, file = "figures/base_map.Rdata")
