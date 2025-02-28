library(terra)
library(dplyr)
library(raster)
library(gdistance)
library(reshape2)
library(ggplot2)

r = raster("../../newguinea_linguisticdiversity/processed_data/data_layers/new_guinea_elevation.tiff")

# Read TNG_metadata and re-load dups
TNG_metadata = read.csv("processed_data/TNGData.csv")
TNG_metadata = TNG_metadata[TNG_metadata$interview_language != "Tok Pisin [tpi]",]
# Baining speakers are unlikely to have come from New Guinea Mainland. This is a TNG_metadata mistake. 
TNG_metadata = TNG_metadata[TNG_metadata$Branch_Family != "Baining",]

## Add in duplicate TNG_metadata
TNG_dups = read.csv("processed_data/genetic_duplicate_data.csv")
TNG_dups = TNG_dups[!TNG_dups$PhyID %in% TNG_metadata$PhyID,] # subset to the duplicates only

TNG_metadata = rbind(TNG_metadata, TNG_dups)

# get variables of interest
TNG_metadata = TNG_metadata %>% select(PhyID, ID, longitude_hybrid, latitude_hybrid)
# Remove missing values
TNG_metadata = TNG_metadata[complete.cases(TNG_metadata),]
xy = TNG_metadata[,c("longitude_hybrid", "latitude_hybrid")]

sp = SpatialPoints(coords = xy, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

# calculate the height difference between cells
altDiff <- function(x){x[2] - x[1]}
hd <- transition(r$wc2.1_2.5m_elev, altDiff, 8, symm=FALSE)
slope <- geoCorrection(hd, scl=FALSE) # divide transition by distances between cells

# Calculate speed
adj <- adjacent(x=r$wc2.1_2.5m_elev, cells=1:ncell(r), direction=8)
speed <- slope
speed[adj] <-  exp(-3.5 * abs(slope[adj] + 0.05)) # m/s

conduct <- geoCorrection(speed, scl=FALSE) # conductance = speed/distance = 1/travel time

time_s = gdistance::costDistance(conduct, sp)
time_h = time_s / 3600
dimnames(time_h) = list(TNG_metadata$ID, TNG_metadata$ID)

saveRDS(time_h, file = "processed_data/hikingtime_h.RDS")

## Testing
# library(geosphere)
# distance_matrix = geosphere::distm(TNG_metadata[,c("longitude_hybrid", "latitude_hybrid")], fun = distHaversine) / 1000
# 
# A = unlist(TNG_metadata[TNG_metadata$ID == "OV0561",c("longitude_hybrid", "latitude_hybrid")])
# B = unlist(TNG_metadata[TNG_metadata$ID == "EGAN00001279046",c("longitude_hybrid", "latitude_hybrid")])
# C = unlist(TNG_metadata[TNG_metadata$ID == "PNG111",c("longitude_hybrid", "latitude_hybrid")])
# 
# csA = gdistance::accCost(conduct, A)
# csB = gdistance::accCost(conduct, B)
# 
# AtoB <- shortestPath(conduct, A, B, output="SpatialLines")
# BtoA <- shortestPath(conduct, B, A, output="SpatialLines")
# AtoC <- shortestPath(conduct, A, C, output="SpatialLines")
# 
# plot(r, xlab="x coordinate (m)", ylab="y coordinate (m)",
#      legend.lab="Altitude (masl)", ylim = c(-7, -5), xlim = c(142, 146))
# lines(AtoB, col="red", lwd=2)
# lines(BtoA, col="blue")
# lines(AtoC, col="darkgreen")
# points(x = TNG_metadata$longitude_hybrid, TNG_metadata$latitude_hybrid)

