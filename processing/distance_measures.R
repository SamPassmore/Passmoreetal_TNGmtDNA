library(terra)
library(raster)
library(gdistance)
library(reshape2)
library(ggplot2)

r = raster("../../newguinea_linguisticdiversity/processed_data/data_layers/new_guinea_elevation.tiff")

metadata = read.csv("processed_data/ArcGIS_LOCATIONS.csv")
# Remove missing values
metadata = metadata[complete.cases(metadata),]
xy = metadata[,c("longitude_hybrid", "latitude_hybrid")]

sp = SpatialPoints(coords = xy, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))


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
dimnames(time_h) = list(metadata$ID, metadata$ID)

saveRDS(time_h, file = "processed_data/hikingtime_h.RDS")

## Testing
library(geosphere)
distance_matrix = geosphere::distm(metadata[,c("longitude_hybrid", "latitude_hybrid")], fun = distHaversine) / 1000

A = unlist(metadata[metadata$ID == "OV0561",3:2])
B = unlist(metadata[metadata$ID == "EGAN00001279046",3:2])
C = unlist(metadata[metadata$ID == "PNG111",3:2])

csA = gdistance::accCost(conduct, A)
csB = gdistance::accCost(conduct, B)

AtoB <- shortestPath(conduct, A, B, output="SpatialLines")
BtoA <- shortestPath(conduct, B, A, output="SpatialLines")
AtoC <- shortestPath(conduct, A, C, output="SpatialLines")

plot(r, xlab="x coordinate (m)", ylab="y coordinate (m)",
     legend.lab="Altitude (masl)", ylim = c(-7, -5), xlim = c(142, 146))
lines(AtoB, col="red", lwd=2)
lines(BtoA, col="blue")
lines(AtoC, col="darkgreen")
points(x = metadata$longitude_hybrid, metadata$latitude_hybrid)

