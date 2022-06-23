#-------------------------------------------------------------------------------
# Name:         forest_type.R
# Author:       Florian Franz
# Description:  script calculates canopy cover and derives areas
#               of open, dense forest and gaps for three different locations:
#               Reinhardshagen, Neukirchen_8 and Neukirchen_9
# Data          nDSM tif files from platforms 
#               aircraft (0.5 m resolution) and drone (0.1 m resolution)
#-------------------------------------------------------------------------------

# Load packages
#---------------
library(terra)
library(raster)



# Data import
#-------------
# Logical switch
# Select 'DSM = TRUE' and 'nDSM = FALSE' to import DSM files
# Select 'DSM = FALSE' and 'nDSM = TRUE' to import nDSM files
# Select 'DSM = TRUE' and 'nDSM = TRUE' to import both
DSM = FALSE
nDSM = TRUE

# Source script for nDSM import and preprocessing
source("J:/scripts/src/process_nDSMs.R", local = TRUE)



# Calculate canopy cover
#------------------------
# Reclassify nDSMs:
# all pixels with a value above 3 m are defined as "forest",
# pixel value >= 3 --> 1 (forest), pixel value < 3 --> 0 (non-forest)

matr <- c(0, 3, 0,
          3, 100, 1)

rclmatr <- matrix(matr, ncol = 3, byrow = TRUE)


##############
# Improved workflow (focal part yet not implemented)
ndsms_list <- list(ndsm_drone_reinhardshagen_1, ndsm_drone_reinhardshagen_1_resampled, ndsm_aircraft_reinhardshagen_1,
                   ndsm_drone_reinhardshagen_2, ndsm_drone_reinhardshagen_2_resampled, ndsm_aircraft_reinhardshagen_2,
                   ndsm_drone_neukirchen8_1, ndsm_drone_neukirchen8_1_resampled, ndsm_aircraft_neukirchen8_1,
                   ndsm_drone_neukirchen8_2, ndsm_drone_neukirchen8_2_resampled, ndsm_aircraft_neukirchen8_2,
                   ndsm_drone_neukirchen9_1, ndsm_drone_neukirchen9_1_resampled, ndsm_aircraft_neukirchen9_1,
                   ndsm_drone_neukirchen9_2, ndsm_drone_neukirchen9_2_resampled, ndsm_aircraft_neukirchen9_2)


ndsms_rcl_list <- lapply(ndsms_list, FUN = function(x) terra::classify(x, rclmatr, right = FALSE))

# Calculate focal mean:
# Circular moving window with radius = 25 m

#ndsms_focal_list <- lapply(test_list, FUN = function(x) terra::focal(x, w = ndsms_mov_wind_list, fun = "mean", na.policy = "omit"))

ndsm_drone_reinhardshagen_1_mov_wind <- terra::focalMat(ndsms_rcl_list[[1]],
                                                        d = 25, type = "circle")

ndsm_drone_reinhardshagen_1_resampled_mov_wind <- terra::focalMat(ndsms_rcl_list[[2]],
                                                                  d = 25, type = "circle")

ndsm_aircraft_reinhardshagen_1_mov_wind <- terra::focalMat(ndsms_rcl_list[[3]],
                                                           d = 25, type = "circle")

ndsm_drone_reinhardshagen_2_mov_wind <- terra::focalMat(ndsms_rcl_list[[4]],
                                                        d = 25, type = "circle")

ndsm_drone_reinhardshagen_2_resampled_mov_wind <- terra::focalMat(ndsms_rcl_list[[5]],
                                                                  d = 25, type = "circle")

ndsm_aircraft_reinhardshagen_2_mov_wind <- terra::focalMat(ndsms_rcl_list[[6]],
                                                           d = 25, type = "circle")

ndsm_drone_neukirchen8_1_mov_wind <- terra::focalMat(ndsms_rcl_list[[7]],
                                                     d = 25, type = "circle")

ndsm_drone_neukirchen8_1_resampled_mov_wind <- terra::focalMat(ndsms_rcl_list[[8]],
                                                               d = 25, type = "circle")

ndsm_aircraft_neukirchen8_1_mov_wind <- terra::focalMat(ndsms_rcl_list[[9]],
                                                        d = 25, type = "circle")

ndsm_drone_neukirchen8_2_mov_wind <- terra::focalMat(ndsms_rcl_list[[10]],
                                                     d = 25, type = "circle")

ndsm_drone_neukirchen8_2_resampled_mov_wind <- terra::focalMat(ndsms_rcl_list[[11]],
                                                               d = 25, type = "circle")

ndsm_aircraft_neukirchen8_2_mov_wind <- terra::focalMat(ndsms_rcl_list[[12]],
                                                        d = 25, type = "circle")

ndsm_drone_neukirchen9_1_mov_wind <- terra::focalMat(ndsms_rcl_list[[13]],
                                                     d = 25, type = "circle")

ndsm_drone_neukirchen9_1_resampled_mov_wind <- terra::focalMat(ndsms_rcl_list[[14]],
                                                               d = 25, type = "circle")

ndsm_aircraft_neukirchen9_1_mov_wind <- terra::focalMat(ndsms_rcl_list[[15]],
                                                        d = 25, type = "circle")

ndsm_drone_neukirchen9_2_mov_wind <- terra::focalMat(ndsms_rcl_list[[16]],
                                                     d = 25, type = "circle")

ndsm_drone_neukirchen9_2_resampled_mov_wind <- terra::focalMat(ndsms_rcl_list[[17]],
                                                               d = 25, type = "circle")

ndsm_aircraft_neukirchen9_2_mov_wind <- terra::focalMat(ndsms_rcl_list[[18]],
                                                        d = 25, type = "circle")

mov_wind_list <- list(ndsm_drone_reinhardshagen_1_mov_wind, ndsm_drone_reinhardshagen_1_resampled_mov_wind, ndsm_aircraft_reinhardshagen_1_mov_wind,
                      ndsm_drone_reinhardshagen_2_mov_wind, ndsm_drone_reinhardshagen_2_resampled_mov_wind, ndsm_aircraft_reinhardshagen_2_mov_wind,
                      ndsm_drone_neukirchen8_1_mov_wind, ndsm_drone_neukirchen8_1_resampled_mov_wind, ndsm_aircraft_neukirchen8_1_mov_wind,
                      ndsm_drone_neukirchen8_2_mov_wind, ndsm_drone_neukirchen8_2_resampled_mov_wind, ndsm_aircraft_neukirchen8_2_mov_wind,
                      ndsm_drone_neukirchen9_1_mov_wind, ndsm_drone_neukirchen9_1_resampled_mov_wind, ndsm_aircraft_neukirchen9_1_mov_wind,
                      ndsm_drone_neukirchen9_2_mov_wind, ndsm_drone_neukirchen9_2_resampled_mov_wind, ndsm_aircraft_neukirchen9_2_mov_wind)



ndsm_drone_reinhardshagen_1_foc <- terra::focal(ndsms_rcl_list[[1]],
                                                w = mov_wind_list[[1]],
                                                fun = "mean", na.policy = "omit")

ndsm_drone_reinhardshagen_1_resampled_foc <- terra::focal(ndsms_rcl_list[[2]],
                                                          w = mov_wind_list[[2]],
                                                          fun = "mean", na.policy = "omit")

ndsm_aircraft_reinhardshagen_1_foc <- terra::focal(ndsms_rcl_list[[3]],
                                                   w = mov_wind_list[[3]],
                                                   fun = "mean", na.policy = "omit")

ndsm_drone_reinhardshagen_2_foc <- terra::focal(ndsms_rcl_list[[4]],
                                                w = mov_wind_list[[4]],
                                                fun = "mean", na.policy = "omit")

ndsm_drone_reinhardshagen_2_resampled_foc <- terra::focal(ndsms_rcl_list[[5]],
                                                          w = mov_wind_list[[5]],
                                                          fun = "mean", na.policy = "omit")

ndsm_aircraft_reinhardshagen_2_foc <- terra::focal(ndsms_rcl_list[[6]],
                                                   w = mov_wind_list[[6]],
                                                   fun = "mean", na.policy = "omit")

ndsm_drone_neukirchen8_1_foc <- terra::focal(ndsms_rcl_list[[7]],
                                             w = mov_wind_list[[7]],
                                             fun = "mean", na.policy = "omit")

ndsm_drone_neukirchen8_1_resampled_foc <- terra::focal(ndsms_rcl_list[[8]],
                                                       w = mov_wind_list[[8]],
                                                       fun = "mean", na.policy = "omit")

ndsm_aircraft_neukirchen8_1_foc <- terra::focal(ndsms_rcl_list[[9]],
                                                w = mov_wind_list[[9]],
                                                fun = "mean", na.policy = "omit")

ndsm_drone_neukirchen8_2_foc <- terra::focal(ndsms_rcl_list[[10]],
                                             w = mov_wind_list[[10]],
                                             fun = "mean", na.policy = "omit")

ndsm_drone_neukirchen8_2_resampled_foc <- terra::focal(ndsms_rcl_list[[11]],
                                                       w = mov_wind_list[[11]],
                                                       fun = "mean", na.policy = "omit")

ndsm_aircraft_neukirchen8_2_foc <- terra::focal(ndsms_rcl_list[[12]],
                                                w = mov_wind_list[[12]],
                                                fun = "mean", na.policy = "omit")

ndsm_drone_neukirchen9_1_foc <- terra::focal(ndsms_rcl_list[[13]],
                                             w = mov_wind_list[[13]],
                                             fun = "mean", na.policy = "omit")

ndsm_drone_neukirchen9_1_resampled_foc <- terra::focal(ndsms_rcl_list[[14]],
                                                       w = mov_wind_list[[14]],
                                                       fun = "mean", na.policy = "omit")

ndsm_aircraft_neukirchen9_1_foc <- terra::focal(ndsms_rcl_list[[15]],
                                                w = mov_wind_list[[15]],
                                                fun = "mean", na.policy = "omit")

ndsm_drone_neukirchen9_2_foc <- terra::focal(ndsms_rcl_list[[16]],
                                             w = mov_wind_list[[16]],
                                             fun = "mean", na.policy = "omit")

ndsm_drone_neukirchen9_2_resampled_foc <- terra::focal(ndsms_rcl_list[[17]],
                                                       w = mov_wind_list[[17]],
                                                       fun = "mean", na.policy = "omit")

ndsm_aircraft_neukirchen9_2_foc <- terra::focal(ndsms_rcl_list[[18]],
                                                w = mov_wind_list[[18]],
                                                fun = "mean", na.policy = "omit")




terra::writeRaster(ndsm_drone_reinhardshagen_1_foc, paste0("J:/output/canopy_cover/", "ndsm_drone_reinhardshagen_1_foc.tif"))
terra::writeRaster(ndsm_drone_reinhardshagen_1_resampled_foc, paste0("J:/output/canopy_cover/", "ndsm_drone_reinhardshagen_1_resampled_foc.tif"))
terra::writeRaster(ndsm_aircraft_reinhardshagen_1_foc, paste0("J:/output/canopy_cover/", "ndsm_aircraft_reinhardshagen_1_foc.tif"))

terra::writeRaster(ndsm_drone_reinhardshagen_2_foc, paste0("J:/output/canopy_cover/", "ndsm_drone_reinhardshagen_2_foc.tif"))
terra::writeRaster(ndsm_drone_reinhardshagen_2_resampled_foc, paste0("J:/output/canopy_cover/", "ndsm_drone_reinhardshagen_2_resampled_foc.tif"))
terra::writeRaster(ndsm_aircraft_reinhardshagen_2_foc, paste0("J:/output/canopy_cover/", "ndsm_aircraft_reinhardshagen_2_foc.tif"))

terra::writeRaster(ndsm_drone_neukirchen8_1_foc, paste0("J:/output/canopy_cover/", "ndsm_drone_neukirchen8_1_foc.tif"))
terra::writeRaster(ndsm_drone_neukirchen8_1_resampled_foc, paste0("J:/output/canopy_cover/", "ndsm_drone_neukirchen8_1_resampled_foc.tif"))
terra::writeRaster(ndsm_aircraft_neukirchen8_1_foc, paste0("J:/output/canopy_cover/", "ndsm_aircraft_neukirchen8_1_foc.tif"))

terra::writeRaster(ndsm_drone_neukirchen8_2_foc, paste0("J:/output/canopy_cover/", "ndsm_drone_neukirchen8_2_foc.tif"))
terra::writeRaster(ndsm_drone_neukirchen8_2_resampled_foc, paste0("J:/output/canopy_cover/", "ndsm_drone_neukirchen8_2_resampled_foc.tif"))
terra::writeRaster(ndsm_aircraft_neukirchen8_2_foc, paste0("J:/output/canopy_cover/", "ndsm_aircraft_neukirchen8_2_foc.tif"))

terra::writeRaster(ndsm_drone_neukirchen9_1_foc, paste0("J:/output/canopy_cover/", "ndsm_drone_neukirchen9_1_foc.tif"))
terra::writeRaster(ndsm_drone_neukirchen9_1_resampled_foc, paste0("J:/output/canopy_cover/", "ndsm_drone_neukirchen9_1_resampled_foc.tif"))
terra::writeRaster(ndsm_aircraft_neukirchen9_1_foc, paste0("J:/output/canopy_cover/", "ndsm_aircraft_neukirchen9_1_foc.tif"))

terra::writeRaster(ndsm_drone_neukirchen9_2_foc, paste0("J:/output/canopy_cover/", "ndsm_drone_neukirchen9_2_foc.tif"))
terra::writeRaster(ndsm_drone_neukirchen9_2_resampled_foc, paste0("J:/output/canopy_cover/", "ndsm_drone_neukirchen9_2_resampled_foc.tif"))
terra::writeRaster(ndsm_aircraft_neukirchen9_2_foc, paste0("J:/output/canopy_cover/", "ndsm_aircraft_neukirchen9_2_foc.tif"))



##############

# Derive areas of open and dense forest
#--------------------------------------
# Reclassify canopy cover maps:
# all pixels with a value above 0.6 (60 %) are defined as "dense forest",
# all pixels with a value below 0.6 (60 %) are defined as "open forest",
# pixel value >= 0.6 --> 1 (dense forest), pixel value < 0.6 --> 0 (open forest)

matr <- c(0, 0.6, 0,
          0.6, 1, 1)

rclmatr <- matrix(matr, ncol = 3, byrow = TRUE)

canopy_cover_list <- list(ndsm_drone_reinhardshagen_1_resampled_foc, ndsm_aircraft_reinhardshagen_1_foc)

canopy_cover_rcl_list <- lapply(canopy_cover_list, FUN = function(x) terra::classify(x, rclmatr, right = FALSE))

# Group regions:
# adjacent pixels with same values in a neighborhood of 8 cells are identified  
patches <- terra::patches(canopy_cover_rcl_list[[1]], directions = 8,
                          zeroAsNA = TRUE, allowGaps = FALSE)

# Remove small patches (patches < 0.5 ha)
rz <- terra::zonal(terra::cellSize(patches, unit = "ha"), patches, fun = "sum", as.raster = TRUE)
s <- terra::ifel(rz < 0.5, NA, patches)

test <- terra::ifel(s >=1, 1, terra::ifel(is.na(s), 0, s))

par(mfrow = c(1,3))
terra::plot(patches)
terra::plot(test)
terra::plot(canopy_cover_rcl_list[[1]])

###
# Test gap detection
matr <- c(0, 1,
          1, 3)

rclmatr <- matrix(matr, ncol = 2, byrow = TRUE)

canopy_cover_rcl_new_list <- lapply(canopy_cover_rcl_list, FUN = function(x) terra::classify(x, rclmatr))

map_four_classes <- canopy_cover_rcl_new_list[[2]] - ndsms_rcl_list[[3]]

matr <- c(0, 0,
          1, 0,
          2, 1,
          3, 2)

rclmatr <- matrix(matr, ncol = 2, byrow = TRUE)

map_four_classes_rcl <- terra::classify(map_four_classes, rclmatr)

patches <- terra::patches(map_four_classes_rcl, directions = 8,
                          zeroAsNA = TRUE, allowGaps = FALSE)

rz <- terra::zonal(terra::cellSize(patches, unit = "m"), patches, fun = "sum", as.raster = TRUE)
s <- terra::ifel(rz < 10, NA, patches)

# Testing with package "ForestGapR"
library(ForestGapR)

raster_drone <- raster::raster(ndsms_list[[2]])
raster_aircraft <- raster::raster(ndsms_list[[3]])

gaps_drone <- getForestGaps(raster_drone, threshold = 3, size = c(10, 20^4))
gaps_aircraft <- getForestGaps(raster_aircraft, threshold = 3, size = c(10, 20^4))

par_org <- par()
par(mfrow = c(1,2))
raster::plot(raster_drone, col = viridis::viridis(10))
raster::plot(gaps_drone, col = "red", add = TRUE, legend = FALSE)
raster::plot(raster_aircraft, col = viridis::viridis(10))
raster::plot(gaps_aircraft, col = "red", add = TRUE, legend = FALSE)

gaps_stat_drone <- GapStats(gaps_drone, raster_drone)

###


par_org <- par()
par(mfrow = c(1,2))
terra::plot(s)
terra::plot(test)
par(par_org)


par_org <- par()
par(mfrow = c(1,2))
terra::plot(patches, col = grDevices::hcl.colors(10, palette = "Set2"))
terra::plot(canopy_cover_rcl_list[[2]])


par_org <- par()
par(mfow = c(1,2))
terra::plot(s)
terra::plot(canopy_cover_rcl_list[[1]])
par(par_org)


par_org <- par()
par(mfrow = c(1,3))
terra::plot(ndsm_drone_reinhardshagen_1_rcl)
terra::plot(ndsm_drone_reinhardshagen_1_resampled_rcl)
terra::plot(ndsm_aircraft_reinhardshagen_1_rcl)
par(par_org)



par(mfrow = c(1,2))
terra::plot(canopy_cover_rcl_list[[1]])
terra::plot(canopy_cover_rcl_list[[2]])



par_org <- par()
par(mfrow = c(1,2))
#terra::plot(ndsm_drone_reinhardshagen_1_foc)
terra::plot(ndsm_drone_reinhardshagen_1_resampled_foc)
terra::plot(ndsm_aircraft_reinhardshagen_1_foc)
par(par_org)

par(mfrow = c(1,1))
terra::plot(test)







