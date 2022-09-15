#-------------------------------------------------------------------------------
# Name:         forest_type.R
# Author:       Florian Franz
# Description:  script calculates canopy cover and derives areas
#               of open, dense forest and gaps for three different locations:
#               Reinhardshagen, Neukirchen_8 and Neukirchen_9,
#               workflow is similar to the method used in the 'F3 project'
#               (by FVA and NW-FVA) to derive forest structural parameters
#               (https://www.waldwissen.net/de/technik-und-planung/waldinventur/ableitung-von-waldtyp)
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

# Calculate focal mean:
# Circular moving window with radius = 25 m
ndsms_list <- list(ndsm_drone_reinhardshagen_1, ndsm_drone_reinhardshagen_1_resampled, ndsm_aircraft_reinhardshagen_1,
                   ndsm_drone_reinhardshagen_2, ndsm_drone_reinhardshagen_2_resampled, ndsm_aircraft_reinhardshagen_2,
                   ndsm_drone_neukirchen8_1, ndsm_drone_neukirchen8_1_resampled, ndsm_aircraft_neukirchen8_1,
                   ndsm_drone_neukirchen8_2, ndsm_drone_neukirchen8_2_resampled, ndsm_aircraft_neukirchen8_2,
                   ndsm_drone_neukirchen9_1, ndsm_drone_neukirchen9_1_resampled, ndsm_aircraft_neukirchen9_1,
                   ndsm_drone_neukirchen9_2, ndsm_drone_neukirchen9_2_resampled, ndsm_aircraft_neukirchen9_2)


ndsms_rcl_list <- lapply(ndsms_list, FUN = function(x) terra::classify(x, rclmatr, right = FALSE))

mov_wind_list <- lapply(ndsms_rcl_list, FUN = function(x) terra::focalMat(x, d = 25, type = "circle"))


# For testing with a smaller list
#ndsms_rcl_list_test <- ndsms_rcl_list[c(17,18)]


# Try (it has to be iterated over mov_wind_list for argument w)
#ndsms_focal_mean_list <- lapply(ndsms_rcl_list_test, FUN = function(x) terra::focal(x, w = , fun = "mean", na.policy = "omit"))

# Define output path
out_path <- "J:/output/canopy_cover/"

ndsm_drone_reinhardshagen_1_foc_mean <- terra::focal(ndsms_rcl_list[[1]],
                                                     w = mov_wind_list[[1]],
                                                     fun = "mean", na.policy = "omit",
                                                     filename = paste0(out_path, "drone_reinhardshagen_1.tif"))

ndsm_drone_reinhardshagen_1_resampled_foc_mean <- terra::focal(ndsms_rcl_list[[2]],
                                                               w = mov_wind_list[[2]],
                                                               fun = "mean", na.policy = "omit",
                                                               filename = paste0(out_path, "drone_reinhardshagen_1_resampled.tif"))

ndsm_aircraft_reinhardshagen_1_foc_mean <- terra::focal(ndsms_rcl_list[[3]],
                                                        w = mov_wind_list[[3]],
                                                        fun = "mean", na.policy = "omit",
                                                        filename = paste0(out_path, "aircraft_reinhardshagen_1.tif"))

ndsm_drone_reinhardshagen_2_foc_mean <- terra::focal(ndsms_rcl_list[[4]],
                                                     w = mov_wind_list[[4]],
                                                     fun = "mean", na.policy = "omit",
                                                     filename = paste0(out_path, "drone_reinhardshagen_2.tif"))

ndsm_drone_reinhardshagen_2_resampled_foc_mean <- terra::focal(ndsms_rcl_list[[5]],
                                                               w = mov_wind_list[[5]],
                                                               fun = "mean", na.policy = "omit",
                                                               filename = paste0(out_path, "drone_reinhardshagen_2_resampled.tif"))

ndsm_aircraft_reinhardshagen_2_foc_mean <- terra::focal(ndsms_rcl_list[[6]],
                                                        w = mov_wind_list[[6]],
                                                        fun = "mean", na.policy = "omit",
                                                        filename = paste0(out_path, "aircraft_reinhardshagen_2.tif"))

ndsm_drone_neukirchen8_1_foc_mean <- terra::focal(ndsms_rcl_list[[7]],
                                                  w = mov_wind_list[[7]],
                                                  fun = "mean", na.policy = "omit",
                                                  filename = paste0(out_path, "drone_neukirchen8_1.tif"))

ndsm_drone_neukirchen8_1_resampled_foc_mean <- terra::focal(ndsms_rcl_list[[8]],
                                                            w = mov_wind_list[[8]],
                                                            fun = "mean", na.policy = "omit",
                                                            filename = paste0(out_path, "drone_neukirchen8_1_resampled.tif"))

ndsm_aircraft_neukirchen8_1_foc_mean <- terra::focal(ndsms_rcl_list[[9]],
                                                     w = mov_wind_list[[9]],
                                                     fun = "mean", na.policy = "omit",
                                                     filename = paste0(out_path, "aircraft_neukirchen8_1.tif"))

ndsm_drone_neukirchen8_2_foc_mean <- terra::focal(ndsms_rcl_list[[10]],
                                                  w = mov_wind_list[[10]],
                                                  fun = "mean", na.policy = "omit",
                                                  filename = paste0(out_path, "drone_neukirchen8_2.tif"))

ndsm_drone_neukirchen8_2_resampled_foc_mean <- terra::focal(ndsms_rcl_list[[11]],
                                                            w = mov_wind_list[[11]],
                                                            fun = "mean", na.policy = "omit",
                                                            filename = paste0(out_path, "drone_neukirchen8_2_resampled.tif"))

ndsm_aircraft_neukirchen8_2_foc_mean <- terra::focal(ndsms_rcl_list[[12]],
                                                     w = mov_wind_list[[12]],
                                                     fun = "mean", na.policy = "omit",
                                                     filename = paste0(out_path, "aircraft_neukirchen8_2.tif"))

ndsm_drone_neukirchen9_1_foc_mean <- terra::focal(ndsms_rcl_list[[13]],
                                                  w = mov_wind_list[[13]],
                                                  fun = "mean", na.policy = "omit",
                                                  filename = paste0(out_path, "drone_neukirchen9_1.tif"))

ndsm_drone_neukirchen9_1_resampled_foc_mean <- terra::focal(ndsms_rcl_list[[14]],
                                                            w = mov_wind_list[[14]],
                                                            fun = "mean", na.policy = "omit",
                                                            filename = paste0(out_path, "drone_neukirchen9_1_resampled.tif"))

ndsm_aircraft_neukirchen9_1_foc_mean <- terra::focal(ndsms_rcl_list[[15]],
                                                     w = mov_wind_list[[15]],
                                                     fun = "mean", na.policy = "omit",
                                                     filename = paste0(out_path, "aircraft_neukirchen9_1.tif"))

ndsm_drone_neukirchen9_2_foc_mean <- terra::focal(ndsms_rcl_list[[16]],
                                                  w = mov_wind_list[[16]],
                                                  fun = "mean", na.policy = "omit",
                                                  filename = paste0(out_path, "drone_neukirchen9_2.tif"))

ndsm_drone_neukirchen9_2_resampled_foc_mean <- terra::focal(ndsms_rcl_list[[17]],
                                                            w = mov_wind_list[[17]],
                                                            fun = "mean", na.policy = "omit",
                                                            filename = paste0(out_path, "drone_neukirchen9_2_resampled.tif"))

ndsm_aircraft_neukirchen9_2_foc_mean <- terra::focal(ndsms_rcl_list[[18]],
                                                     w = mov_wind_list[[18]],
                                                     fun = "mean", na.policy = "omit",
                                                     filename = paste0(out_path, "aircraft_neukirchen9_2.tif"))

# Load data
file_path_canopy_cover <- "D:/output/canopy_cover/"

canopy_cover_files <- list.files(file_path_canopy_cover,
                                 pattern = glob2rx("*.tif"),
                                 full.names = TRUE)

canopy_cover_files_list <- lapply(canopy_cover_files, FUN = function(x, i) terra::rast(x[i]))

canopy_cover_aircraft <- canopy_cover_files_list[c(1:6)]

canopy_cover_drone <- canopy_cover_files_list[c(7,9,11,13,15,17)]

# Calculate mean
canopy_cover_means_aircraft <- lapply(canopy_cover_aircraft,
                                      function(x) round(mean(values(x, na.rm = TRUE)), 2))

canopy_cover_means_drone <- lapply(canopy_cover_drone,
                                   function(x) round(mean(values(x, na.rm = TRUE)), 2))

canopy_cover_means_aircraft_df <- do.call(rbind, canopy_cover_means_aircraft)

canopy_cover_means_drone_df <- do.call(rbind, canopy_cover_means_drone)

rownames(canopy_cover_means_aircraft_df) <- c("Neukirchen8_1", "Neukirchen8_2",
                                              "Neukirchen9_1", "Neukirchen9_2",
                                              "Reinhardshagen_1", "Reinhardshagen_2")

rownames(canopy_cover_means_drone_df) <- c("Neukirchen8_1", "Neukirchen8_2",
                                           "Neukirchen9_1", "Neukirchen9_2",
                                           "Reinhardshagen_1", "Reinhardshagen_2")

canopy_cover_means_df <- data.frame(cbind(canopy_cover_means_drone_df,
                                          canopy_cover_means_aircraft_df))

colnames(canopy_cover_means_df) <- c("Überschirmung Drohne", "Überschirmung Flugzeug")

# Test plots
par(mfrow = c(1,2))
terra::plot(canopy_cover_drone[[5]],
            col = grDevices::hcl.colors(50, palette = "Greens", rev = TRUE))
terra::plot(canopy_cover_aircraft[[5]],
            col = grDevices::hcl.colors(50, palette = "Greens", rev = TRUE))


# Derive areas of open and dense forest
#--------------------------------------
# Reclassify canopy cover maps:
# all pixels with a value above 0.6 (60 %) are defined as "dense forest",
# all pixels with a value below 0.6 (60 %) are defined as "open forest",
# pixel value >= 0.6 --> 1 (dense forest), pixel value < 0.6 --> 0 (open forest)

matr <- c(0, 0.6, 0,
          0.6, 1, 1)

rclmatr <- matrix(matr, ncol = 3, byrow = TRUE)

canopy_cover_list <- list(ndsm_drone_reinhardshagen_1_foc_mean, ndsm_drone_reinhardshagen_1_resampled_foc_mean, ndsm_aircraft_reinhardshagen_1_foc_mean,
                          ndsm_drone_reinhardshagen_2_foc_mean, ndsm_drone_reinhardshagen_2_resampled_foc_mean, ndsm_aircraft_reinhardshagen_2_foc_mean,
                          ndsm_drone_neukirchen8_1_foc_mean, ndsm_drone_neukirchen8_1_resampled_foc_mean, ndsm_aircraft_neukirchen8_1_foc_mean,
                          ndsm_drone_neukirchen8_2_foc_mean, ndsm_drone_neukirchen8_2_resampled_foc_mean, ndsm_aircraft_neukirchen8_2_foc_mean,
                          ndsm_drone_neukirchen9_1_foc_mean, ndsm_drone_neukirchen9_1_resampled_foc_mean, ndsm_aircraft_neukirchen9_1_foc_mean,
                          ndsm_drone_neukirchen9_2_foc_mean, ndsm_drone_neukirchen9_2_resampled_foc_mean, ndsm_aircraft_neukirchen9_2_foc_mean)

canopy_cover_rcl_list <- lapply(canopy_cover_list, FUN = function(x) terra::classify(x, rclmatr, right = FALSE))

# Define output path
out_path <- "J:/output/forest_type/"

for (i in seq(canopy_cover_rcl_list)){

  terra::writeRaster(canopy_cover_rcl_list[[i]],
                     filename = paste0(out_path, substr(terra::sources(canopy_cover_list[[i]]), 24, nchar(terra::sources(canopy_cover_list[[i]])))))
  
}

# Test plots
par(mfrow = c(1,3))
terra::plot(canopy_cover_rcl_list[[16]])
terra::plot(canopy_cover_rcl_list[[17]])
terra::plot(canopy_cover_rcl_list[[18]])

### This workflow doesn't work
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



# Gap detection
#--------------------------------------
# Reclassify maps of open and dense forest:
# all pixels with value 0 get value 1,
# all pixels with value 1 get value 3,
# then, subtract the reclassified nDSMs from step 1
# from the new reclassified open and dense forest maps
# to add all pixels with heights < 3 m to the classified dense forest

matr <- c(0, 1,
          1, 3)

rclmatr <- matrix(matr, ncol = 2, byrow = TRUE)

canopy_cover_rcl_new_list <- lapply(canopy_cover_rcl_list, FUN = function(x) terra::classify(x, rclmatr))

map_four_classes_list <- mapply('-', canopy_cover_rcl_new_list, ndsms_rcl_list, SIMPLIFY = FALSE)

# combine classes 0 (open forest) and 1 (gaps into open forest)
# into one class = open forest
# result: map with three classes:
# 1 (open forest), 2 (dense forest), 3 (gaps into dense forest)

matr <- c(0, 1,
          1, 1,
          2, 2,
          3, 3)

rclmatr <- matrix(matr, ncol = 2, byrow = TRUE)

map_three_classes_list <- lapply(map_four_classes_list, FUN = function(x) terra::classify(x, rclmatr))

# Define output path
out_path <- "J:/output/forest_type_with_gaps/"

for (i in seq(map_three_classes_list)){
  
  terra::writeRaster(map_three_classes_list[[i]],
                     filename = paste0(out_path, substr(terra::sources(canopy_cover_list[[i]]), 24, nchar(terra::sources(canopy_cover_list[[i]])))))
  
}

### This workflow doesn't work
patches <- terra::patches(map_four_classes_rcl, directions = 8,
                          zeroAsNA = TRUE, allowGaps = FALSE)

rz <- terra::zonal(terra::cellSize(patches, unit = "m"), patches, fun = "sum", as.raster = TRUE)
s <- terra::ifel(rz < 10, NA, patches)
###

# Load data
file_path_forest_type_with_gaps <- "D:/output/forest_type_with_gaps/"

forest_type_with_gaps_files <- list.files(file_path_forest_type_with_gaps,
                                          pattern = glob2rx("*.tif"),
                                          full.names = TRUE)

forest_type_with_gaps_files_list <- lapply(forest_type_with_gaps_files,
                                           FUN = function(x, i) terra::rast(x[i]))

forest_type_with_gaps_aircraft <- forest_type_with_gaps_files_list[c(1:6)]

forest_type_with_gaps_drone <- forest_type_with_gaps_files_list[c(7,9,11,13,15,17)]


# Test plot
par(mfrow = c(1,2))
terra::plot(forest_type_with_gaps_drone[[5]])
#terra::plot(forest_type_with_gaps_drone_resampled[[1]])
terra::plot(forest_type_with_gaps_aircraft[[5]])

# Test plot
par(mfrow = c(1,2))
#colors <- c("palegoldenrod", "palegreen4", "palegreen2")
colors <- c("burlywood", "forestgreen", "palegreen2")
terra::plot(forest_type_with_gaps_drone[[5]], legend = FALSE,
            col = colors, main = "Drohne")
legend("topright", legend = c("offener Bestand", "geschlossener Bestand", "Lücke"),
       fill = colors, border = FALSE, bty = "n")
terra::plot(forest_type_with_gaps_aircraft[[5]], legend = FALSE,
            col = colors, main = "Flugzeug")
legend("topright", legend = c("offener Bestand", "geschlossener Bestand", "Lücke"),
       fill = colors, border = FALSE, bty = "n")


### This plotting method doesn't work yet
levels(map_three_classes_list[[1]]) <- c("offen", "geschlossen", "luecke")
levels(map_three_classes_list[[1]])
###


# Testing with package "ForestGapR"
library(ForestGapR)

raster_drone <- raster::raster(ndsms_list[[1]])
raster_aircraft <- raster::raster(ndsms_list[[3]])

gaps_drone <- getForestGaps(raster_drone, threshold = 3, size = c(10, 10^4))
gaps_aircraft <- getForestGaps(raster_aircraft, threshold = 3, size = c(10, 10^4))

par_org <- par()
par(mfrow = c(1,2))
raster::plot(raster_drone, col = viridis::viridis(10))
raster::plot(gaps_drone, col = "red", add = TRUE, legend = FALSE)
raster::plot(raster_aircraft, col = viridis::viridis(10))
raster::plot(gaps_aircraft, col = "red", add = TRUE, legend = FALSE)

gaps_stat_drone <- GapStats(gaps_drone, raster_drone)
