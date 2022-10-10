#-------------------------------------------------------------------------------
# Name:         process_nDSMs.R
# Author:       Florian Franz
# Description:  script prepares nDSMs from platforms aircraft and drone
#               for further derivation of forest structural parameters,
#               steps: 
#               1) merging aircraft tiles at one location into one
#               2) cropping these merged nDSMs to the extent of the respective drone nDSMs
#               3) resampling drone nDSMs to the structure of the cropped aircraft nDSMs
#               4) replacing negative values in drone nDSMs by zero
#               5) masking the aircraft nDSMs based on NA values in the resampled drone nDSMs
#               three different locations: Reinhardshagen, Neukirchen_8 and Neukirchen_9
# Data          nDSM tif files
#-------------------------------------------------------------------------------

# Load packages
#---------------
library(terra)
library(raster)



# Data import
#-------------
# Source script for filepath definitions
source("J:/scripts/src/set_filepaths.R", local = TRUE)

# Aircraft #
ndsm_reinhardshagen_tile1 <- terra::rast(paste0(filepath_reinhardshagen_ndsm_aircraft,
                                                "ndsm_325415698.tif"))

ndsm_reinhardshagen_tile2 <- terra::rast(paste0(filepath_reinhardshagen_ndsm_aircraft,
                                                "ndsm_325425698.tif"))

ndsm_reinhardshagen_tile3 <- terra::rast(paste0(filepath_reinhardshagen_ndsm_aircraft,
                                                "ndsm_325415699.tif"))

ndsm_neukirchen8_tile1 <- terra::rast(paste0(filepath_neukirchen8_ndsm_aircraft,
                                             "ndsm_325255644.tif"))

ndsm_neukirchen8_tile2 <- terra::rast(paste0(filepath_neukirchen8_ndsm_aircraft,
                                             "ndsm_325255645.tif"))

ndsm_neukirchen8_tile3 <- terra::rast(paste0(filepath_neukirchen8_ndsm_aircraft,
                                             "ndsm_325275644.tif"))

ndsm_neukirchen8_tile4 <- terra::rast(paste0(filepath_neukirchen8_ndsm_aircraft,
                                             "ndsm_325275645.tif"))

ndsm_neukirchen9_tile1 <- terra::rast(paste0(filepath_neukirchen9_ndsm_aircraft,
                                             "ndsm_325225639.tif"))

ndsm_neukirchen9_tile2 <- terra::rast(paste0(filepath_neukirchen9_ndsm_aircraft,
                                             "ndsm_325235639.tif"))

# Drone #
ndsm_drone_reinhardshagen_1 <- terra::rast(paste0(filepath_reinhardshagen_ndsm_drone,
                                                  "ndsm_32541_25698.tif"))

ndsm_drone_reinhardshagen_2 <- terra::rast(paste0(filepath_reinhardshagen_ndsm_drone,
                                                  "ndsm_325415698_9.tif"))

ndsm_drone_neukirchen8_1 <- terra::rast(paste0(filepath_neukirchen8_ndsm_drone,
                                               "ndsm_325255644_5.tif"))

ndsm_drone_neukirchen8_2 <- terra::rast(paste0(filepath_neukirchen8_ndsm_drone,
                                               "ndsm_325275644_5.tif"))

ndsm_drone_neukirchen9_1 <- terra::rast(paste0(filepath_neukirchen9_ndsm_drone,
                                               "ndsm_325225639.tif"))

ndsm_drone_neukirchen9_2 <- terra::rast(paste0(filepath_neukirchen9_ndsm_drone,
                                               "ndsm_325235639.tif"))



# Data preprocessing
#--------------------
# Aircraft #

# Merge tiles (nDSMs) at one location into one
ndsm_aircraft_reinhardshagen_merged <- terra::merge(ndsm_reinhardshagen_tile1,
                                                    ndsm_reinhardshagen_tile2,
                                                    ndsm_reinhardshagen_tile3)

ndsm_aircraft_neukirchen8_1_merged <- terra::merge(ndsm_neukirchen8_tile1,
                                                   ndsm_neukirchen8_tile2)

ndsm_aircraft_neukirchen8_2_merged <- terra::merge(ndsm_neukirchen8_tile3,
                                                   ndsm_neukirchen8_tile4)

ndsm_aircraft_neukirchen9_merged <- terra::merge(ndsm_neukirchen9_tile1,
                                                 ndsm_neukirchen9_tile2)

# Crop out subsets of the merged aircraft nDSMs
# corresponding to the extent of the respective drone nDSMs

# Compare CRS
raster::compareCRS(ndsm_aircraft_reinhardshagen_merged, ndsm_drone_reinhardshagen_1)
raster::compareCRS(ndsm_aircraft_reinhardshagen_merged, ndsm_drone_reinhardshagen_2)

raster::compareCRS(ndsm_aircraft_neukirchen8_1_merged, ndsm_drone_neukirchen8_1)
raster::compareCRS(ndsm_aircraft_neukirchen8_2_merged, ndsm_drone_neukirchen8_2)

raster::compareCRS(ndsm_aircraft_neukirchen9_merged, ndsm_drone_neukirchen9_1)
raster::compareCRS(ndsm_aircraft_neukirchen9_merged, ndsm_drone_neukirchen9_2)

# Crop out subsets
ndsm_aircraft_reinhardshagen_1_subset <- terra::crop(ndsm_aircraft_reinhardshagen_merged,
                                                     ndsm_drone_reinhardshagen_1)

ndsm_aircraft_reinhardshagen_2_subset <- terra::crop(ndsm_aircraft_reinhardshagen_merged,
                                                     ndsm_drone_reinhardshagen_2)

ndsm_aircraft_neukirchen8_1_subset <- terra::crop(ndsm_aircraft_neukirchen8_1_merged,
                                                  ndsm_drone_neukirchen8_1)

ndsm_aircraft_neukirchen8_2_subset <- terra::crop(ndsm_aircraft_neukirchen8_2_merged,
                                                  ndsm_drone_neukirchen8_2)

ndsm_aircraft_neukirchen9_1_subset <- terra::crop(ndsm_aircraft_neukirchen9_merged,
                                                  ndsm_drone_neukirchen9_1)

ndsm_aircraft_neukirchen9_2_subset <- terra::crop(ndsm_aircraft_neukirchen9_merged,
                                                  ndsm_drone_neukirchen9_2)

# Resample drone rasters to the structure of the cropped aircraft rasters
# --> necessary for direct comparisons between them (calculations)
ndsm_drone_reinhardshagen_1_resampled <- terra::resample(ndsm_drone_reinhardshagen_1,
                                                         ndsm_aircraft_reinhardshagen_1_subset)

ndsm_drone_reinhardshagen_2_resampled <- terra::resample(ndsm_drone_reinhardshagen_2,
                                                         ndsm_aircraft_reinhardshagen_2_subset)

ndsm_drone_neukirchen8_1_resampled <- terra::resample(ndsm_drone_neukirchen8_1,
                                                      ndsm_aircraft_neukirchen8_1_subset)

ndsm_drone_neukirchen8_2_resampled <- terra::resample(ndsm_drone_neukirchen8_2,
                                                      ndsm_aircraft_neukirchen8_2_subset)

ndsm_drone_neukirchen9_1_resampled <- terra::resample(ndsm_drone_neukirchen9_1,
                                                      ndsm_aircraft_neukirchen9_1_subset)

ndsm_drone_neukirchen9_2_resampled <- terra::resample(ndsm_drone_neukirchen9_2,
                                                      ndsm_aircraft_neukirchen9_2_subset)

# Replace negative values by zero in drone nDSMs
# (already done for aircraft nDSMs in external script)
ndsm_drone_reinhardshagen_1[ndsm_drone_reinhardshagen_1 < 0] <- 0
ndsm_drone_reinhardshagen_1_resampled[ndsm_drone_reinhardshagen_1_resampled < 0] <- 0

ndsm_drone_reinhardshagen_2[ndsm_drone_reinhardshagen_2 < 0] <- 0
ndsm_drone_reinhardshagen_2_resampled[ndsm_drone_reinhardshagen_2_resampled < 0] <- 0

ndsm_drone_neukirchen8_1[ndsm_drone_neukirchen8_1 < 0] <- 0
ndsm_drone_neukirchen8_1_resampled[ndsm_drone_neukirchen8_1_resampled < 0] <- 0

ndsm_drone_neukirchen8_2[ndsm_drone_neukirchen8_2 < 0] <- 0
ndsm_drone_neukirchen8_2_resampled[ndsm_drone_neukirchen8_2_resampled < 0] <- 0

ndsm_drone_neukirchen9_1[ndsm_drone_neukirchen9_1 < 0] <- 0
ndsm_drone_neukirchen9_1_resampled[ndsm_drone_neukirchen9_1_resampled < 0] <- 0

ndsm_drone_neukirchen9_2[ndsm_drone_neukirchen9_2 < 0] <- 0
ndsm_drone_neukirchen9_2_resampled[ndsm_drone_neukirchen9_2_resampled < 0] <- 0

# Mask out the aircraft nDSMs based on NA values in the resampled drone nDSMs
ndsm_aircraft_reinhardshagen_1 <- terra::mask(ndsm_aircraft_reinhardshagen_1_subset,
                                              ndsm_drone_reinhardshagen_1_resampled)

ndsm_aircraft_reinhardshagen_2 <- terra::mask(ndsm_aircraft_reinhardshagen_2_subset,
                                              ndsm_drone_reinhardshagen_2_resampled)

ndsm_aircraft_neukirchen8_1 <- terra::mask(ndsm_aircraft_neukirchen8_1_subset,
                                           ndsm_drone_neukirchen8_1_resampled)

ndsm_aircraft_neukirchen8_2 <- terra::mask(ndsm_aircraft_neukirchen8_2_subset,
                                           ndsm_drone_neukirchen8_2_resampled)

ndsm_aircraft_neukirchen9_1 <- terra::mask(ndsm_aircraft_neukirchen9_1_subset,
                                           ndsm_drone_neukirchen9_1_resampled)

ndsm_aircraft_neukirchen9_2 <- terra::mask(ndsm_aircraft_neukirchen9_2_subset,
                                           ndsm_drone_neukirchen9_2_resampled)