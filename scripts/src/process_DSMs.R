#-------------------------------------------------------------------------------
# Name:         process_DSMs.R
# Author:       Florian Franz
# Description:  script prepares DSMs from platform aircraft and drone
#               for further derivation of forest structural parameters,
#               steps: 
#               1) merging aircraft tiles at one location into one
#               2) cropping these merged DSMs to the extent of the respective drone DSMs
#               3) resampling drone DSMs to the structure of the cropped aircraft DSMs
#               4) replacing negative values in drone DSMs by zero
#               5) masking the aircraft DSMs based on NA values in the resampled drone DSMs
#               three different locations: Reinhardshagen, Neukirchen_8 and Neukirchen_9
# Data          DSM tif files
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
dsm_reinhardshagen_tile1 <- terra::rast(paste0(filepath_reinhardshagen_dsm_aircraft,
                                                "dsm_325415698.tif"))

dsm_reinhardshagen_tile2 <- terra::rast(paste0(filepath_reinhardshagen_dsm_aircraft,
                                                "dsm_325425698.tif"))

dsm_reinhardshagen_tile3 <- terra::rast(paste0(filepath_reinhardshagen_dsm_aircraft,
                                                "dsm_325415699.tif"))

dsm_neukirchen8_tile1 <- terra::rast(paste0(filepath_neukirchen8_dsm_aircraft,
                                             "dsm_325255644.tif"))

dsm_neukirchen8_tile2 <- terra::rast(paste0(filepath_neukirchen8_dsm_aircraft,
                                             "dsm_325255645.tif"))

dsm_neukirchen8_tile3 <- terra::rast(paste0(filepath_neukirchen8_dsm_aircraft,
                                             "dsm_325275644.tif"))

dsm_neukirchen8_tile4 <- terra::rast(paste0(filepath_neukirchen8_dsm_aircraft,
                                             "dsm_325275645.tif"))

dsm_neukirchen9_tile1 <- terra::rast(paste0(filepath_neukirchen9_dsm_aircraft,
                                             "dsm_325225639.tif"))

dsm_neukirchen9_tile2 <- terra::rast(paste0(filepath_neukirchen9_dsm_aircraft,
                                             "dsm_325235639.tif"))

# Drone #
dsm_drone_reinhardshagen_1 <- terra::rast(paste0(filepath_reinhardshagen_dsm_drone,
                                                  "Reinhardshagen_550_A_2_bDOM10_neu.tif"))

dsm_drone_reinhardshagen_2 <- terra::rast(paste0(filepath_reinhardshagen_dsm_drone,
                                                  "Reinhardshagen_560_A_4_bDOM10_neu.tif"))

dsm_drone_neukirchen8_1 <- terra::rast(paste0(filepath_neukirchen8_dsm_drone,
                                               "Neukirchen_1213_C_1_bDOM10_neu.tif"))

dsm_drone_neukirchen8_2 <- terra::rast(paste0(filepath_neukirchen8_dsm_drone,
                                               "Neukirchen_1071_B_1_bDOM10_neu.tif"))

dsm_drone_neukirchen9_1 <- terra::rast(paste0(filepath_neukirchen9_dsm_drone,
                                               "Neukirchen_521_C_1_bDOM10_neu.tif"))

dsm_drone_neukirchen9_2 <- terra::rast(paste0(filepath_neukirchen9_dsm_drone,
                                               "Neukirchen_525_B_1_bDOM10_neu.tif"))



# Data preprocessing
#--------------------
# Aircraft #

# Merge tiles (DSMs) at one location into one
dsm_aircraft_reinhardshagen_merged <- terra::merge(dsm_reinhardshagen_tile1,
                                                   dsm_reinhardshagen_tile2,
                                                   dsm_reinhardshagen_tile3)

dsm_aircraft_neukirchen8_1_merged <- terra::merge(dsm_neukirchen8_tile1,
                                                  dsm_neukirchen8_tile2)

dsm_aircraft_neukirchen8_2_merged <- terra::merge(dsm_neukirchen8_tile3,
                                                  dsm_neukirchen8_tile4)

dsm_aircraft_neukirchen9_merged <- terra::merge(dsm_neukirchen9_tile1,
                                                dsm_neukirchen9_tile2)

# Crop out subsets of the merged aircraft DSMs
# corresponding to the extent of the respective drone DSMs

# Compare CRS
raster::compareCRS(dsm_aircraft_reinhardshagen_merged, dsm_drone_reinhardshagen_1)
raster::compareCRS(dsm_aircraft_reinhardshagen_merged, dsm_drone_reinhardshagen_2)

raster::compareCRS(dsm_aircraft_neukirchen8_1_merged, dsm_drone_neukirchen8_1)
raster::compareCRS(dsm_aircraft_neukirchen8_2_merged, dsm_drone_neukirchen8_2)

raster::compareCRS(dsm_aircraft_neukirchen9_merged, dsm_drone_neukirchen9_1)
raster::compareCRS(dsm_aircraft_neukirchen9_merged, dsm_drone_neukirchen9_2)

# Adjust drone extents (delete the "32" before x coordinates)
terra::ext(dsm_drone_reinhardshagen_1) <- terra::ext(541261.11176, 542399.4807,
                                                     5698193.6554, 5698869.39952)

terra::ext(dsm_drone_reinhardshagen_2) <- terra::ext(541258.70713, 541461.83012,
                                                     5698996.45757, 5699284.04854)

terra::ext(dsm_drone_neukirchen8_1) <- terra::ext(525125.91429, 525547.74882,
                                                  5644903.61659, 5645566.05673)

terra::ext(dsm_drone_neukirchen8_2) <- terra::ext(527660.89647, 527827.33144,
                                                  5644953.32855, 5645106.86856)

terra::ext(dsm_drone_neukirchen9_1) <- terra::ext(522260.51421, 522545.70185,
                                                  5639230.99173, 5639473.19631)

terra::ext(dsm_drone_neukirchen9_2) <- terra::ext(522746.20346, 523065.67767,
                                                  5639684.50608, 5639935.40729)

# Crop out subsets
dsm_aircraft_reinhardshagen_1_subset <- terra::crop(dsm_aircraft_reinhardshagen_merged,
                                                    dsm_drone_reinhardshagen_1)

dsm_aircraft_reinhardshagen_2_subset <- terra::crop(dsm_aircraft_reinhardshagen_merged,
                                                    dsm_drone_reinhardshagen_2)

dsm_aircraft_neukirchen8_1_subset <- terra::crop(dsm_aircraft_neukirchen8_1_merged,
                                                 dsm_drone_neukirchen8_1)

dsm_aircraft_neukirchen8_2_subset <- terra::crop(dsm_aircraft_neukirchen8_2_merged,
                                                 dsm_drone_neukirchen8_2)

dsm_aircraft_neukirchen9_1_subset <- terra::crop(dsm_aircraft_neukirchen9_merged,
                                                 dsm_drone_neukirchen9_1)

dsm_aircraft_neukirchen9_2_subset <- terra::crop(dsm_aircraft_neukirchen9_merged,
                                                 dsm_drone_neukirchen9_2)

# Calculate focal mean for small band where NA values exist
# based on moving window with size 3
# --> not need for "dsm_aircraft_neukirchen9_1_subset"
dsm_aircraft_reinhardshagen_1_subset <- terra::focal(dsm_aircraft_reinhardshagen_1_subset,
                                                     w = 3, fun = "mean",
                                                     na.policy = "only", na.rm = TRUE)

dsm_aircraft_reinhardshagen_2_subset <- terra::focal(dsm_aircraft_reinhardshagen_2_subset,
                                                     w = 3, fun = "mean",
                                                     na.policy = "only", na.rm = TRUE)

dsm_aircraft_neukirchen8_1_subset <- terra::focal(dsm_aircraft_neukirchen8_1_subset,
                                                  w = 3, fun = "mean",
                                                  na.policy = "only", na.rm = TRUE)

dsm_aircraft_neukirchen8_2_subset <- terra::focal(dsm_aircraft_neukirchen8_2_subset,
                                                  w = 3, fun = "mean",
                                                  na.policy = "only", na.rm = TRUE)

dsm_aircraft_neukirchen9_2_subset <- terra::focal(dsm_aircraft_neukirchen9_2_subset,
                                                  w = 3, fun = "mean",
                                                  na.policy = "only", na.rm = TRUE)

# Resample drone rasters to the structure of the cropped aircraft rasters
# --> necessary for direct comparisons between them (calculations)
dsm_drone_reinhardshagen_1_resampled <- terra::resample(dsm_drone_reinhardshagen_1,
                                                        dsm_aircraft_reinhardshagen_1_subset)

dsm_drone_reinhardshagen_2_resampled <- terra::resample(dsm_drone_reinhardshagen_2,
                                                        dsm_aircraft_reinhardshagen_2_subset)

dsm_drone_neukirchen8_1_resampled <- terra::resample(dsm_drone_neukirchen8_1,
                                                     dsm_aircraft_neukirchen8_1_subset)

dsm_drone_neukirchen8_2_resampled <- terra::resample(dsm_drone_neukirchen8_2,
                                                     dsm_aircraft_neukirchen8_2_subset)

dsm_drone_neukirchen9_1_resampled <- terra::resample(dsm_drone_neukirchen9_1,
                                                     dsm_aircraft_neukirchen9_1_subset)

dsm_drone_neukirchen9_2_resampled <- terra::resample(dsm_drone_neukirchen9_2,
                                                     dsm_aircraft_neukirchen9_2_subset)

# Mask out the aircraft DSMs based on NA values in the resampled drone DSMs
dsm_aircraft_reinhardshagen_1 <- terra::mask(dsm_aircraft_reinhardshagen_1_subset,
                                             dsm_drone_reinhardshagen_1_resampled)

dsm_aircraft_reinhardshagen_2 <- terra::mask(dsm_aircraft_reinhardshagen_2_subset,
                                             dsm_drone_reinhardshagen_2_resampled)

dsm_aircraft_neukirchen8_1 <- terra::mask(dsm_aircraft_neukirchen8_1_subset,
                                          dsm_drone_neukirchen8_1_resampled)

dsm_aircraft_neukirchen8_2 <- terra::mask(dsm_aircraft_neukirchen8_2_subset,
                                          dsm_drone_neukirchen8_2_resampled)

dsm_aircraft_neukirchen9_1 <- terra::mask(dsm_aircraft_neukirchen9_1_subset,
                                          dsm_drone_neukirchen9_1_resampled)

dsm_aircraft_neukirchen9_2 <- terra::mask(dsm_aircraft_neukirchen9_2_subset,
                                          dsm_drone_neukirchen9_2_resampled)