#-------------------------------------------------------------------------------
# Name:         tree_heights.R
# Author:       Florian Franz
# Description:  script derives tree heights using nDSMs from three different locations:
#               Reinhardshagen, Neukirchen_8 and Neukirchen_9
# Data          nDSM tif files from platforms aircraft (1 m resolution) and drone (0.1 m resolution)
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


ndsm_neukirchen9_tile1_rast <- raster::raster(paste0(filepath_neukirchen9_ndsm_aircraft,
                                             "ndsm_325225639.tif"))

ndsm_neukirchen9_tile2_rast <- raster::raster(paste0(filepath_neukirchen9_ndsm_aircraft,
                                             "ndsm_325235639.tif"))

# Drone --> currently still DSM, has to be changed to nDSM!!! #
dsm_drone_reinhardshagen_1 <- terra::rast("J:/Drohne/nwe_flaechen_hartwich/daten/Neukirchen/521_C_1/Agisoft/Neukirchen_521_C_1_DEM10.tif")
dsm_drone_reinhardshagen_2 <- terra::rast("J:/data/drone/reinhardshagen_e_str1_bis_4/DSM/Reinhardshagen_560_A_4_bDOM10_neu.tif")

dsm_drone_neukirchen8_1 <- terra::rast("J:/data/drone/neukirchen_e_str8/DSM/Neukirchen_1071_B_1_bDOM10_neu.tif")
dsm_drone_neukirchen8_2 <- terra::rast("J:/data/drone/neukirchen_e_str8/DSM/Neukirchen_1213_C_1_bDOM10_neu.tif")

dsm_drone_neukirchen9_1 <- terra::rast("J:/Drohne/nwe_flaechen_hartwich/daten/Neukirchen/521_C_1/Agisoft/Neukirchen_521_C_1_DEM10.tif")
dsm_drone_neukirchen9_2 <- terra::rast("J:/data/drone/neukirchen_e_str9/DSM/Neukirchen_525_B_1_bDOM10_neu.tif")


ndsm_drone_neukirchen8_1_1 <- terra::rast(paste0(filepath_neukirchen8_ndsm_drone,
                                                 "ndsm_325255644.tif"))

ndsm_drone_neukirchen8_1_2 <- terra::rast(paste0(filepath_neukirchen8_ndsm_drone,
                                                 "ndsm_325255645.tif"))

ndsm_drone_neukirchen8_2_1 <- terra::rast(paste0(filepath_neukirchen8_ndsm_drone,
                                                 "ndsm_325275644.tif"))

ndsm_drone_neukirchen8_2_2 <- terra::rast(paste0(filepath_neukirchen8_ndsm_drone,
                                                 "ndsm_325275645.tif"))

ndsm_drone_neukirchen9_1 <- terra::rast(paste0(filepath_neukirchen9_ndsm_drone,
                                               "ndsm_325225639.tif"))

ndsm_drone_neukirchen9_1_rast <- raster::raster(paste0(filepath_neukirchen9_ndsm_drone,
                                                       "ndsm_325225639.tif"))

ndsm_drone_neukirchen9_2_1 <- terra::rast(paste0(filepath_neukirchen9_ndsm_drone,
                                               "ndsm_325235639_1.tif"))

ndsm_drone_neukirchen9_2_2 <- terra::rast(paste0(filepath_neukirchen9_ndsm_drone,
                                                 "ndsm_325235639_2.tif"))



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

ndsm_aircraft_neukirchen9_merged_rast <- raster::merge(ndsm_neukirchen9_tile1_rast,
                                                 ndsm_neukirchen9_tile2_rast)

# Drone #

# Merge parts (nDSMs) at one location into one (where its necessary)
ndsm_drone_neukirchen8_1 <- terra::merge(ndsm_drone_neukirchen8_1_1,
                                         ndsm_drone_neukirchen8_1_2)

ndsm_drone_neukirchen8_2 <- terra::merge(ndsm_drone_neukirchen8_2_1,
                                         ndsm_drone_neukirchen8_2_2)

ndsm_drone_neukirchen9_2 <- terra::merge(ndsm_drone_neukirchen9_2_1,
                                         ndsm_drone_neukirchen9_2_2)

# Crop out subsets of the merged aircraft nDSMs
# corresponding to the extent of the respective drone nDSMs

# Compare CRS
raster::compareCRS(ndsm_aircraft_reinhardshagen_merged, dsm_drone_reinhardshagen_1)
raster::compareCRS(ndsm_aircraft_reinhardshagen_merged, dsm_drone_reinhardshagen_2)

raster::compareCRS(ndsm_aircraft_neukirchen8_1_merged, ndsm_drone_neukirchen8_1)
raster::compareCRS(ndsm_aircraft_neukirchen8_2_merged, ndsm_drone_neukirchen8_2)

raster::compareCRS(ndsm_aircraft_neukirchen9_merged, ndsm_drone_neukirchen9_1)
raster::compareCRS(ndsm_aircraft_neukirchen9_merged, ndsm_drone_neukirchen9_2)

# Crop out subsets
ndsm_aircraft_reinhardshagen_1 <- terra::crop(ndsm_aircraft_reinhardshagen_merged,
                                              dsm_drone_reinhardshagen_1)

ndsm_aircraft_reinhardshagen_2 <- terra::crop(ndsm_aircraft_reinhardshagen_merged,
                                              dsm_drone_reinhardshagen_2)

ndsm_aircraft_neukirchen8_1 <- terra::crop(ndsm_aircraft_neukirchen8_1_merged,
                                           ndsm_drone_neukirchen8_1)

ndsm_aircraft_neukirchen8_2 <- terra::crop(ndsm_aircraft_neukirchen8_2_merged,
                                           ndsm_drone_neukirchen8_2)

ndsm_aircraft_neukirchen9_1 <- terra::crop(ndsm_aircraft_neukirchen9_merged,
                                           ndsm_drone_neukirchen9_1)

ndsm_aircraft_neukirchen9_1_rast <- raster::crop(ndsm_aircraft_neukirchen9_merged_rast,
                                           ndsm_drone_neukirchen9_1_rast)

ndsm_aircraft_neukirchen9_2 <- terra::crop(ndsm_aircraft_neukirchen9_merged,
                                           ndsm_drone_neukirchen9_2)

# Resample drone rasters to the structure of the cropped aircraft rasters
# --> necessary for direct comparisons between them (calculations)
ndsm_drone_neukirchen8_1_resampled <- terra::resample(ndsm_drone_neukirchen8_1,
                                                      ndsm_aircraft_neukirchen8_1)

ndsm_drone_neukirchen8_2_resampled <- terra::resample(ndsm_drone_neukirchen8_2,
                                                      ndsm_aircraft_neukirchen8_2)

ndsm_drone_neukirchen9_1_resampled <- terra::resample(ndsm_drone_neukirchen9_1,
                                                      ndsm_aircraft_neukirchen9_1)

ndsm_drone_neukirchen9_2_resampled <- terra::resample(ndsm_drone_neukirchen9_2,
                                                      ndsm_aircraft_neukirchen9_2)


# Replace negative values by zero in drone nDSMs
# (already done for aircraft nDSMs in external script)
ndsm_drone_neukirchen8_1[ndsm_drone_neukirchen8_1 < 0] <- 0
ndsm_drone_neukirchen8_1_resampled[ndsm_drone_neukirchen8_1_resampled < 0] <- 0

ndsm_drone_neukirchen8_2[ndsm_drone_neukirchen8_2 < 0] <- 0
ndsm_drone_neukirchen8_2_resampled[ndsm_drone_neukirchen8_2_resampled < 0] <- 0

ndsm_drone_neukirchen9_1[ndsm_drone_neukirchen9_1 < 0] <- 0
ndsm_drone_neukirchen9_1_resampled[ndsm_drone_neukirchen9_1_resampled < 0] <- 0

ndsm_drone_neukirchen9_2[ndsm_drone_neukirchen9_2 < 0] <- 0
ndsm_drone_neukirchen9_2_resampled[ndsm_drone_neukirchen9_2_resampled < 0] <- 0







# Parallel Plots of aircraft and drone nDSMs
par_org <- par()
par(mfrow = c(1,3))
terra::plot(ndsm_drone_neukirchen8_1)
terra::plot(ndsm_drone_neukirchen8_1_resampled)
terra::plot(ndsm_aircraft_neukirchen8_1)
par(par_org)

par_org <- par()
par(mfrow = c(1,2))
raster::plot(test)
raster::plot(ndsm_aircraft_neukirchen9_1_rast)
par(par_org)


test <- ndsm_aircraft_neukirchen9_1 - resampled_drone_rast

test1 <- ndsm_aircraft_neukirchen9_1_rast - resampled_drone_rast

terra::plot(ndsm_drone_neukirchen8_1)


