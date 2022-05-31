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
                                             "ndsm_325255645.tif"))

ndsm_neukirchen8_tile2 <- terra::rast(paste0(filepath_neukirchen8_ndsm_aircraft,
                                             "ndsm_325275645.tif"))

ndsm_neukirchen9_tile1 <- terra::rast(paste0(filepath_neukirchen9_ndsm_aircraft,
                                             "ndsm_325225639.tif"))

ndsm_neukirchen9_tile2 <- terra::rast(paste0(filepath_neukirchen9_ndsm_aircraft,
                                             "ndsm_325235639.tif"))

# Drone --> currently still DSM, has to be changed to nDSM!!! #
dsm_drone_reinhardshagen1 <- terra::rast("J:/data/drone/reinhardshagen_e_str1_bis_4/DSM/Reinhardshagen_550_A_2_DEM10.tif")
dsm_drone_reinhardshagen2 <- terra::rast("J:/data/drone/reinhardshagen_e_str1_bis_4/DSM/Reinhardshagen_560_A_4_DEM10.tif")

dsm_drone_neukirchen91 <- terra::rast("J:/data/drone/neukirchen_e_str9/DSM/Neukirchen_521_C_1_DEM10.tif")
dsm_drone_neukirchen92 <- terra::rast("J:/data/drone/neukirchen_e_str9/DSM/Neukirchen_525_B_1_DEM10.tif")


# Data preprocessing
#--------------------
# Aircraft #

# Merge raster files (nDSMs) at one location into one
ndsm_aircraft_reinhardshagen <- terra::merge(ndsm_reinhardshagen_tile1,
                                             ndsm_reinhardshagen_tile2,
                                             ndsm_reinhardshagen_tile3)

ndsm_aircraft_neukirchen9 <- terra::merge(ndsm_neukirchen9_tile1,
                                          ndsm_neukirchen9_tile2)

# Crop out subsets of the aircraft nDSMs
# corresponding to the extent of the respective drone nDSMs

# Compare CRS
raster::compareCRS(ndsm_aircraft_reinhardshagen, dsm_drone_reinhardshagen1)
raster::compareCRS(ndsm_aircraft_reinhardshagen, dsm_drone_reinhardshagen2)

raster::compareCRS(ndsm_aircraft_neukirchen9, dsm_drone_neukirchen91)
raster::compareCRS(ndsm_aircraft_neukirchen9, dsm_drone_neukirchen92)

# Get extent of drone nDSMs
terra::ext(dsm_drone_reinhardshagen1)
terra::ext(dsm_drone_reinhardshagen2)

terra::ext(dsm_drone_neukirchen91)
terra::ext(dsm_drone_neukirchen92)

# Define this extent new (UTM zone-number 32 has to be removed)
terra::ext(dsm_drone_reinhardshagen1) <- terra::ext(541261.1117562, 542399.480701,
                                                    5698193.6553981, 5698869.39952059)

terra::ext(dsm_drone_reinhardshagen2) <- terra::ext(541258.7071286, 541461.8301163,
                                                    5698996.4575736, 5699284.04853613)

terra::ext(dsm_drone_neukirchen91) <- terra::ext(522260.5142122, 522545.7018505,
                                                 5639230.99173351, 5639473.19630677)

terra::ext(dsm_drone_neukirchen92) <- terra::ext(522746.2034559, 523065.677675,
                                                 5639684.5060757, 5639935.40729279)

# Crop out subsets
ndsm_aircraft_reinhardshagen1 <- terra::crop(ndsm_aircraft_reinhardshagen,
                                             dsm_drone_reinhardshagen1)

ndsm_aircraft_reinhardshagen2 <- terra::crop(ndsm_aircraft_reinhardshagen,
                                             dsm_drone_reinhardshagen2)

ndsm_aircraft_neukirchen91 <- terra::crop(ndsm_aircraft_neukirchen9,
                                          dsm_drone_neukirchen91)

ndsm_aircraft_neukirchen92 <- terra::crop(ndsm_aircraft_neukirchen9,
                                          dsm_drone_neukirchen92)




















# Plots
par_org <- par()
par(mfrow = c(1,2))
terra::plot(ndsm_aircraft_neukirchen91)
terra::plot(dsm_drone_neukirchen91)
par(par_org)








