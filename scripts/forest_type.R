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
# All pixels with a value above 3 m are defined as "forest",
# pixel value >= 3 --> 1 (forest), pixel value < 3 --> 0 (non-forest)

matr <- c(0, 3, 0,
          3, 100, 1)

rclmatr <- matrix(matr, ncol = 3, byrow = TRUE)

ndsm_drone_reinhardshagen_1_rcl <- terra::classify(values(ndsm_drone_reinhardshagen_1, na.rm = TRUE),
                                                   rclmatr, right = FALSE)

ndsm_drone_reinhardshagen_1_rcl <- terra::classify(ndsm_drone_reinhardshagen_1,
                                                   rclmatr, right = FALSE)

ndsm_drone_reinhardshagen_1_resampled_rcl <- terra::classify(ndsm_drone_reinhardshagen_1_resampled,
                                                             rclmatr, right = FALSE)

ndsm_aircraft_reinhardshagen_1_rcl <- terra::classify(ndsm_aircraft_reinhardshagen_1,
                                                      rclmatr, right = FALSE)


par_org <- par()
par(mfrow = c(1,3))
terra::plot(ndsm_drone_reinhardshagen_1_rcl)
terra::plot(ndsm_drone_reinhardshagen_1_resampled_rcl)
terra::plot(ndsm_aircraft_reinhardshagen_1_rcl)
par(par_org)










