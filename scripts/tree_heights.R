#-------------------------------------------------------------------------------
# Name:         tree_heights.R
# Author:       Florian Franz
# Description:  script derives tree heights using nDSMs from three different locations:
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



# Quick looks
#-------------
# Parallel plots of all aircraft and drone nDSMs (raw and resampled)
par_org <- par()
par(mfrow = c(6,3))
terra::plot(ndsm_drone_reinhardshagen_1)
terra::plot(ndsm_drone_reinhardshagen_1_resampled)
terra::plot(ndsm_aircraft_reinhardshagen_1)
terra::plot(ndsm_drone_reinhardshagen_2)
terra::plot(ndsm_drone_reinhardshagen_2_resampled)
terra::plot(ndsm_aircraft_reinhardshagen_2)
terra::plot(ndsm_drone_neukirchen8_1)
terra::plot(ndsm_drone_neukirchen8_1_resampled)
terra::plot(ndsm_aircraft_neukirchen8_1)
terra::plot(ndsm_drone_neukirchen8_2)
terra::plot(ndsm_drone_neukirchen8_2_resampled)
terra::plot(ndsm_aircraft_neukirchen8_2)
terra::plot(ndsm_drone_neukirchen9_1)
terra::plot(ndsm_drone_neukirchen9_1_resampled)
terra::plot(ndsm_aircraft_neukirchen9_1)
terra::plot(ndsm_drone_neukirchen9_2)
terra::plot(ndsm_drone_neukirchen9_2_resampled)
terra::plot(ndsm_aircraft_neukirchen9_2)
par(par_org)

# Parallel plots of aircraft and drone nDSMs (raw and resampled)
# for single locations
par_org <- par()
par(mfrow = c(1,3))

terra::plot(ndsm_drone_reinhardshagen_1)
terra::plot(ndsm_drone_reinhardshagen_1_resampled)
terra::plot(ndsm_aircraft_reinhardshagen_1)

terra::plot(ndsm_drone_reinhardshagen_2)
terra::plot(ndsm_drone_reinhardshagen_2_resampled)
terra::plot(ndsm_aircraft_reinhardshagen_2)

terra::plot(ndsm_drone_neukirchen8_1)
terra::plot(ndsm_drone_neukirchen8_1_resampled)
terra::plot(ndsm_aircraft_neukirchen8_1)

terra::plot(ndsm_drone_neukirchen8_2)
terra::plot(ndsm_drone_neukirchen8_2_resampled)
terra::plot(ndsm_aircraft_neukirchen8_2)

terra::plot(ndsm_drone_neukirchen9_1)
terra::plot(ndsm_drone_neukirchen9_1_resampled)
terra::plot(ndsm_aircraft_neukirchen9_1)

terra::plot(ndsm_drone_neukirchen9_2)
terra::plot(ndsm_drone_neukirchen9_2_resampled)
terra::plot(ndsm_aircraft_neukirchen9_2)

par(par_org)

# Single plots
par_org <- par()
terra::plot(ndsm_aircraft_neukirchen8_1)
par(par_org)


# Basic statistics
#------------------
terra::global(ndsm_drone_neukirchen8_1_resampled,
              fun = "mean", na.rm = TRUE)

terra::global(ndsm_aircraft_neukirchen8_1,
              fun = "mean", na.rm = TRUE)




test <- ndsm_aircraft_neukirchen8_1 - ndsm_drone_neukirchen8_1_resampled


terra::plot(test)

