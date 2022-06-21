#-------------------------------------------------------------------------------
# Name:         tree_heights.R
# Author:       Florian Franz
# Description:  script derives tree heights using nDSMs from three different locations:
#               Reinhardshagen, Neukirchen_8 and Neukirchen_9
#               basic statistics and difference between aircraft and drone nDSMs are calculated
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
# Mean height
ndsms_list <- list(ndsm_drone_reinhardshagen_1, ndsm_drone_reinhardshagen_1_resampled, ndsm_aircraft_reinhardshagen_1,
                   ndsm_drone_reinhardshagen_2, ndsm_drone_reinhardshagen_2_resampled, ndsm_aircraft_reinhardshagen_2,
                   ndsm_drone_neukirchen8_1, ndsm_drone_neukirchen8_1_resampled, ndsm_aircraft_neukirchen8_1,
                   ndsm_drone_neukirchen8_2, ndsm_drone_neukirchen8_2_resampled, ndsm_aircraft_neukirchen8_2,
                   ndsm_drone_neukirchen9_1, ndsm_drone_neukirchen9_1_resampled, ndsm_aircraft_neukirchen9_1,
                   ndsm_drone_neukirchen9_2, ndsm_drone_neukirchen9_2_resampled, ndsm_aircraft_neukirchen9_2)

# Mean and maximum heights
list_means <- c()
list_maxs <- c()

for (i in ndsms_list){
  
  list_means <- append(list_means, terra::global(i, fun = "mean", na.rm = TRUE))
  list_maxs <- append(list_maxs, terra::global(i, fun = "max", na.rm = TRUE))
  
}

df_means <- do.call(cbind, list_means)
df_maxs <- do.call(cbind, list_maxs)

#test_df <- do.call(cbind, Map(rbind, Name = seq_along(test), test))

#test_df <- dplyr::bind_cols(test)

locations <- c("Reinhardshagen_1", "Reinhardshagen_1", "Reinhardshagen_1",
               "Reinhardshagen_2", "Reinhardshagen_2", "Reinhardshagen_2",
               "Neukirchen8_1", "Neukirchen8_1", "Neukirchen8_1",
               "Neukirchen8_2","Neukirchen8_2", "Neukirchen8_2",
               "Neukirchen9_1", "Neukirchen9_1", "Neukirchen9_1",
               "Neukirchen9_2", "Neukirchen9_2", "Neukirchen9_2")

df_means_split <- split(df_means, locations)
df_maxs_split <- split(df_maxs, locations)

df_mean_heights <- do.call(rbind, df_means_split)
df_max_heights <- do.call(rbind, df_maxs_split)

colnames(df_mean_heights) <- c("Drone", "Drone_resampled", "Aircraft")
colnames(df_max_heights) <- c("Drone", "Drone_resampled", "Aircraft")



# Difference
#------------------
diff_reinhardshagen_1 <- ndsm_aircraft_reinhardshagen_1 - ndsm_drone_reinhardshagen_1_resampled
diff_reinhardshagen_2 <- ndsm_aircraft_reinhardshagen_2 - ndsm_drone_reinhardshagen_2_resampled
diff_neukirchen8_1 <- ndsm_aircraft_neukirchen8_1 - ndsm_drone_neukirchen8_1_resampled
diff_neukirchen8_2 <- ndsm_aircraft_neukirchen8_2 - ndsm_drone_neukirchen8_2_resampled
diff_neukirchen9_1 <- ndsm_aircraft_neukirchen9_1 - ndsm_drone_neukirchen9_1_resampled
diff_neukirchen9_2 <- ndsm_aircraft_neukirchen9_2 - ndsm_drone_neukirchen9_2_resampled

par(par_org)
par(mfrow = c(1,3))

terra::plot(diff_reinhardshagen_1, col = grDevices::hcl.colors(50, palette = "blue-red 3"))
terra::plot(ndsm_drone_reinhardshagen_1_resampled)
terra::plot(ndsm_aircraft_reinhardshagen_1)

terra::plot(diff_reinhardshagen_2, col = grDevices::hcl.colors(50, palette = "blue-red 3"))
terra::plot(ndsm_drone_reinhardshagen_2_resampled)
terra::plot(ndsm_aircraft_reinhardshagen_2)

terra::plot(diff_neukirchen8_1, col = grDevices::hcl.colors(50, palette = "blue-red 3"))
terra::plot(ndsm_drone_neukirchen8_1_resampled)
terra::plot(ndsm_aircraft_neukirchen8_1)

terra::plot(diff_neukirchen8_2, col = grDevices::hcl.colors(50, palette = "blue-red 3"))
terra::plot(ndsm_drone_neukirchen8_2_resampled)
terra::plot(ndsm_aircraft_neukirchen8_2)

terra::plot(diff_neukirchen9_1, col = grDevices::hcl.colors(50, palette = "blue-red 3"))
terra::plot(ndsm_drone_neukirchen9_1_resampled)
terra::plot(ndsm_aircraft_neukirchen9_1)

terra::plot(diff_neukirchen9_2, col = grDevices::hcl.colors(50, palette = "blue-red 3"))
terra::plot(ndsm_drone_neukirchen9_2_resampled)
terra::plot(ndsm_aircraft_neukirchen9_2)

