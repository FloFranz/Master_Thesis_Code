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
library(lidR)



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



# Write to disk
#---------------
# Define ouput path
out_path <- "J:/output/tree_heights/"


ndsms_list <- list(ndsm_aircraft_reinhardshagen_1, ndsm_aircraft_reinhardshagen_2,
                   ndsm_aircraft_neukirchen8_1, ndsm_aircraft_neukirchen8_2,
                   ndsm_aircraft_neukirchen9_1, ndsm_aircraft_neukirchen9_2,
                   ndsm_drone_reinhardshagen_1, ndsm_drone_reinhardshagen_2,
                   ndsm_drone_neukirchen8_1, ndsm_drone_neukirchen8_2,
                   ndsm_drone_neukirchen9_1, ndsm_drone_neukirchen9_2,
                   ndsm_drone_reinhardshagen_1_resampled, ndsm_drone_reinhardshagen_2_resampled,
                   ndsm_drone_neukirchen8_1_resampled, ndsm_drone_neukirchen8_2_resampled,
                   ndsm_drone_neukirchen9_1_resampled, ndsm_drone_neukirchen9_2_resampled)

for (file in seq(ndsms_list)){
    
  terra::writeRaster(ndsms_list[[file]],
                     filename = paste0(out_path, 
                                       substr(names(ndsms_list[[file]]), 6,
                                              nchar(names(ndsms_list[[file]]))), ".tif"),
                     overwrite = TRUE)
  
}


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



# Comparing aircraft and drone nDSMs (tree heights)
# with LiDAR derived CHMs
#---------------------------------------------------
# Read data
file_path_ndsms <- "J:/output/tree_heights/"
file_path_chms <- "J:/output/lidar_CHMs/"

ndsms_aircraft_drone_files <- list.files(file_path_ndsms,
                                         pattern = glob2rx("*.tif"),
                                         full.names = TRUE)

chm_lidar_files <- list.files(file_path_chms,
                              pattern = glob2rx("*.tif"),
                              full.names = TRUE)

# Read all aircraft, drone (raw and resampled) and LiDAR files into one list
chm_lidar_files_list <- lapply(chm_lidar_files, FUN = function(x, i) terra::rast(x[i]))

ndsms_aircraft_drone_files_list <- lapply(ndsms_aircraft_drone_files, FUN = function(x, i) terra::rast(x[i]))

ndsms_aircraft_files_list <- ndsms_aircraft_drone_files_list[c(1:6)]

ndsms_drone_files_list <- ndsms_aircraft_drone_files_list[c(7,9,11,13,15,17)]

ndsms_drone_resampled_files_list <- ndsms_aircraft_drone_files_list[c(8,10,12,14,16,18)]


# Plots
par_org <- par()
par(mfrow = c(1,3))
terra::plot(chm_lidar_files_list[[1]])              # LiDAR CHM
#terra::plot(ndsms_drone_files_list[[5]])            # Drone
terra::plot(ndsms_drone_resampled_files_list[[1]])  # Drone resampled
terra::plot(ndsms_aircraft_files_list[[1]])         # Aircraft
par(par_org)

# Calculate differences
diff_lidar_aircraft <- mapply('-', chm_lidar_files_list, ndsms_aircraft_files_list,
                              SIMPLIFY = FALSE)

diff_mean_lidar_aircraft <- lapply(diff_lidar_aircraft,
                                   function(x) mean(values(x, na.rm = TRUE)))

diff_lidar_drone <- mapply('-', chm_lidar_files_list, ndsms_drone_resampled_files_list,
                           SIMPLIFY = FALSE)

diff_mean_lidar_drone <- lapply(diff_lidar_drone,
                                function(x) mean(values(x, na.rm = TRUE)))

# Create data frame with differences results
diff_mean_lidar_aircraft_df <- do.call(rbind, diff_mean_lidar_aircraft)

diff_mean_lidar_drone_df <- do.call(rbind, diff_mean_lidar_drone)

rownames(diff_mean_lidar_aircraft_df) <- c("Neukirchen8_1", "Neukirchen8_2",
                                           "Neukirchen9_1", "Neukirchen9_2",
                                           "Reinhardshagen_1", "Reinhardshagen_2")

rownames(diff_mean_lidar_drone_df) <- c("Neukirchen8_1", "Neukirchen8_2",
                                        "Neukirchen9_1", "Neukirchen9_2",
                                        "Reinhardshagen_1", "Reinhardshagen_2")


diff_means_df <- data.frame(cbind(diff_mean_lidar_drone_df, diff_mean_lidar_aircraft_df))

colnames(diff_means_df) <- c("Differenz LiDAR Drohne", "Differenz LiDAR Flugzeug")

# Plots of differences
par_org <- par()
par(mfrow = c(2,3))
terra::plot(chm_lidar_files_list[[5]],
            main = "CHM LiDAR")
terra::plot(ndsms_drone_resampled_files_list[[5]],
            main = "nDSM Drohne")
terra::plot(ndsms_aircraft_files_list[[5]],
            main = "nDSM Flugzeug")
terra::plot(diff_lidar_drone[[5]],
            col = grDevices::hcl.colors(50, palette = "blue-red 3"),
            main = "Differenz LiDAR - Drohne")
terra::plot(diff_lidar_aircraft[[5]],
            col = grDevices::hcl.colors(50, palette = "blue-red 3"),
            main = "Differenz LiDAR - Flugzeug")
par(par_org)


# Calculate RMSE between LiDAR derived CHMS
# and aircraft and drone derived nDSMs

rmse_lidar_aircraft <- mapply(function(x, y) sqrt(mean((x - y)^2)),
                              chm_lidar_files_list, ndsms_aircraft_files_list,
                              SIMPLIFY = FALSE)

rmse_mean_lidar_aircraft <- lapply(rmse_lidar_aircraft,
                                   function(x) mean(values(x, na.rm = TRUE)))


rmse_lidar_drone <- mapply(function(x, y) sqrt(mean((x - y)^2)),
                           chm_lidar_files_list, ndsms_drone_resampled_files_list,
                           SIMPLIFY = FALSE)

rmse_mean_lidar_drone <- lapply(rmse_lidar_drone,
                                function(x) mean(values(x, na.rm = TRUE)))

# Create data frame with RMSE results
rmse_mean_lidar_aircraft_df <- do.call(rbind, rmse_mean_lidar_aircraft)

rmse_mean_lidar_drone_df <- do.call(rbind, rmse_mean_lidar_drone)

rownames(rmse_mean_lidar_aircraft_df) <- c("Neukirchen8_1", "Neukirchen8_2",
                                           "Neukirchen9_1", "Neukirchen9_2",
                                           "Reinhardshagen_1", "Reinhardshagen_2")

rownames(rmse_mean_lidar_drone_df) <- c("Neukirchen8_1", "Neukirchen8_2",
                                        "Neukirchen9_1", "Neukirchen9_2",
                                        "Reinhardshagen_1", "Reinhardshagen_2")


rmse_means_df <- data.frame(cbind(rmse_mean_lidar_drone_df, rmse_mean_lidar_aircraft_df))

colnames(rmse_means_df) <- c("RMSE LiDAR Drohne", "RMSE LiDAR Flugzeug")

# Plots of RMSE
par_org <- par()
par(mfrow = c(2,3))
terra::plot(chm_lidar_files_list[[5]],
            main = "CHM LiDAR")
terra::plot(ndsms_drone_resampled_files_list[[5]],
            main = "nDSM Drohne")
terra::plot(ndsms_aircraft_files_list[[5]],
            main = "nDSM Flugzeug")
terra::plot(rmse_lidar_drone[[5]],
            col = grDevices::hcl.colors(50, palette = "RdYlBu", rev = TRUE),
            main = "RMSE LiDAR - Drohne")
terra::plot(rmse_lidar_aircraft[[5]],
            col = grDevices::hcl.colors(50, palette = "RdYlBu", rev = TRUE),
            main = "RMSE LiDAR - Flugzeug")
par(par_org)

