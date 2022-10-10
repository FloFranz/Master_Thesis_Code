#-------------------------------------------------------------------------------
# Name:         tree_heights.R
# Author:       Florian Franz
# Description:  script derives tree heights using nDSMs from three different locations:
#               Reinhardshagen, Neukirchen_8 and Neukirchen_9,
#               basic statistics and difference between aircraft and drone nDSMs are calculated,
#               aircraft and drone nDSMs are further compared with LiDAR dervid CHMs and
#               mean error (ME) and root-mean-square error (RMSE) are calculated
# Data          nDSM tif files from platforms 
#               aircraft (0.5 m resolution) and drone (0.1 m resolution),
#               LiDAR derived CHM tif files (0.5 m resolution) 
#-------------------------------------------------------------------------------

# Load packages
#---------------
library(terra)



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
terra::plot(ndsm_drone_neukirchen9_2)
par(par_org)



# Write to disk
#---------------
# Define output path
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
terra::plot(chm_lidar_files_list[[4]])              # LiDAR CHM
terra::plot(ndsms_drone_files_list[[4]])            # Drone
terra::plot(ndsms_aircraft_files_list[[4]])         # Aircraft
par(par_org)



# Mean and maximum heights

# Aircraft
list_means_aircraft <- c()
list_maxs_aircraft <- c()

for (i in ndsms_aircraft_files_list){
  
  list_means_aircraft <- append(list_means_aircraft,
                                round(terra::global(i, fun = "mean", na.rm = TRUE), 1))
  
  list_maxs_aircraft <- append(list_maxs_aircraft,
                               round(terra::global(i, fun = "max", na.rm = TRUE), 1))
  
}

# Drone
list_means_drone <- c()
list_maxs_drone <- c()

for (i in ndsms_drone_files_list){
  
  list_means_drone <- append(list_means_drone,
                             round(terra::global(i, fun = "mean", na.rm = TRUE), 1))
  
  list_maxs_drone <- append(list_maxs_drone,
                            round(terra::global(i, fun = "max", na.rm = TRUE), 1))
  
}

# LiDAR
list_means_lidar <- c()
list_maxs_lidar <- c()

for (i in chm_lidar_files_list){
  
  list_means_lidar <- append(list_means_lidar,
                             round(terra::global(i, fun = "mean", na.rm = TRUE), 1))
  
  list_maxs_lidar <- append(list_maxs_lidar,
                            round(terra::global(i, fun = "max", na.rm = TRUE), 1))
  
}

# Create data frame of the results
df_means <- do.call(cbind, c(list_means_aircraft, list_means_drone, list_means_lidar))
df_maxs <- do.call(cbind, c(list_maxs_aircraft, list_maxs_drone, list_maxs_lidar))

locations <- c("Neukirchen8_1", "Neukirchen8_2", "Neukirchen9_1",
               "Neukirchen9_2","Reinhardshagen_1", "Reinhardshagen_2",
               "Neukirchen8_1", "Neukirchen8_2", "Neukirchen9_1",
               "Neukirchen9_2","Reinhardshagen_1", "Reinhardshagen_2",
               "Neukirchen8_1", "Neukirchen8_2", "Neukirchen9_1",
               "Neukirchen9_2","Reinhardshagen_1", "Reinhardshagen_2")

df_means_split <- split(df_means, locations)
df_maxs_split <- split(df_maxs, locations)

df_mean_heights <- do.call(rbind, df_means_split)
df_max_heights <- do.call(rbind, df_maxs_split)

colnames(df_mean_heights) <- c("Aircraft", "Drone", "LiDAR")
colnames(df_max_heights) <- c("Aircraft", "Drone", "LiDAR")



# Difference between aircraft and drone nDSMs
# For this, resampled drone nDSMs must be used 
diff_aircraft_drone <-  mapply(function(x, y) (x - y),
                               ndsms_aircraft_files_list,
                               ndsms_drone_resampled_files_list,
                               SIMPLIFY = FALSE)

# Write to disk
out_path <- "J:/output/tree_heights_diff/"

for (file in seq(diff_aircraft_drone)){
  
  terra::writeRaster(diff_aircraft_drone[[file]],
                     filename = paste0(out_path, 
                                       substr(names(diff_aircraft_drone[[file]]), 6,
                                              nchar(names(diff_aircraft_drone[[file]]))), ".tif"),
                     overwrite = TRUE)
  
}

# Plot difference maps
par(par_org)
par(mfrow = c(1,3))

terra::plot(diff_aircraft_drone[[5]], col = grDevices::hcl.colors(50, palette = "blue-red 3"))
terra::plot(ndsms_drone_resampled_files_list[[5]])
terra::plot(ndsms_aircraft_files_list[[5]])

terra::plot(diff_aircraft_drone[[6]], col = grDevices::hcl.colors(50, palette = "blue-red 3"))
terra::plot(ndsms_drone_resampled_files_list[[6]])
terra::plot(ndsms_aircraft_files_list[[6]])

terra::plot(diff_aircraft_drone[[1]], col = grDevices::hcl.colors(50, palette = "blue-red 3"))
terra::plot(ndsms_drone_resampled_files_list[[1]])
terra::plot(ndsms_aircraft_files_list[[1]])

terra::plot(diff_aircraft_drone[[2]], col = grDevices::hcl.colors(50, palette = "blue-red 3"))
terra::plot(ndsms_drone_resampled_files_list[[2]])
terra::plot(ndsms_aircraft_files_list[[2]])

terra::plot(diff_aircraft_drone[[3]], col = grDevices::hcl.colors(50, palette = "blue-red 3"))
terra::plot(ndsms_drone_resampled_files_list[[3]])
terra::plot(ndsms_aircraft_files_list[[3]])

terra::plot(diff_aircraft_drone[[4]], col = grDevices::hcl.colors(50, palette = "blue-red 3"))
terra::plot(ndsms_drone_resampled_files_list[[4]])
terra::plot(ndsms_aircraft_files_list[[4]])



# Calculate ME between LiDAR derived CHMS
# and aircraft and drone derived nDSMs
me_lidar_aircraft <-  mapply(function(x, y) mean(x - y),
                             ndsms_aircraft_files_list, chm_lidar_files_list,
                             SIMPLIFY = FALSE)

me_mean_lidar_aircraft <- lapply(me_lidar_aircraft,
                                 function(x) round(mean(values(x, na.rm = TRUE)), 2))

me_lidar_drone <-  mapply(function(x, y) mean(x - y),
                          ndsms_drone_resampled_files_list, chm_lidar_files_list,
                          SIMPLIFY = FALSE)

me_mean_lidar_drone <- lapply(me_lidar_drone,
                              function(x) round(mean(values(x, na.rm = TRUE)), 2))


me_mean_lidar_aircraft_df <- do.call(rbind, me_mean_lidar_aircraft)

me_mean_lidar_drone_df <- do.call(rbind, me_mean_lidar_drone)

rownames(me_mean_lidar_aircraft_df) <- c("Neukirchen8_1", "Neukirchen8_2",
                                         "Neukirchen9_1", "Neukirchen9_2",
                                         "Reinhardshagen_1", "Reinhardshagen_2")

rownames(me_mean_lidar_drone_df) <- c("Neukirchen8_1", "Neukirchen8_2",
                                      "Neukirchen9_1", "Neukirchen9_2",
                                      "Reinhardshagen_1", "Reinhardshagen_2")

me_means_df <- data.frame(cbind(me_mean_lidar_drone_df, me_mean_lidar_aircraft_df))

colnames(me_means_df) <- c("ME LiDAR Drohne", "ME LiDAR Flugzeug")

# Write to disk
out_path <- "J:/output/tree_heights_me/"

names(me_lidar_aircraft[[1]]) <- "aircraft_neukirchen8_1"
names(me_lidar_aircraft[[2]]) <- "aircraft_neukirchen8_2"
names(me_lidar_aircraft[[3]]) <- "aircraft_neukirchen9_1"
names(me_lidar_aircraft[[4]]) <- "aircraft_neukirchen9_2"
names(me_lidar_aircraft[[5]]) <- "aircraft_reinhardshagen_1"
names(me_lidar_aircraft[[6]]) <- "aircraft_reinhardshagen_2"

for (file in seq(me_lidar_aircraft)){
  
  terra::writeRaster(me_lidar_aircraft[[file]],
                     filename = paste0(out_path, 
                                       substr(names(me_lidar_aircraft[[file]]), 1,
                                              nchar(names(me_lidar_aircraft[[file]]))), ".tif"),
                     overwrite = TRUE)
  
}

names(me_lidar_drone[[1]]) <- "drone_neukirchen8_1"
names(me_lidar_drone[[2]]) <- "drone_neukirchen8_2"
names(me_lidar_drone[[3]]) <- "drone_neukirchen9_1"
names(me_lidar_drone[[4]]) <- "drone_neukirchen9_2"
names(me_lidar_drone[[5]]) <- "drone_reinhardshagen_1"
names(me_lidar_drone[[6]]) <- "drone_reinhardshagen_2"

for (file in seq(me_lidar_drone)){
  
  terra::writeRaster(me_lidar_drone[[file]],
                     filename = paste0(out_path, 
                                       substr(names(me_lidar_drone[[file]]), 1,
                                              nchar(names(me_lidar_drone[[file]]))), ".tif"),
                     overwrite = TRUE)
  
}

# Plots of ME
par_org <- par()
par(mfrow = c(2,3))
terra::plot(chm_lidar_files_list[[3]],
            main = "CHM LiDAR")
terra::plot(ndsms_drone_resampled_files_list[[3]],
            main = "nDSM Drohne")
terra::plot(ndsms_aircraft_files_list[[3]],
            main = "nDSM Flugzeug")
terra::plot(me_lidar_drone[[3]],
            col = grDevices::hcl.colors(50, palette = "blue-red 3"),
            main = "Mean Error LiDAR - Drohne")
terra::plot(me_lidar_aircraft[[3]],
            col = grDevices::hcl.colors(50, palette = "blue-red 3"),
            main = "Mean Error LiDAR - Flugzeug")
par(par_org)

# Calculate RMSE between LiDAR derived CHMS
# and aircraft and drone derived nDSMs
rmse_lidar_aircraft <- mapply(function(x, y) sqrt(mean((x - y)^2)),
                              ndsms_aircraft_files_list, chm_lidar_files_list,
                              SIMPLIFY = FALSE)

rmse_mean_lidar_aircraft <- lapply(rmse_lidar_aircraft,
                                   function(x) round(mean(values(x, na.rm = TRUE)), 2))

rmse_lidar_drone <- mapply(function(x, y) sqrt(mean((x - y)^2)),
                           ndsms_drone_resampled_files_list, chm_lidar_files_list,
                           SIMPLIFY = FALSE)

rmse_mean_lidar_drone <- lapply(rmse_lidar_drone,
                                function(x) round(mean(values(x, na.rm = TRUE)), 2))

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

# Write to disk
out_path <- "J:/output/tree_heights_rmse/"

names(rmse_lidar_aircraft[[1]]) <- "aircraft_neukirchen8_1"
names(rmse_lidar_aircraft[[2]]) <- "aircraft_neukirchen8_2"
names(rmse_lidar_aircraft[[3]]) <- "aircraft_neukirchen9_1"
names(rmse_lidar_aircraft[[4]]) <- "aircraft_neukirchen9_2"
names(rmse_lidar_aircraft[[5]]) <- "aircraft_reinhardshagen_1"
names(rmse_lidar_aircraft[[6]]) <- "aircraft_reinhardshagen_2"

for (file in seq(rmse_lidar_aircraft)){
  
  terra::writeRaster(rmse_lidar_aircraft[[file]],
                     filename = paste0(out_path, 
                                       substr(names(rmse_lidar_aircraft[[file]]), 1,
                                              nchar(names(rmse_lidar_aircraft[[file]]))), ".tif"),
                     overwrite = TRUE)
  
}

names(rmse_lidar_drone[[1]]) <- "drone_neukirchen8_1"
names(rmse_lidar_drone[[2]]) <- "drone_neukirchen8_2"
names(rmse_lidar_drone[[3]]) <- "drone_neukirchen9_1"
names(rmse_lidar_drone[[4]]) <- "drone_neukirchen9_2"
names(rmse_lidar_drone[[5]]) <- "drone_reinhardshagen_1"
names(rmse_lidar_drone[[6]]) <- "drone_reinhardshagen_2"

for (file in seq(rmse_lidar_drone)){
  
  terra::writeRaster(rmse_lidar_drone[[file]],
                     filename = paste0(out_path, 
                                       substr(names(rmse_lidar_drone[[file]]), 1,
                                              nchar(names(rmse_lidar_drone[[file]]))), ".tif"),
                     overwrite = TRUE)
  
}

# Plots of RMSE
par_org <- par()
par(mfrow = c(2,3))
terra::plot(chm_lidar_files_list[[3]],
            main = "CHM LiDAR")
terra::plot(ndsms_drone_resampled_files_list[[3]],
            main = "nDSM Drohne")
terra::plot(ndsms_aircraft_files_list[[3]],
            main = "nDSM Flugzeug")
terra::plot(rmse_lidar_drone[[3]],
            col = grDevices::hcl.colors(50, palette = "RdYlBu", rev = TRUE),
            main = "RMSE LiDAR - Drohne")
terra::plot(rmse_lidar_aircraft[[3]],
            col = grDevices::hcl.colors(50, palette = "RdYlBu", rev = TRUE),
            main = "RMSE LiDAR - Flugzeug")
par(par_org)