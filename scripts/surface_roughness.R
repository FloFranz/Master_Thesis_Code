#-------------------------------------------------------------------------------
# Name:         surface_roughness.R
# Author:       Florian Franz
# Description:  script calculates surface roughness for three different locations:
#               Reinhardshagen, Neukirchen_8 and Neukirchen_9,
#               two metrics are calculated: standard deviation and
#               percentile difference between the 5 % and 95 % percentile
#               workflow is similar to the method used in the 'F3 project'
#               (by FVA and NW-FVA) to derive forest structural parameters
#               (https://www.waldwissen.net/de/technik-und-planung/waldinventur/ableitung-von-kronendachrauigkeit)
# Data          DSM tif files from platforms 
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
DSM = TRUE
nDSM = FALSE

# Source script for DSM import and preprocessing
source("D:/scripts/src/process_DSMs.R", local = TRUE)



# Quick looks
#-------------
# Parallel plots of aircraft and drone DSMs (raw and resampled)
par_org <- par()
par(mfrow = c(1,3))
terra::plot(dsm_drone_reinhardshagen_1)
terra::plot(dsm_drone_reinhardshagen_1_resampled)
terra::plot(dsm_aircraft_reinhardshagen_1)
par(par_org)



# Write to disk
#---------------
# Define output path
out_path <- "D:/output/dsm_heights/"

# Rename DSMs
names(dsm_aircraft_reinhardshagen_1)  <- "aircraft_reinhardshagen_1"
names(dsm_aircraft_reinhardshagen_2)  <- "aircraft_reinhardshagen_2"
names(dsm_aircraft_neukirchen8_1)     <- "aircraft_neukirchen8_1"
names(dsm_aircraft_neukirchen8_2)     <- "aircraft_neukirchen8_2"
names(dsm_aircraft_neukirchen9_1)     <- "aircraft_neukirchen9_1"
names(dsm_aircraft_neukirchen9_2)     <- "aircraft_neukirchen9_2"

names(dsm_drone_reinhardshagen_1)  <- "drone_reinhardshagen_1"
names(dsm_drone_reinhardshagen_2)  <- "drone_reinhardshagen_2"
names(dsm_drone_neukirchen8_1)     <- "drone_neukirchen8_1"
names(dsm_drone_neukirchen8_2)     <- "drone_neukirchen8_2"
names(dsm_drone_neukirchen9_1)     <- "drone_neukirchen9_1"
names(dsm_drone_neukirchen9_2)     <- "drone_neukirchen9_2"

dsms_list <- list(dsm_aircraft_reinhardshagen_1, dsm_aircraft_reinhardshagen_2,
                  dsm_aircraft_neukirchen8_1, dsm_aircraft_neukirchen8_2,
                  dsm_aircraft_neukirchen9_1, dsm_aircraft_neukirchen9_2,
                  dsm_drone_reinhardshagen_1, dsm_drone_reinhardshagen_2,
                  dsm_drone_neukirchen8_1, dsm_drone_neukirchen8_2,
                  dsm_drone_neukirchen9_1, dsm_drone_neukirchen9_2)

# Change CRS of drone DSMs to ETRS89 / UTM zone 32N (EPSG:25832)
for (drone_dsm in 7:length(dsms_list)){
  
  terra::crs(dsms_list[[drone_dsm]]) <- "epsg:25832"
  
}

for (file in seq(dsms_list)){
  
  terra::writeRaster(dsms_list[[file]],
                     filename = paste0(out_path, 
                                       substr(names(dsms_list[[file]]), 1,
                                              nchar(names(dsms_list[[file]]))), ".tif"),
                     overwrite = TRUE)
  
}



# Calculate standard deviation of the height
#--------------------------------------------
# Focal moving window with size 41 x 41 pixels
# --> 40 surrounding neighbors, 0.5 m pixel size
# --> = 20.5 m size of the moving window
# Calculation considers only resampled drone DSMs
# as the moving window size varies with pixel size

dsms_list <- list(dsm_drone_reinhardshagen_1_resampled, dsm_aircraft_reinhardshagen_1,
                  dsm_drone_reinhardshagen_2_resampled, dsm_aircraft_reinhardshagen_2,
                  dsm_drone_neukirchen8_1_resampled, dsm_aircraft_neukirchen8_1,
                  dsm_drone_neukirchen8_2_resampled, dsm_aircraft_neukirchen8_2,
                  dsm_drone_neukirchen9_1_resampled, dsm_aircraft_neukirchen9_1,
                  dsm_drone_neukirchen9_2_resampled, dsm_aircraft_neukirchen9_2)

focal_sd_list <- lapply(dsms_list, 
                        FUN = function(x) terra::focal(x, w = 41, fun = "sd",
                                                       na.policy = "omit"))

# Write to disk
# Define output path
out_path <- "J:/output/surface_roughness_sd/"

names(focal_sd_list[[1]]) <- "drone_resampled_reinhardshagen_1.tif"
names(focal_sd_list[[2]]) <- "aircraft_reinhardshagen_1.tif"
names(focal_sd_list[[3]]) <- "drone_resampled_reinhardshagen_2.tif"
names(focal_sd_list[[4]]) <- "aircraft_reinhardshagen_2.tif"
names(focal_sd_list[[5]]) <- "drone_resampled_neukirchen8_1.tif"
names(focal_sd_list[[6]]) <- "aircraft_neukirchen8_1.tif"
names(focal_sd_list[[7]]) <- "drone_resampled_neukirchen8_2.tif"
names(focal_sd_list[[8]]) <- "aircraft_neukirchen8_2.tif"
names(focal_sd_list[[9]]) <- "drone_resampled_neukirchen9_1.tif"
names(focal_sd_list[[10]]) <- "aircraft_neukirchen9_1.tif"
names(focal_sd_list[[11]]) <- "drone_resampled_neukirchen9_2.tif"
names(focal_sd_list[[12]]) <- "aircraft_neukirchen9_2.tif"

for (i in seq(focal_sd_list)){
  
  terra::writeRaster(focal_sd_list[[i]],
                     filename = paste0(out_path, substr(focal_sd_list[[i]]@ptr[["names"]], 1,
                                                        nchar(focal_sd_list[[i]]@ptr[["names"]]))))
  
}

# Calculate mean
focal_sd_aircraft_list <- focal_sd_list[c(2,4,6,8,10,12)]

focal_sd_drone_list <- focal_sd_list[c(1,3,5,7,9,11)]

sd_means_aircraft <- lapply(focal_sd_aircraft_list,
                            function(x) round(mean(values(x, na.rm = TRUE)), 2))

sd_means_drone <- lapply(focal_sd_drone_list,
                         function(x) round(mean(values(x, na.rm = TRUE)), 2))

sd_means_aircraft_df <- do.call(rbind, sd_means_aircraft)

sd_means_drone_df <- do.call(rbind, sd_means_drone)

rownames(sd_means_aircraft_df) <- c("Reinhardshagen_1", "Reinhardshagen_2",
                                    "Neukirchen8_1", "Neukirchen8_2",
                                    "Neukirchen9_1", "Neukirchen9_2")

rownames(sd_means_drone_df) <- c("Reinhardshagen_1", "Reinhardshagen_2",
                                 "Neukirchen8_1", "Neukirchen8_2",
                                 "Neukirchen9_1", "Neukirchen9_2")

sd_means_df <- data.frame(cbind(sd_means_drone_df, sd_means_aircraft_df))

colnames(sd_means_df) <- c("SD Drohne", "SD Flugzeug")

# Test plots
par_org <- par()
par(mfrow = c(1,2))
terra::plot(focal_sd_list[[1]], col = viridis::viridis(50))
terra::plot(focal_sd_list[[2]], col = viridis::viridis(50))
par(par_org)



### Remains for single testing
test <- terra::focal(dsm_aircraft_neukirchen8_2, w = 41,
                     fun = "sd", na.policy = "omit")

test1 <- terra::focal(dsm_drone_neukirchen8_2_resampled, w = 41,
                     fun = "sd", na.policy = "omit")

test2 <- terra::focal(dsm_drone_neukirchen8_2, w = 41,
                      fun = "sd", na.policy = "omit")


mov_wind <- terra::focalMat(dsm_aircraft_reinhardshagen_1, d = c(20,20), type = "rect")

test <- terra::focal(dsm_aircraft_reinhardshagen_1, w = mov_wind,
                     fun = "sd", na.policy = "omit")


par(mfrow = c(1,3))
terra::plot(test, col = viridis::viridis(50))
terra::plot(test1, col = viridis::viridis(50))
terra::plot(test2, col = viridis::viridis(50))
################################################

# Calculate difference of percentile of height
#---------------------------------------------
# First, the 5 % and 95 % percentile are calculated
# Second, the interval between these is calculated
# Focal moving window with size 41 x 41 pixels

focal_p5_list <- pbapply::pblapply(dsms_list, 
                                   FUN = function(x)
                                     terra::focal(x, w = 41, fun = function(i) 
                                       quantile(i, 0.05, na.rm = TRUE),
                                       na.policy = "omit"))

focal_p59_list <- pbapply::pblapply(dsms_list, 
                                    FUN = function(x)
                                      terra::focal(x, w = 41, fun = function(i)
                                        quantile(i, 0.95, na.rm = TRUE),
                                        na.policy = "omit"))

p_diff_list <- mapply('-', focal_p59_list, focal_p5_list, SIMPLIFY = FALSE)

# Write to disk
# Define output path
out_path <- "J:/output/surface_roughness_percentile_diff/"

names(p_diff_list[[1]]) <- "drone_resampled_reinhardshagen_1.tif"
names(p_diff_list[[2]]) <- "aircraft_reinhardshagen_1.tif"
names(p_diff_list[[3]]) <- "drone_resampled_reinhardshagen_2.tif"
names(p_diff_list[[4]]) <- "aircraft_reinhardshagen_2.tif"
names(p_diff_list[[5]]) <- "drone_resampled_neukirchen8_1.tif"
names(p_diff_list[[6]]) <- "aircraft_neukirchen8_1.tif"
names(p_diff_list[[7]]) <- "drone_resampled_neukirchen8_2.tif"
names(p_diff_list[[8]]) <- "aircraft_neukirchen8_2.tif"
names(p_diff_list[[9]]) <- "drone_resampled_neukirchen9_1.tif"
names(p_diff_list[[10]]) <- "aircraft_neukirchen9_1.tif"
names(p_diff_list[[11]]) <- "drone_resampled_neukirchen9_2.tif"
names(p_diff_list[[12]]) <- "aircraft_neukirchen9_2.tif"

for (i in seq(p_diff_list)){
  
  terra::writeRaster(p_diff_list[[i]],
                     filename = paste0(out_path, substr(p_diff_list[[i]]@ptr[["names"]], 1,
                                                        nchar(p_diff_list[[i]]@ptr[["names"]]))))
  
}

# Load data
file_path_p_diff <- "J:/output/surface_roughness_percentile_diff/"

p_diff_files <- list.files(file_path_p_diff,
                           pattern = glob2rx("*.tif"),
                           full.names = TRUE)

p_diff_files_list <- lapply(p_diff_files,
                            FUN = function(x, i) terra::rast(x[i]))

p_diff_aircraft <- p_diff_files_list[c(1:6)]

p_diff_drone <- p_diff_files_list[c(7:12)]

# Calculate mean
p_diff_means_aircraft <- lapply(p_diff_aircraft,
                                function(x) mean(values(x, na.rm = TRUE)))

p_diff_means_drone <- lapply(p_diff_drone,
                             function(x) mean(values(x, na.rm = TRUE)))

p_diff_means_aircraft_df <- do.call(rbind, p_diff_means_aircraft)

p_diff_means_drone_df <- do.call(rbind, p_diff_means_drone)

rownames(p_diff_means_aircraft_df) <- c("Neukirchen8_1", "Neukirchen8_2",
                                        "Neukirchen9_1", "Neukirchen9_2",
                                        "Reinhardshagen_1", "Reinhardshagen_2")

rownames(p_diff_means_drone_df) <- c("Neukirchen8_1", "Neukirchen8_2",
                                     "Neukirchen9_1", "Neukirchen9_2",
                                     "Reinhardshagen_1", "Reinhardshagen_2")

p_diff_means_df <- data.frame(cbind(p_diff_means_drone_df, 
                                    p_diff_means_aircraft_df))

colnames(p_diff_means_df) <- c("p_diff Drohne", "p_diff Flugzeug")

# Test plots
par_org <- par()
par(mfrow = c(1,2))
terra::plot(p_diff_drone[[2]], col = grDevices::hcl.colors(50, "Zissou 1"))
terra::plot(p_diff_aircraft[[2]], col = grDevices::hcl.colors(50, "Zissou 1"))
par(par_org)



# Plots of drone and aircraft roughness metrics
# with respective DSM
par_org <- par()
par(mfrow = c(3,2))
terra::plot(dsms_list[[1]], main = "DSM Drohne")
terra::plot(dsms_list[[2]], main = "DSM Flugzeug")
terra::plot(focal_sd_list[[1]], col = viridis::viridis(50),
            main = "Standardabweichung der H?he - Drohne")
terra::plot(focal_sd_list[[2]], col = viridis::viridis(50),
            main = "Standardabweichung der H?he - Flugzeug")
terra::plot(p_diff_drone[[5]], col = grDevices::hcl.colors(50, "Zissou 1"),
            main = "Perzentilabstand der H?he - Drohne")
terra::plot(p_diff_aircraft[[5]], col = grDevices::hcl.colors(50, "Zissou 1"),
            main = "Perzentilabstand der H?he - Flugzeug")
par(par_org)



### Remains for single testing
quant_5 <- terra::focal(dsm_aircraft_reinhardshagen_1, w = 41,
                        fun = function(x) quantile(x, 0.05, na.rm = TRUE),
                        na.policy = "omit")


quant_95 <- terra::focal(dsm_aircraft_reinhardshagen_1, w = 41,
                         fun = function(x) quantile(x, c(0.95), na.rm = TRUE),
                         na.policy = "omit")

p_diff <- quant_95 - quant_5
##################################

