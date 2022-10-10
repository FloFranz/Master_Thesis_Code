#-------------------------------------------------------------------------------
# Name:         lidar.R
# Author:       Florian Franz
# Description:  script processes LiDAR data used as reference for
#               nDSMs from platforms aircraft and drone generated via image processing,
#               CHMs are calculated for the corresponding locations in 0.5 m resolution
# Data          las files,
#               nDSM tif files from platform aircraft (0.5 m resolution)
#-------------------------------------------------------------------------------

# Load packages
#---------------
library(terra)
library(lidR)



# Read las files
#----------------
# Define file path
file_path <- "N:/ALS-NW-FVA/"

las_files <- list.files(file_path,
                        pattern = glob2rx("*.las"),
                        full.names = TRUE)


# The las files of interest (of the investigated locations)
las_files_numbers <- c(10, 14, 27, 28, 33, 34, 78, 79, 86)

# Read all las files of interest in one list 
las_files_list <- c()

for (las_file in las_files_numbers){
  
  las_files_list <- append(las_files_list, lidR::readLAS(las_files[las_file]))
  
}

# Assign CRS
epsg_number <- 25832

for (i in seq(las_files_list)){
  
  crs(las_files_list[[i]]) <- epsg_number
  
}

# Check individual las files
lidR::las_check(las_files_list[[1]])



# This part is only required if location Reinhardshagen should be processed,
# which is not done here.
#
# Filter point-cloud of location Reinhardshagen,
# there are some outliers in the Z-dimension

# Visible in the plot
lidR::plot(las_files_list[[7]], size = 3)

min(las_files_list[[7]]@data[["Z"]])  # 278.1 --> this is correct
max(las_files_list[[7]]@data[["Z"]])  # 760.824

min(las_files_list[[8]]@data[["Z"]])  # -2581.703
max(las_files_list[[8]]@data[["Z"]])  # 783.501

min(las_files_list[[9]]@data[["Z"]])  # -13.979
max(las_files_list[[9]]@data[["Z"]])  # 790.651

# Filter all points where Z is smaller than 420 m,
# --> Max in aircraft nDSM is approx. 416 m,
# Minimum is selected individual for each tile to get correct value
las_files_list[[7]] <- lidR::filter_poi(las_files_list[[7]], Z >= 0, Z <= 420)

max(las_files_list[[7]]@data[["Z"]])  # 412.927

las_files_list[[8]] <- lidR::filter_poi(las_files_list[[8]], Z >= 53, Z <= 420)

min(las_files_list[[8]]@data[["Z"]])  # 350.223
max(las_files_list[[8]]@data[["Z"]])  # 419.309

las_files_list[[9]] <- lidR::filter_poi(las_files_list[[9]], Z >= 46, Z <= 420)

min(las_files_list[[9]]@data[["Z"]])  # 271.275
max(las_files_list[[9]]@data[["Z"]])  # 416.242



# Create nDSMs (CHMs)
#---------------------
# Height normalization within the point cloud
norm_las_list <- lapply(las_files_list,
                        FUN = function(x) lidR::normalize_height(x, knnidw()))

# Check if all ground points are 0
for (norm_las in norm_las_list){
  
  hist(filter_ground(norm_las)$Z, breaks = seq(-0.5, 0.5, 0.01), main = "", xlab = "Elevation")
  
}

# Calculate CHM
chm_list <- lapply(norm_las_list,
                   FUN = function(x) 
                     lidR::rasterize_canopy(x, res = 0.5,
                                            algorithm = pitfree(thresholds = c(0, 10, 20), 
                                                                max_edge = c(0, 1.5))))

# Test plots
lidR::plot(chm_list[[1]], col = lidR::height.colors(50))

# Merge tiles into one
chm_reinhardshagen_merged <- terra::merge(chm_list[[7]],
                                          chm_list[[8]],
                                          chm_list[[9]])

chm_neukirchen8_1_merged <- terra::merge(chm_list[[3]],
                                         chm_list[[4]])

chm_neukirchen8_2_merged <- terra::merge(chm_list[[5]],
                                         chm_list[[6]])

chm_neukirchen9_merged <- terra::merge(chm_list[[1]],
                                       chm_list[[2]])

# Crop out subsets of the merged LiDAR CHMs
# corresponding to the extent of the respective aircraft and drone nDSMs,
# here, aircraft nDSMs extent is used

# Read in aircraft nDSMs
file_path <- "J:/output/tree_heights"

aircraft_tif_files <- list.files(file_path,
                                 pattern = glob2rx("*.tif"),
                                 full.names = TRUE)

ndsms_aircraft_list <- lapply(aircraft_tif_files, FUN = function(x) terra::rast(x))

# Crop out subsets
chm_reinhardshagen_1_subset <- terra::crop(chm_reinhardshagen_merged,
                                           ndsms_aircraft_list[[5]])

chm_reinhardshagen_2_subset <- terra::crop(chm_reinhardshagen_merged,
                                           ndsms_aircraft_list[[6]])

chm_neukirchen8_1_subset <- terra::crop(chm_neukirchen8_1_merged,
                                        ndsms_aircraft_list[[1]])

chm_neukirchen8_2_subset <- terra::crop(chm_neukirchen8_2_merged,
                                        ndsms_aircraft_list[[2]])

chm_neukirchen9_1_subset <- terra::crop(chm_neukirchen9_merged,
                                        ndsms_aircraft_list[[3]])

chm_neukirchen9_2_subset <- terra::crop(chm_neukirchen9_merged,
                                        ndsms_aircraft_list[[4]])

# Replace negative values by zero in CHMs
chm_reinhardshagen_1_subset[chm_reinhardshagen_1_subset < 0] <- 0
chm_reinhardshagen_2_subset[chm_reinhardshagen_2_subset < 0] <- 0

chm_neukirchen8_1_subset[chm_neukirchen8_1_subset < 0] <- 0
chm_neukirchen8_2_subset[chm_neukirchen8_2_subset < 0] <- 0

chm_neukirchen9_1_subset[chm_neukirchen9_1_subset < 0] <- 0
chm_neukirchen9_2_subset[chm_neukirchen9_2_subset < 0] <- 0

# Mask out the CHMs based on NA values in the aircraft nDSMs
out_path <- "J:/output/lidar_CHMs/"

chm_reinhardshagen_1 <- terra::mask(chm_reinhardshagen_1_subset,
                                    ndsms_aircraft_list[[5]],
                                    filename = paste0(out_path, "reinhardshagen_1.tif"),
                                    overwrite = TRUE)

chm_reinhardshagen_2 <- terra::mask(chm_reinhardshagen_2_subset,
                                    ndsms_aircraft_list[[6]],
                                    filename = paste0(out_path, "reinhardshagen_2.tif"),
                                    overwrite = TRUE)

chm_neukirchen8_1 <- terra::mask(chm_neukirchen8_1_subset,
                                 ndsms_aircraft_list[[1]],
                                 filename = paste0(out_path, "neukirchen8_1.tif"),
                                 overwrite = TRUE)

chm_neukirchen8_2 <- terra::mask(chm_neukirchen8_2_subset,
                                 ndsms_aircraft_list[[2]],
                                 filename = paste0(out_path, "neukirchen8_2.tif"),
                                 overwrite = TRUE)

chm_neukirchen9_1 <- terra::mask(chm_neukirchen9_1_subset,
                                 ndsms_aircraft_list[[3]],
                                 filename = paste0(out_path, "neukirchen9_1.tif"),
                                 overwrite = TRUE)

chm_neukirchen9_2 <- terra::mask(chm_neukirchen9_2_subset,
                                 ndsms_aircraft_list[[4]],
                                 filename = paste0(out_path, "neukirchen9_2.tif"),
                                 overwrite = TRUE)

# Plots
par_org <- par()
par(mfrow = c(1,2))
terra::plot(ndsms_aircraft_list[[1]])
lidR::plot(chm_neukirchen8_1)
par(par_org)