#-------------------------------------------------------------------------------
# Name:         lidar.R
# Author:       Florian Franz
# Description:  script processes LiDAR data used as reference for
#               nDSMs from platforms aircraft and drone generated via image processing
# Data          nDSM tif files from platforms 
#               aircraft (0.5 m resolution) and drone (0.1 m resolution)
#-------------------------------------------------------------------------------

# Load packages
#---------------
library(terra)
library(raster)
library(lidR)



# Read las files
#----------------
# Define file path
file_path <- "K:/ALS-NW-FVA/"

las_files <- list.files(file_path,
                        pattern = glob2rx("*.las"),
                        full.names = TRUE)


###
# The las files of interest (of the investigated locations)
las_files_numbers <- c(10, 14, 27, 28, 33, 34, 78, 79, 86)

# Use these numbers to exclude Reinhardshagen, which doesn't work yet
las_files_numbers <- c(10, 14, 27, 28, 33, 34)

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

# Cross section plot
p1 <- c(541000, 5698000)
p2 <- c(542000, 5699000)
las_tr <- clip_transect(las, p1, p2, width = 4, xz = TRUE)

ggplot(las_tr@data, aes(X,Z, color = Z)) + 
  geom_point(size = 0.5) + 
  coord_equal() + 
  theme_minimal() +
  scale_color_gradientn(colours = height.colors(50))




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
lidR::plot(chm_list[[4]], col = lidR::height.colors(50))

# Merge tiles into one
chm_neukirchen8_1_merged <- terra::merge(chm_list[[3]],
                                         chm_list[[4]])

chm_neukirchen8_2_merged <- terra::merge(chm_list[[5]],
                                         chm_list[[6]])

chm_neukirchen9_merged <- terra::merge(chm_list[[1]],
                                       chm_list[[2]])

# Test plots
lidR::plot(chm_neukirchen9_merged, col = lidR::height.colors(50))

# Crop out subsets of the merged LiDAR CHMs
# corresponding to the extent of the respective aircraft and drone nDSMs,
# here, aircraft nDSMs extent is used

file_path <- "J:/output/tree_heights"

aircraft_tif_files <- list.files(file_path,
                                 pattern = glob2rx("*.tif"),
                                 full.names = TRUE)

ndsm_aircraft9_1 <- terra::rast(aircraft_tif_files[3])

chm_neukirchen9_1_subset <- terra::crop(chm_neukirchen9_merged,
                                 ndsm_aircraft9_1)

# Mask out the CHM based on NA values in the aircraft nDSMs
chm_neukirchen9_1 <- terra::mask(chm_neukirchen9_1_subset,
                                 ndsm_aircraft9_1)

par(mfrow = c(1,2))
terra::plot(ndsm_aircraft9_1)
lidR::plot(chm_neukirchen9_1)

diff <- chm_neukirchen9_1 - ndsm_aircraft9_1

par(mfrow = c(1,1))
terra::plot(diff)












################################################
# Remains for single read in
las <- lidR::readLAS(las_files[86])
las_right <- lidR::readLAS(las_files[28])

epsg_number <- 25832

crs(las) <- epsg_number

lidR::las_check(las)

las

#####
# Trying to solve the problem with location Reinhardshagen
gnd <- filter_ground(las)
plot(gnd, size = 3, bg = "white", color = "Classification") 



test <- lidR::filter_duplicates(las)
plot(test, size = 3, bg = "white", color = "Classification")
######

# Height normalization within the point cloud
norm_las <- lidR::normalize_height(las, lidR::knnidw())

# Check if all ground points are 0
hist(filter_ground(norm_las)$Z, breaks = seq(-0.45, 0.45, 0.01), main = "", xlab = "Elevation")


dsm <- lidR::rasterize_canopy(las, res = 1, algorithm = p2r())

# Calculate CHM
chm <- lidR::rasterize_canopy(norm_las, res = 0.5, algorithm = pitfree(thresholds = c(0, 10, 20), max_edge = c(0, 1.5)))
chm

lidR::plot(chm, col = lidR::height.colors(50))
################################################












