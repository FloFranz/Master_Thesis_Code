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

las <- lidR::readLAS(las_files[78])

epsg_number <- 25832

crs(las) <- epsg_number

las



# Create nDSM (CHM)
#---------------------
# Height normalisation within the point cloud
norm_las <- lidR::normalize_height(las, lidR::knnidw())

# Check if all ground points are 0
hist(filter_ground(norm_las)$Z, breaks = seq(-0.45, 0.45, 0.01), main = "", xlab = "Elevation")

# Calculate CHM
chm <- lidR::grid_canopy(norm_las, res = 1.0, lidR::p2r())

plot(chm, col = lidR::height.colors(10))













