#-------------------------------------------------------------------------------
# Name:         set_filepaths.R
# Author:       Florian Franz
# Description:  script defines file paths for data import,
#               it sets paths for DSM and nDSM depending on what is selected in
#               the main script generated from platforms aircraft and drone,
#               three different locations: Reinhardshagen, Neukirchen_8 and Neukirchen_9
# Data          DSM and nDSM tif files
#-------------------------------------------------------------------------------

# Define root directory
root_dir <- "J:/data/"

# Define platform
aircraft <- "aircraft/"
drone <- "drone/"

# Define location
reinhardshagen <- "Reinhardshagen_e_str1_bis_4/"
neukirchen8 <- "neukirchen_e_str8/"
neukirchen9 <- "neukirchen_e_str9/"

if (DSM == TRUE & nDSM == FALSE) {
  
  # Define path for DSM
  dsm <- "DSM/"
  
  # Aircraft
  #-------------
  # Location Reinhardshagen
  filepath_reinhardshagen_dsm_aircraft <- paste0(root_dir, aircraft, reinhardshagen, dsm)
  
  # Location Neukirchen_8
  filepath_neukirchen8_dsm_aircraft <- paste0(root_dir, aircraft, neukirchen8, dsm)
  
  # Location Neukirchen_9
  filepath_neukirchen9_dsm_aircraft <- paste0(root_dir, aircraft, neukirchen9, dsm)
  
  # Drone
  #--------
  # Location Reinhardshagen
  filepath_reinhardshagen_dsm_drone <- paste0(root_dir, drone, reinhardshagen, dsm)
  
  # Location Neukirchen_8
  filepath_neukirchen8_dsm_drone <- paste0(root_dir, drone, neukirchen8, dsm)
  
  # Location Neukirchen_9
  filepath_neukirchen9_dsm_drone <- paste0(root_dir, drone, neukirchen9, dsm)
  
} else if (DSM == FALSE & nDSM == TRUE) {
  
  # Define path for nDSM
  ndsm <- "nDSM/"
  
  # Aircraft
  #-------------
  # Location Reinhardshagen
  filepath_reinhardshagen_ndsm_aircraft <- paste0(root_dir, aircraft, reinhardshagen, ndsm)
  
  # Location Neukirchen_8
  filepath_neukirchen8_ndsm_aircraft <- paste0(root_dir, aircraft, neukirchen8, ndsm)
  
  # Location Neukirchen_9
  filepath_neukirchen9_ndsm_aircraft <- paste0(root_dir, aircraft, neukirchen9, ndsm)
  
  # Drone
  #--------
  # Location Reinhardshagen
  filepath_reinhardshagen_ndsm_drone <- paste0(root_dir, drone, reinhardshagen, ndsm)
  
  # Location Neukirchen_8
  filepath_neukirchen8_ndsm_drone <- paste0(root_dir, drone, neukirchen8, ndsm)
  
  # Location Neukirchen_9
  filepath_neukirchen9_ndsm_drone <- paste0(root_dir, drone, neukirchen9, ndsm)
  
} else if (DSM == TRUE & nDSM == TRUE) {
  
  # Define path for DSM and nDSM
  dsm <- "DSM/"
  ndsm <- "nDSM/"
  
  # Aircraft
  #-------------
  # Location Reinhardshagen
  filepath_reinhardshagen_dsm_aircraft <- paste0(root_dir, aircraft, reinhardshagen, dsm)
  filepath_reinhardshagen_ndsm_aircraft <- paste0(root_dir, aircraft, reinhardshagen, ndsm)
  
  # Location Neukirchen_8
  filepath_neukirchen8_dsm_aircraft <- paste0(root_dir, aircraft, neukirchen8, dsm)
  filepath_neukirchen8_ndsm_aircraft <- paste0(root_dir, aircraft, neukirchen8, ndsm)
  
  # Location Neukirchen_9
  filepath_neukirchen9_dsm_aircraft <- paste0(root_dir, aircraft, neukirchen9, dsm)
  filepath_neukirchen9_ndsm_aircraft <- paste0(root_dir, aircraft, neukirchen9, ndsm)
  
  # Drone
  #--------
  # Location Reinhardshagen
  filepath_reinhardshagen_dsm_drone <- paste0(root_dir, drone, reinhardshagen, dsm)
  
  # Location Neukirchen_8
  filepath_neukirchen8_dsm_drone <- paste0(root_dir, drone, neukirchen8, dsm)
  
  # Location Neukirchen_9
  filepath_neukirchen9_dsm_drone <- paste0(root_dir, drone, neukirchen9, dsm)
  
} else {
  
  print("Nothing selected. Select TRUE for DSM, nDSM, or both.")
  
}