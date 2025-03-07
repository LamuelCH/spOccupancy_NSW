library(raster)
library(magrittr)
library(terra)
library(slga)
library(ggplot2)
library(sf)
library(geodata)
library(tidyr)

###########################################################################
# Create species range mask ---------------------------------------------------
library(terra)
library(dplyr)

options(mc.core = parallel::detectCores())

# Define study area
aoi_3577 = ext(684147.3, 2225320.5, -4597511.6, -2556437.4)
aoi_32755 <- ext(-100000, 1225969, 5500000, 7400000)
# aoi_polygon <- st_as_sfc(st_bbox(c(xmin = -100000, xmax = 1225969, ymin = 5500000, ymax = 7400000), crs = 32755))

# Step 3: Transform the polygon to EPSG:3577
# aoi_polygon_3577 <- st_transform(aoi_polygon, crs = 3577)

# Step 4: Extract the transformed extent
# aoi_3577 <- st_bbox(aoi_polygon_3577)
# Print the transformed extent
# print(aoi_3577)


# Read the species IUCN range shapefile
# range = vect("spatial_data/Range data/Dasyurus Maculatus/data_0.shp")
# range = project(range, "EPSG:32755")   # reproject to Easting and Northing
# range = crop(range,aoi)

# Create a raster mask with 1 km resolution
# Generate an empty raster with the defined extent and resolution
raster_3577 = rast(extent = aoi_3577, resolution = 1000, crs = "EPSG:3577")
raster_32755 = rast(extent = aoi_32755, resolution = 1000, crs = "EPSG:32755")

# mask = rasterize(range, raster, field = 1)  # Field = 1 assigns a value of 1 to the rasterized area

# Save the raster mask (optional)
# writeRaster(mask, "RStudio/spOccupancy_NSW_20241219/input/mask_1km.tif", overwrite = TRUE)

# Plot the results (optional)
# plot(mask, main = "1 km Resolution Raster Mask")


###########################################################################
# Raster processing -------------------------------------------------------
# CHELSA bioclimatic data -------------------------------------------------
## Define the directory containing the CHELSA rasters
input_dir <- "data/beta/CHELSA/climatologies/2011-2040/GFDL-ESM4/ssp126/bio/"
output_dir <- "input/"

# List all raster files in the directory
bio <- list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)

# Extract layer names from file names
layer_names <- sub(".*CHELSA_([^_]+)_.*", "\\1", basename(bio))

bio = rast(bio)

# Assign names to the layers in the stack
names(bio) <- layer_names

# Reproject the raster using bilinear method
bio <- terra::project(bio, raster_3577, method = "bilinear")

# Save the stacked raster to a file
output_file = file.path(output_dir, "CHELSA_bio_2011-2040_gfdl-esm4_ssp126_V.2.1_EPSG3577.tif")
writeRaster(bio,output_file, overwrite = TRUE)


# pr ----------------------------------------------------------------------
## Define the directory containing the CHELSA rasters
input_dir <- "spatial_data/CHELSA/climatologies/2011-2040/GFDL-ESM4/ssp126/pr/"
output_dir <- "RStudio/spOccupancy_NSW/data/beta/"

# List all raster files in the directory
pr <- list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)

# Extract layer names using the updated regex
layer_names <- sub(".*CHELSA_.*_([a-z]+_[0-9]{2}).*", "\\1", basename(pr))

pr = rast(pr)
pr = mean(pr)
# Assign names to the layers in the stack
names(pr) <- "pr_mean"

# Reproject the raster using bilinear method
pr <- terra::project(pr, raster, method = "bilinear")

# Save the stacked raster to a file
output_file = file.path(output_dir, "CHELSA_gfdl-esm4_r1i1p1f1_w5e5_ssp126_pr_mean_2011_2040_norm_EPSG3577.tif")
writeRaster(pr,output_file, overwrite = TRUE)


# tas ---------------------------------------------------------------------
## Define the directory containing the CHELSA rasters
input_dir <- "spatial_data/CHELSA/climatologies/2011-2040/GFDL-ESM4/ssp126/tas/"
output_dir <- "RStudio/spOccupancy_NSW_20241219/input/"

# List all raster files in the directory
tas <- list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)

# Extract layer names using the updated regex
layer_names <- sub(".*CHELSA_.*_([a-z]+_[0-9]{2}).*", "\\1", basename(tas))

tas = rast(tas)
tas = mean(tas)
# Assign names to the layers in the stack
names(tas) <- "tas_mean"

# Reproject the raster using bilinear method
tas <- terra::project(tas, raster, method = "bilinear")

# Save the stacked raster to a file
output_file = file.path(output_dir, "CHELSA_gfdl-esm4_r1i1p1f1_w5e5_ssp126_tas_mean_2011_2040_EPSG32755.tif")
writeRaster(tas,output_file, overwrite = TRUE)


# tasmax ------------------------------------------------------------------
## Define the directory containing the CHELSA rasters
input_dir <- "spatial_data/CHELSA/climatologies/2011-2040/GFDL-ESM4/ssp126/tasmax/"
output_dir <- "RStudio/spOccupancy_NSW_20241219/input/"

# List all raster files in the directory
tasmax <- list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)

# Extract layer names using the updated regex
layer_names <- sub(".*CHELSA_.*_([a-z]+_[0-9]{2}).*", "\\1", basename(tasmax))

tasmax = rast(tasmax)
tasmax = mean(tasmax)

# Assign names to the layers in the stack
names(tasmax) <- "tasmax_mean"

# Reproject the raster using bilinear method
tasmax <- terra::project(tasmax, raster, method = "bilinear")

# Save the stacked raster to a file
output_file = file.path(output_dir, "CHELSA_gfdl-esm4_r1i1p1f1_w5e5_ssp126_tasmax_mean_2011_2040_EPSG32755.tif")
writeRaster(tasmax,output_file, overwrite = TRUE)

# tasmin ------------------------------------------------------------------
## Define the directory containing the CHELSA rasters
input_dir <- "spatial_data/CHELSA/climatologies/2011-2040/GFDL-ESM4/ssp126/tasmin/"
output_dir <- "RStudio/spOccupancy_NSW_20241219/input/"

# List all raster files in the directory
tasmin <- list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)

tasmin = rast(tasmin)
tasmin = mean(tasmin)

# Assign names to the layers in the stack
names(tasmin) <- "tasmin_mean"

# Reproject the raster using bilinear method
tasmin <- terra::project(tasmin, raster, method = "bilinear")

# Save the stacked raster to a file
output_file = file.path(output_dir, "CHELSA_gfdl-esm4_r1i1p1f1_w5e5_ssp126_tasmin_mean_2011_2040_EPSG32755.tif")
writeRaster(tasmin,output_file, overwrite = TRUE)



# Ensure the output directory exists
# if (!dir.exists(output_dir)) {
# dir.create(output_dir, recursive = TRUE)
# }

# List all CHELSA raster files with the naming pattern
# raster_files <- list.files(input_dir, pattern = "CHELSA_swb_\\d{4}_V\\.2\\.1\\.tif$", full.names = TRUE)

# Load the mask (assumed already created)
# Loop through each raster file
# for (file in raster_files) {
# Load the raster
# bioclim <- rast(file)

# Reproject and apply the mask
# bioclim_processed <- project(bioclim, mask, method = "bilinear") * mask

# Create the output file path
# output_file <- file.path(output_dir, basename(gsub(".tif", "_masked.tif", file)))

# Save the processed raster
# writeRaster(bioclim_processed, output_file, overwrite = TRUE)

# Print progress message (optional)
# message("Processed and saved: ", output_file)
# }
# Optional: Combine into a single multi-layer raster stack
# raster_stack <- rast(processed_rasters)

# plot(bioclim_processed)
# 
###########################################################################
# FOREST COVER##########################################################################
#In this section we create the forest cover layer. Because there is obvious pixel bandings and quality issue for 2005, 2010 and 2015 layer,
#we will only use 2000 forest cover 
# WE COMBINE ALL FOREST COVER FROM 2000-2015 AND TAKE THE MEAN VALUE TO EACH PIXEL
# Get a list of all subdirectories ending with "TC_2010"
subdirs <- list.dirs("data/RawData.RAW/Forest Cover/", 
                     recursive = TRUE, 
                     full.names = TRUE)
subdirs <- subdirs[grep("TC_2000$", subdirs)]

# Initialize an empty list to store raster files
GFCC_list <- list()

# Loop through each subdirectory
for (subdir in subdirs) {
  # Get a list of raster files ending with "TC_2010.tif" in this subdirectory
  files <- list.files(path = subdir, pattern = "\\.tif$", full.names = TRUE)
  files <- files[!grepl("_err\\.tif$|_idx\\.tif$|_idx\\.txt$", files)]
  
  # Read the rasters into a list of SpatRaster objects
  GFCC_list <- c(GFCC_list, lapply(files, terra::rast))
}

# REPROJET TO SAME CRS AND RESAMPLE TO 0.01 DEGREE
reprojected_GFCC <- lapply(GFCC_list, function(r) project(r, 
                                                          raster_1km, 
                                                          method = "bilinear")) %>% 
  sprc() %>% 
  mosaic() #%>% 

mask(., env_auMask_1km)

# Filter out FOREST COVER PERCENTAGE > 100%
forestCover_1km <- ifel(reprojected_GFCC > 100,
                        NA, 
                        reprojected_GFCC)

names(forestCover_1km) <- "forestCover"

plot(forestCover_1km)

writeRaster(x = forestCover_1km, filename = "data/env_forestCover_1km.tif", overwrite = TRUE)

###########################################################################
# FOLIAGE PROJECTIVE COVER ####
foliage <- rast("data/raw/raster/Foliage Projective Cover/Woody vegetation cover - Landsat, JRSRP, Australian coverage, 2000-2010/lztmre_aus_y20002011_dm7a2_d20050630_r500m.tif")

#resample to 100m resolution
foliage_1km <- terra::project(foliage, 
                              env_tasMask_1km, 
                              method = "bilinear") * env_tasMask_1km
foliage_1km = clamp((foliage_1km - 100) / (200 - 100) * 100, 0, 100)

names(foliage_1km) = "foliageCover"

plot(foliage_1km)

writeRaster(x = foliage_1km, filename = "data/env_foliage_1km.tif", overwrite = TRUE)
###########################################################################
# DEPTH OF REGOLITH ####
der = rast("data/beta/Soil and Landscape Grid National Soil Attribute Maps/data/DER_000_999_EV_N_P_AU_NAT_C_20150601.tif")
der <- terra::project(der, raster_3577, method = "bilinear")

names(der) = "der"

writeRaster(der, filename = "input/env_der_EPSG3577.tif", overwrite = TRUE)
###########################################################################
# ROAD DENSITY####
#Road Density
roads <- vect("data/beta/road/2024_12_NationalRoads/NationalRoads2024_12.shp") 
roads = project(roads, y = "EPSG:3577")

# Convert the roads to raster with 1km resolution
roads <- terra::rasterize(roads,raster_3577, field = "Shape_Leng", fun = sum)
# roads = crop(roads, raster) 

# Calculate the total length of road within a 10 km moving window for each cell
roads_density <- terra::focal(roads, 
                              w = 9, #9 cells roughly equal to 10 km
                              fun = sum,
                              na.rm = TRUE)

# Replace NA values in roads_10km_1km with 0
roads_density[is.na(roads_density)] <- 0


# Mask roads_10km_1km with the mask raster
roads_density <- roads_density

plot(roads_density)

names(roads_density) = "roadLength"

writeRaster(x = roads_density,filename = "input/env_roadDensity_EPSG3577.tif", overwrite = TRUE)
###########################################################################
# Population Density ------------------------------------------------------
input_dir <- "data/beta/pd/2010-2020/"
output_dir <- "input/"

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define years to process
years <- 2010:2020

# Loop through each year and process the corresponding file
for (year in years) {
  # Construct file name
  input_file <- file.path(input_dir, paste0("aus_pd_", year, "_1km_UNadj.tif"))
  
  # Check if the file exists
  if (file.exists(input_file)) {
    # Load the raster
    pd <- rast(input_file)
    
    # Reproject the raster using bilinear method
    pd_reprojected <- terra::project(pd, raster_3577, method = "bilinear")
    
    # Save the reprojected raster to the output directory
    output_file <- file.path(output_dir, paste0("aus_pd_", year, "_1km_UNadj_EPSG3577.tif"))
    writeRaster(pd_reprojected, output_file, overwrite = TRUE)
    
    # Log the progress
    cat("Processed year:", year, "-> Saved to:", output_file, "\n")
  } else {
    # Log if file is missing
    cat("File missing for year:", year, "\n")
  }
}

cat("Processing complete!")

# Create disturbance index

# DATA READY FOR ANALYSIS


# Human disturbance index -------------------------------------------------
## Define the directory containing the CHELSA rasters
## Define the directory containing the CHELSA rasters
input_dir <- "data/beta/pd/2010-2020/"
output_dir <- "input/"

# List all raster files in the directory
pd <- list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)

pd = rast(pd)
pd = mean(pd)


# Reproject the raster using bilinear method
pd <- terra::project(pd, raster_3577, method = "bilinear")

# Calculate the total density of population within a 10 km moving window for each cell
pd_smooth <- terra::focal(pd, w = 9, fun = sum, na.rm = TRUE)


# Assign names to the layers in the stack
names(pd_smooth) <- "pd_mean"


# Save the stacked raster to a file
output_file = file.path(output_dir, "aus_pd_2010-2020_1km_UNadj_EPSG3577.tif")
writeRaster(pd_smooth,output_file, overwrite = TRUE)



# disturbance = roads_density*pd_smooth

# Pearson Correlation Analysis --------------------------------------------


###########################################################################
bio = rast("RStudio/spOccupancy_NSW_20241219/input/CHELSA_bio_2011-2040_gfdl-esm4_ssp126_V.2.1_EPSG32755.tif")
bio = bio[[1:19]]
pr_mean = rast("RStudio/spOccupancy_NSW_20241219/input/CHELSA_gfdl-esm4_r1i1p1f1_w5e5_ssp126_pr_mean_2011_2040_norm_EPSG32755.tif")
tas_mean = rast("RStudio/spOccupancy_NSW_20241219/input/CHELSA_gfdl-esm4_r1i1p1f1_w5e5_ssp126_tas_mean_2011_2040_EPSG32755.tif")
tasmax_mean = rast("RStudio/spOccupancy_NSW_20241219/input/CHELSA_gfdl-esm4_r1i1p1f1_w5e5_ssp126_tasmax_mean_2011_2040_EPSG32755.tif")
tasmin_mean = rast("RStudio/spOccupancy_NSW_20241219/input/CHELSA_gfdl-esm4_r1i1p1f1_w5e5_ssp126_tasmin_mean_2011_2040_EPSG32755.tif")
roads_density = rast("RStudio/spOccupancy_NSW_20241219/input/env_roadDensity_EPSG32755.tif")
pd = rast("RStudio/spOccupancy_NSW_20241219/input/aus_pd_2010-2020_1km_UNadj_EPSG32755.tif")
der = rast("RStudio/spOccupancy_NSW_20241219/input/env_der_EPSG32755.tif")


env_stack <- c(bio, pr_mean, tas_mean, tasmax_mean, tasmin_mean, roads_density, pd, der)
correlation_matrix <- layerCor(env_stack, "pearson", na.rm = TRUE)

write.csv(correlation_matrix, file = "RStudio/spOccupancy_NSW_20241219/output/correlation_matrix.csv", row.names = TRUE)

###########################################################################
# install.packages("corrplot")
library(corrplot)

corrplot(correlation_matrix$correlation, method = "circle", type = "lower", insig = "blank", addCoef.col = "black", number.cex = 0.6)
# Assuming your correlation matrix is stored as a matrix or data.frame
# For example: correlation_matrix <- layerCor(your_raster, "pearson", na.rm = TRUE)$pearson

# Convert correlation matrix to a dissimilarity matrix (1 - |correlation|)
dissimilarity <- as.dist(1 - abs(correlation_matrix$correlation))

# Perform hierarchical clustering
hc <- hclust(dissimilarity, method = "average")  # "average" or "complete" linkage works well

# Cut the dendrogram to group variables with |r| > 0.7
groups <- cutree(hc, h = 0.3)  # h = 1 - 0.7 = 0.3 threshold

# View the groups
correlation_groups <- split(names(groups), groups)
print(correlation_groups)

corrplot(correlation_matrix$correlation, method = "square", order = "hclust", hclust.method = "average", 
         addrect = length(correlation_groups), tl.cex = 0.7)


env_stack = env_stack[[c("bio5","bio12", "roadLength", "pd_mean", "der")]]
correlation_matrix <- layerCor(env_stack, "pearson", na.rm = TRUE)
corrplot(correlation_matrix$correlation, method = "circle", type = "lower", insig = "blank", addCoef.col = "black", number.cex = 0.6)



# Add all your raster layers here

# Repeat this for all your raster layers
###########################################################################
