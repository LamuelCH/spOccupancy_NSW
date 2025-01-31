# Load required libraries
library(dplyr)
library(lubridate)
library(tidyr)
library(ggmap)
library(spOccupancy)
library(terra)

# Setting Project Environment ---------------------------------------------
# Set the path to the root 'data' directory
setwd("D:/")
root_path <- "RStudio/spOccupancy_NSW_20241219/data"

# Fetch species occurrence records ------------------------------------------------------------
# List all files that start with 'Atlas_records' in the 'data/BioNet*/' folder structure
file_list <- list.files(path = root_path, 
                        pattern = "^Atlas_records", 
                        recursive = TRUE, 
                        full.names = TRUE)

# Read all files into a list of DataFrames and ensure SpeciesCode is a character
data_list <- lapply(file_list, function(file) {
  df <- read.delim(file, 
                   header = TRUE, 
                   sep = "\t", 
                   stringsAsFactors = FALSE)
  df$SpeciesCode <- as.character(df$SpeciesCode) # Convert SpeciesCode to character
  return(df)
})

# Combine all the files into a single df
df = bind_rows(data_list)

# Convert df into standard format -------------------------
df = df %>%
  filter(DatasetName == "Wild Count Fauna") %>%
  select(CommonName, FamilyName, ScientificName, DateLast, LocationKey, Latitude_GDA94, Longitude_GDA94, Zone, Easting, Northing) %>% 
  mutate(
    # Parse the DateLast column into a datetime object
    Date = dmy_hms(DateLast), 
    
    # Extract components
    year = year(Date),
    month = month(Date),
    day = day(Date),
    hour = hour(Date),
    minute = minute(Date)
  ) 

df = df %>% 
  filter(year %in% 2012:2021) %>% 
  select(-DateLast, -Zone, -Date, -day, -hour, -minute) %>% # Drop unwanted columns
  rename(
    Species = CommonName,      # Rename CommonName to Species
    plot = LocationKey,
    x = Easting, # Rename Latitude_GDA94 to Latitude
    y = Northing # Rename Longitude_GDA94 to Longitude
  ) %>% 
  na.omit()

# Assign sequential plotID for unique Easting-Northing combinations
df <- df %>%
  arrange(x, y) %>%  # Arrange by Easting and Northing to ensure order
  mutate(plotID = dense_rank(paste(x, y, sep = "_")))

y.long = df %>% 
  mutate(replicate =  ((year - 2012) %/% 2) + 1) #Map 2012-2021 to 1-10


y.long <- y.long %>%
  group_by(plotID, replicate, Species) %>%
  summarize(occ = ifelse(n() > 0, 1, 0)) %>%
  ungroup() %>%
  glimpse()

# y.long = y.long %>% 
  # mutate(replicate =  ((year - 2012) %/% 2) + 1) #Map 2012-2021 to 1-10

# [optional] Visualise species occurrence data ----------------------------
# visualize the data so to have a glimpse on how species distributed
# recall online map registration
# register_stadiamaps("1794b86f-c7ad-4aab-82b4-23ae3e2f039d", write = TRUE)
# occ.plot = qmplot(Longitude_GDA94, Latitude_GDA94, data = df, zoom = 1, maptype = "stamen_terrain", darken = c(0.5, "black"), color = FamilyName)
# occ.plot


# Extract Detection-nondetection data ------------------------------------------
# Species codes, adjust accordingly
sp.codes <- sort(unique(y.long$Species))

# Plot (site) codes.
plot.codes <- sort(unique(y.long$plotID))

# Number of species
N <- length(sp.codes)

# Maximum number of replicates (no. of years) at a site
K <- length(unique(y.long$replicate))  # from 2012 to 2021

# Number of sites
J <- length(unique(y.long$plotID))

# Create array for detection-nondetection data. 
y <- array(NA, dim = c(N, J, K))

# Label the dimensions of y (not necessary, but helpful)
dimnames(y)[[1]] <- sp.codes
dimnames(y)[[2]] <- plot.codes

# Look at the structure of our array y
str(y)

# Loop through each plotID and replicate (year)
for (site in plot.codes) {
  for (rep in 1:K) {
    # Subset data for the current plotID and replicate
    sampled_data <- y.long %>%
      filter(plotID == site, replicate == rep)
    
    # Check if the plotID was sampled during this replicate
    if (nrow(sampled_data) > 0) {
      # Loop through each species
      for (species in sp.codes) {
        # Check if this species was observed
        was_observed <- any(sampled_data$Species == species)
        
        # Populate the array
        sp_index <- which(dimnames(y)[[1]] == species)
        site_index <- which(dimnames(y)[[2]] == site)
        y[sp_index, site_index, rep] <- ifelse(was_observed, 1, 0)
      }
    }
  }
}
str(y)

# Check the total number of observations for each species
apply(y, 1, sum, na.rm = TRUE)



# Format detection covariates ---------------------------------------------
month.df = df %>% 
  mutate(replicate =  ((year - 2012) %/% 2) + 1) %>%  #Map 2012-2021 to 1-10
  group_by(plotID, replicate) %>% 
  summarize(effort.months = n_distinct(paste(year, month, sep = "_"))) %>%  # Count unique (year, month) pairs
  ungroup() %>% 
  glimpse

# convert to numeric format
month.df <- month.df %>%
  mutate(across(everything(), as.numeric))

# Check the structure of the updated data frame
str(month.df)

# Detection covariates 1: survey effort
# Get unique plot IDs and years
replicate <- sort(unique(month.df$replicate))         # Years (columns)

# Initialize deployment matrix with NAs
effort.matrix <- matrix(NA, nrow = J, ncol = K, 
                            dimnames = list(plot.codes, replicate))

# Fill the matrix
for (j in 1:J) {
  for (k in 1:K) {
    # Get deployment length for site j and year k
    val <- month.df %>%
      filter(plotID == plot.codes[j], replicate == replicate[k]) %>%
      pull(effort.months)
    
    # Assign value if data exists, otherwise assign NA
    effort.matrix[j, k] <- ifelse(length(val) > 0, val, NA)
    
  }
}

# check the structure of our detection covariates 
str(effort.matrix)

# AS OF 20250130 No active season was included in the detection covariates. A big problem in 
# multi species modelling as all have species-specific active period. 

# det.covs = list(effort = effort.matrix)


# Format site coordinates -------------------------------------------------
coords <- df %>%
  distinct(plotID, .keep_all = TRUE) %>%
  select(x, y)

# Convert coords into matrix
coords = as.matrix(coords)

# Ensure the matrix is numeric
storage.mode(coords) <- "numeric"

# Assign sequential row names (plotID) as characters (e.g., "1", "2", "3", ...)
rownames(coords) <- as.character(1:nrow(coords))

# Rename columns from "x" and "y" to "X" and "Y"
colnames(coords) <- c("X", "Y")

# Verify the structure
str(coords)

# Format beta/occurrence covariates --------------------------------------------------
# Recall the beta covariates and stacked into single raster
bio = rast("RStudio/spOccupancy_NSW_20241219/input/CHELSA_bio_2011-2040_gfdl-esm4_ssp126_V.2.1_EPSG32755.tif")
bio = bio[[c("bio5", "bio12")]] #we will only use bio5 and bio 12 based on pearson correlation analysis in S0

roads_density = rast("RStudio/spOccupancy_NSW_20241219/input/env_roadDensity_EPSG32755.tif")
pd = rast("RStudio/spOccupancy_NSW_20241219/input/aus_pd_2010-2020_1km_UNadj_EPSG32755.tif")
der = rast("RStudio/spOccupancy_NSW_20241219/input/env_der_EPSG32755.tif")

env_stack <- c(bio,roads_density, pd, der)

beta <- terra::extract(env_stack, coords)

str(beta)

# Identify rows with NA values in occ.covs
rows_with_na <- which(rowSums(is.na(beta)) > 0)
print(rows_with_na) # plot 794 contains NA value

# Format detection covariates ---------------------------------------------
# Remove rows in all dataset with NA values 
beta_clean <- beta[-rows_with_na, ]  # Remove rows from occ.covs
effort_clean <- effort.matrix[-rows_with_na, ]
y_clean <- y[ , -rows_with_na, ]                # Remove corresponding rows from y
coords_clean <- coords[-rows_with_na, ]      # Remove corresponding rows from coords

# Pack all things into a list object
y = y_clean
str(y)

det.covs = list(effort = effort_clean)
str(det.covs)

occ.covs = beta_clean
str(occ.covs)

coords = coords_clean
str(coords_clean)

data.sfMsPGOcc <- list(y = y, 
                  occ.covs = occ.covs, 
                  det.covs = det.covs, 
                  coords = coords)
str(data.sfMsPGOcc)

# Test run
out.msom <- spMsPGOcc(occ.formula = ~ scale(bio5) + I(scale(bio5)^2) + scale(bio12) + I(scale(bio12)^2) + # Quadratic for bio12
                        scale(roadLength) +              # Linear for road density
                        scale(pd_mean) +        # Linear for population density
                        scale(der) + I(scale(der)^2),  # Quadratic for depth of regolith, or do we expect expotential relationship?
                      det.formula = ~ scale(effort),
                      data = data.sfMsPGOcc,
                      n.batch = 10, 
                      batch.length = 25, 
                      cov.model = 'exponential', 
                      NNGP = TRUE, 
                      verbose = FALSE) 
summary(out.msom, level = 'community')


save(data.sfMsPGOcc, file = "RStudio/spOccupancy_NSW_20241219/input/data.sfMsPGOcc.RData")






































# Assuming data.msom is your list and y is the array within it
# dimnames(data.sfMsPGOcc$y)[[3]] <- as.character(1:5)

# Now check the structure again
# str(data.sfMsPGOcc)

# save(data.sfMsPGOcc, file = "RStudio/spOccupancy_NSW_20241219/input/data.sfMsPGOcc.RData")

# Test run ---------------------------------------------------------------
# out.msom <- sfMsPGOcc(occ.formula = ~ scale(bio5) + I(scale(bio5)^2) +  # Quadratic for bio5
                        # scale(bio12) + I(scale(bio12)^2) + # Quadratic for bio12
                        # scale(roadLength) +              # Linear for road density
                        # scale(pd_mean) +        # Linear for population density
                        # scale(der) + I(scale(der)^2),  # Quadratic for depth of regolith, or do we expect expotential relationship?
                      # det.formula = ~ scale(det_clean),
                      # data = data.sfMsPGOcc,
                      # n.batch = 10,
                      # n.chains = 4,
                      # n.factors = 1,
                      # batch.length = 25, 
                      # cov.model = 'exponential', 
                      # NNGP = TRUE,
                      # n.neighbors = 15) 

# summary(out.msom)
# names(out.msom)

# summary(out.msom$lambda.samples)

# Data ready to go.