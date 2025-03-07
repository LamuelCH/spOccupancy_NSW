# Load required libraries
library(tidyverse)
library(ggmap)
library(spOccupancy)
library(terra)
library(sf)
# Setting Project Environment ---------------------------------------------
# Set the path to the root 'data' directory
# setwd("D:/")

options(mc.cores = parallel::detectCores()-1)

# VICTORIA VBA DATA WRANGLING ------------------------------------------------------
  # Specify folder path 
  path_VBA = "data/occ/VBA"
  
  # Create a list of all VBA data inside the folder
  list_VBA = list.files(path = path_VBA, 
                          pattern = "^VBA", 
                          recursive = TRUE, 
                          full.names = TRUE) 
  
  # Create VBA df from the list 
  df_VBA = lapply(list_VBA, function(file) {
      df <- read.csv(file, 
                     header = TRUE,
                     skip = 16, # skip the first 16 rows are metadata                   
                     stringsAsFactors = FALSE)
      df$Scientific.Name <- as.character(df$Scientific.Name) # Convert SpeciesCode to character
      df$Project.ID = as.character(df$Project.ID)
        return(df)
    }) %>% 
    bind_rows(.) # combine all rows together
    
  # column operation
    # Creat new columns essential columns
    # Add year month day columns
      df_VBA = df_VBA %>%
        mutate(
        # Parse the Survey date column into a datetime object
        Date = dmy(Survey.Start.Date), 
    
        # Extract components
        year = year(Date),
        month = month(Date),
        day = day(Date),
        
        #Add a source column and label all values as "VIC"
        Source = "VIC"
        ) 
  
    # Add genus column
      df_VBA$genus <- sub(" .*", "", df_VBA$Scientific.Name)
  
   
  # Re-arranging df and retain essentialy information
    df_VBA = df_VBA %>% 
      select(Source = Source,
             ProjectID = Project.ID,
             Year = year,
             Month = month,
             Day = day,
             Genus = genus,
             ScientificName = Scientific.Name,
             Longitude_GDA94 = Longitude.GDA94,
             Latitude_GDA94 = Latitude.GDA94)

  # Final filtering
    df_VBA = df_VBA %>% 
      distinct() %>%
      filter(!is.na(Longitude_GDA94) & !is.na(Latitude_GDA94)) %>% 
      filter(Year %in% 2012:2021) %>% 
      na.omit() # Remove rows with missing values

    
  # Summarize observationsdistinct() # Summarize observations by genus
    summary_genus_VBA <- df_VBA %>%
      group_by(Genus) %>%
      summarise(Observations = n()) %>%
      arrange(desc(Observations))
    
    # View the summary
    print(summary_genus_VBA, n=32)

# Remove rows where Longitude or Latitude are NA
# df_VBA <- df_VBA %>%
  # filter(!is.na(Longitude.GDA94) & !is.na(Latitude.GDA94))

# Convert df_VBA to an sf object with GDA94 (EPSG: 4283)
# VBA_sf <- df_VBA %>%
  # st_as_sf(coords = c("Longitude.GDA94", "Latitude.GDA94"), crs = 4283)  

# Transform coordinates to Australian Albers Equal Area (EPSG: 3577)
# VBA_sf <- st_transform(VBA_sf, crs = 3577)

# Create the 'plot' column by assigning a unique ID to each unique Easting-Northing combination
# df_VBA <- df_VBA %>%
  # group_by(Easting, Northing) %>%
  # mutate(plot = cur_group_id()) %>%
  # ungroup()



  # rename(
    # genus = genus,      # Rename CommonName to Species
    # plot = plot,
    # x = Easting, # Rename Latitude_GDA94 to Latitude
    # y = Northing # Rename Longitude_GDA94 to Longitude
  # ) %>% 
  # na.omit()


# NSW BIONET DATA WRANGLING -----------------------------------------------
  # Specify folder path 
  path_BioNet <- "data/occ/BioNet"
  
  # Create a list of all BioNet data inside the folder
  list_BioNet <- list.files(path = path_BioNet, 
                          pattern = "^BioNet", 
                          recursive = TRUE, 
                          full.names = TRUE)
  
  # Read all files into a list of DataFrames and ensure SpeciesCode is a character
  df_BioNet <- lapply(list_BioNet, function(file) {
    df <- read.delim(file, 
                     header = TRUE, 
                     sep = "\t", 
                     stringsAsFactors = FALSE)
    df$SpeciesCode <- as.character(df$SpeciesCode) # Convert SpeciesCode to character
    return(df)
  }) %>% 
    bind_rows(.) %>% 
    filter(DatasetName == "Wild Count Fauna") # combine all rows together and only retain WildCount 
  
  # column operation
    # Creat new columns essential columns
    # Add year month day columns
    df_BioNet = df_BioNet %>%
      mutate(
        # Parse the DateLast column into a datetime object
        Date = dmy_hms(DateLast), 
        
        # Extract components
        year = year(Date),
        month = month(Date),
        day = day(Date),
        
        # Add a source column and label all values as "NSW"
        Source = "NSW"
        )
  
    # Add genus column
    df_BioNet$genus <- sub(" .*", "", df_BioNet$ScientificName)
    
   
  # Re-arranging df and retain essentialy information
    df_BioNet = df_BioNet %>% 
      select(Source = Source,
             ProjectID = DatasetName,
             Year = year,
             Month = month,
             Day = day,
             Genus = genus,
             ScientificName = ScientificName,
             Longitude_GDA94 = Longitude_GDA94,
             Latitude_GDA94 = Latitude_GDA94) %>% 
      na.omit()
    
    
  # Final filtering
    df_BioNet = df_BioNet %>% distinct() %>% 
      filter(!is.na(Longitude_GDA94) & !is.na(Latitude_GDA94)) %>% 
      filter(Year %in% 2012:2021) 

    
  # Summarize observations by genus
    summary_genus_BioNet <- df_BioNet %>%
      group_by(Genus) %>%
      summarise(Observations = n()) %>%
      arrange(desc(Observations))
  
    # View the summary
    print(summary_genus_BioNet, n=30)

    
    
# COMBINED OPERATION VIC AND NSW DATA ------------------------------------------------
  df = bind_rows(df_VBA, df_BioNet)

 # Visualize on map
    register_stadiamaps("1794b86f-c7ad-4aab-82b4-23ae3e2f039d", write = TRUE)
    qmplot(Longitude_GDA94, Latitude_GDA94, 
                      data = df,
                      zoom = 8, maptype = "stamen_terrain",
                      darken = c(0.5, "black"), color = Genus)
  
    # Convert PCS to EPSG:3577 
    df <- df %>%  
      filter(!is.na(Longitude_GDA94) & !is.na(Latitude_GDA94)) %>% 
      st_as_sf(coords = c("Longitude_GDA94", "Latitude_GDA94"), 
               crs = 4283)  %>% 
      st_transform( , crs = 3577) %>% 
      mutate(
        Easting_3577 = st_coordinates(.)[,1],  # X (Easting)
        Northing_3577 = st_coordinates(.)[,2]  # Y (Northing)
      )  %>%
      st_drop_geometry()  # Remove spatial features to keep as a data frame  
    
  # Assign sequential plotID for unique Easting-Northing combinations
    df <- df %>%
      arrange(Easting_3577, Northing_3577) %>%  
      mutate(plotID = dense_rank(paste(Easting_3577, Northing_3577, sep = "_"))) %>% 
      mutate(replicate =  ((Year - 2012) %/% 2) + 1) #Map 2012-2021 to 1-10
    
  # Exclude plots surveyed only once
    df_filtered <- df %>%
      group_by(plotID) %>%  # Group by plotID
      filter(n_distinct(replicate) >= 2) %>%  # Keep plots with at least 2 unique replicates
      ungroup()  # Ungroup for further operations
    
  # qmplot(Longitude_GDA94, Latitude_GDA94, 
           # data = df_filtered, 
           # zoom = 8, maptype = "stamen_terrain", 
           # darken = c(0.5, "black"), color = Genus)
  
    # Convert df into long format 
    y.long <- df_filtered %>%
      group_by(plotID, replicate, Genus) %>%
      summarize(occ = ifelse(n() > 0, 1, 0)) %>%
      ungroup() %>%
      glimpse()  
    
    

# FORMAT DETECTION NON-DETECTION DATA ------------------------------------------
# Species codes, adjust accordingly
genus.codes <- sort(unique(y.long$Genus))

# Plot (site) codes.
plot.codes <- sort(unique(y.long$plotID))

# Number of species
N <- length(genus.codes)

# Maximum number of replicates (no. of years) at a site
K <- length(unique(y.long$replicate))  # from 2012 to 2021

# Number of sites
J <- length(unique(y.long$plotID))

# Create array for detection-nondetection data. 
y <- array(NA, dim = c(N, J, K))

# Label the dimensions of y (not necessary, but helpful)
dimnames(y)[[1]] <- genus.codes
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
      for (genus in genus.codes) {
        # Check if this species was observed
        was_observed <- any(sampled_data$Genus == genus)
        
        # Populate the array
        genus_index <- which(dimnames(y)[[1]] == genus)
        site_index <- which(dimnames(y)[[2]] == site)
        y[genus_index, site_index, rep] <- ifelse(was_observed, 1, 0)
      }
    } else {
      # If the plotID was not sampled during this replicate, set all species occurrences to NA
      for (genus in genus.codes) {
        genus_index <- which(dimnames(y)[[1]] == genus)
        site_index <- which(dimnames(y)[[2]] == site)
        y[genus_index, site_index, rep] <- NA
      }
    }
  }
}
str(y)

# Check the total number of observations for each species
apply(y, 1, sum, na.rm = TRUE)



# FORMAT DETECTION COVARIATES ---------------------------------------------
# Survey Effort matrix
  month.df = df_filtered %>% 
    mutate(replicate =  ((Year - 2012) %/% 2) + 1) %>%  #Map 2012-2021 to 1-10
    group_by(plotID, replicate) %>%
    summarize(effort.months = n_distinct(paste(Year, Month, sep = "_"))) %>%  # Count unique (year, month) pairs
    ungroup() %>% 
    glimpse
  
  # convert to numeric format
  month.df <- month.df %>%
    mutate(across(everything(), as.numeric))
  
  # Check the structure of the updated data frame
  str(month.df)
  
  # Detection covariates 1: survey effort
  # Get unique plot IDs and years
  replicate <- 1:K # Years (columns)
  
  # Initialize deployment matrix with NAs
  effort.matrix <- matrix(NA, nrow = J, ncol = K, 
                              dimnames = list(plot.codes, replicate))
  
  # Fill the matrix
  for (j in 1:J) {
    for (k in 1:K) {
      # Get effort for plot j and replicate k
      val <- month.df %>%
        filter(plotID == plot.codes[j], replicate == k) %>%  # Directly use k, not replicate[k]
        pull(effort.months)
      
      effort.matrix[j, k] <- ifelse(length(val) > 0, val, NA)
    }
  }
  
  # check the structure of our detection covariates 
  str(effort.matrix)
  
  # det.covs = list(effort = effort.matrix)
  
  # Check if the NA values aligned between y and det
  na_effort <- is.na(effort.matrix)
  na_y <- is.na(y[1,,]) #Here I just use the first species, but it should be identical for all species. 
  
  identical(na_effort, na_y)
  # Identify positions where the NA structure differs
  diff_positions <- which(na_effort != na_y, arr.ind = TRUE)
  
  # Print the differing positions
  print(diff_positions)
  
  str(na_effort)
  str(na_y)

  
# project matrix
  # Identify the ProjectID for each plotID and replicate
  project.df <- df_filtered %>%
    mutate(replicate = ((Year - 2012) %/% 2) + 1) %>%  # Map 2012-2021 to 1-10
    group_by(plotID, replicate) %>%
    summarize(ProjectID = first(ProjectID)) %>%  # Use the first ProjectID for each plotID and replicate
    ungroup() %>%
    glimpse()
  
  # Convert ProjectID to a factor (if it's not already)
  project.df <- project.df %>%
    mutate(ProjectID = as.factor(ProjectID))
  
  # Initialize project matrix with NAs
  project.matrix <- matrix(NA, nrow = length(plot.codes), ncol = K, 
                           dimnames = list(plot.codes, replicate))
  
  # Fill the matrix
  for (j in 1:length(plot.codes)) {
    for (k in 1:K) {
      # Get ProjectID for plot j and replicate k
      val <- project.df %>%
        filter(plotID == plot.codes[j], replicate == k) %>%
        pull(ProjectID)
      
      # Assign ProjectID to the matrix (if it exists)
      project.matrix[j, k] <- ifelse(length(val) > 0, as.character(val), NA)
    }
  }
  
  # Check the structure of the project matrix
  str(project.matrix)
  head(project.matrix)



# Format site coordinates -------------------------------------------------
coords <- df_filtered %>%
  distinct(plotID, .keep_all = TRUE) %>%
  select(Easting_3577, Northing_3577)

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
bio = rast("input/CHELSA_bio_2011-2040_gfdl-esm4_ssp126_V.2.1_EPSG3577.tif")
bio = bio[[c("bio5", "bio12")]] #we will only use bio5 and bio 12 based on pearson correlation analysis in S0

roads_density = rast("input/env_roadDensity_EPSG3577.tif")
pd = rast("input/aus_pd_2010-2020_1km_UNadj_EPSG3577.tif")
der = rast("input/env_der_EPSG3577.tif")

env_stack <- c(bio,roads_density, pd, der)
# target_crs = "EPSG:3577"
# env_stack_3577 = project(env_stack, target_crs)
# crs(env_stack_3577)



beta <- terra::extract(env_stack, coords)

str(beta)

# Identify rows with NA values in occ.covs
rows_with_na <- which(rowSums(is.na(beta)) > 0)
print(rows_with_na) # plot 794 contains NA value

# Format detection covariates ---------------------------------------------
# Remove rows in all dataset with NA values 
beta_clean <- beta[-rows_with_na, ]  # Remove rows from occ.covs
effort_clean <- effort.matrix[-rows_with_na, ]
project_clean <- project.matrix[-rows_with_na, ]
y_clean <- y[ , -rows_with_na, ]                # Remove corresponding rows from y
coords_clean <- coords[-rows_with_na, ]      # Remove corresponding rows from coords

# Pack all things into a list object
y = y_clean
str(y)

det.covs = list(effort = effort_clean, project = project_clean)
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

save(data.sfMsPGOcc, file = "input/list_sfMsPGOcc_2011-2040_gfdl-esm4_ssp126.RData")

# Test run
test.sfMsPGOcc <- spMsPGOcc(occ.formula = ~ scale(bio5) + I(scale(bio5)^2) + scale(bio12) + I(scale(bio12)^2) + # Quadratic for bio12
                        scale(roadLength) +              # Linear for road density
                        scale(pd_mean) +        # Linear for population density
                        scale(der) + I(scale(der)^2),  # Quadratic for depth of regolith, or do we expect expotential relationship?
                        det.formula = ~ scale(effort), #+ scale(project),
                        data = data.sfMsPGOcc,
                        n.batch = 10, 
                        batch.length = 25, 
                        cov.model = 'exponential', 
                        NNGP = TRUE, 
                        verbose = FALSE) 

summary(test.sfMsPGOcc, level = 'community')

# Data ready for modil fitting






































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