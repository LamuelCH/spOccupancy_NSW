
# For Slurm Script Array  -------------------------------------------------
# Read command-line argument for array index
args <- commandArgs(trailingOnly = TRUE)
thinning_index <- as.integer(args[1])  # SLURM_ARRAY_TASK_ID passed here

# Define thinning factors and corresponding parameters
thinning_factors <- c(1, 10, 100, 250, 500, 1000)
n_thin <- thinning_factors[thinning_index + 1]  # +1 because array starts at 0

# Calculate required n.batch (rounded up to ensure enough samples)
n_batch <- ceiling((1000 * n_thin + 3000) / 25)



# Original script ---------------------------------------------------------
set.seed(123) #so to get the same results


library(spOccupancy)
setwd(".")



# Load the data -----------------------------------------------------------
load("input/list_sfMsPGOcc_2011-2040_gfdl-esm4_ssp126.RData")
str(data.sfMsPGOcc)
# Species codes
sp.names = rownames(data.sfMsPGOcc$y)


# Parametize --------------------------------------------------------------
# check the raw probabilities of occurence from detection-non-dectection data
apply(data.sfMsPGOcc$y, 1, mean, na.rm = TRUE)

# Our modelled community is relatively small (N = 26)
# and the species are all quite similar (mammals)
# We have some really rare species and some pretty common species 
# So we try the model with a small number of latent factors (q in statistical notation)
n.factors <- 2



# Ordering species in detection-nondetection array ------------------------
# careful consideration of the ordering of species can lead to:
# (1) increased interpretability of the factors and factor loadings; 
# (2) faster model convergence; and 
# (3) improved mixing.

# Current species (here family) ordering
sp.names

# Canidae      Dasyuridae         Felidae       Leporidae    Macropodidae 
# 0.334317145     0.042929970     0.114032734     0.081298632     0.345854575 
# Muridae     Peramelidae      Petauridae   Phalangeridae      Potoroidae 
# 0.327609337     0.200697612     0.002414811     0.638314999     0.031660853 
# Pseudocheiridae      Vombatidae 
# 0.030587604     0.293265361 

# Reorder species:
# 1.  Place a common species first
# 2. for the remaning q-1 factors, place species that believed will show different occurrence patterns than the first species and the other species placed before it 
sp.ordered = c("Phalangeridae",
               "Dasyuridae",
               "Canidae",
               "Felidae",
               "Leporidae",
               "Macropodidae",
               "Muridae",
               "Peramelidae",
               "Petauridae",
               "Potoroidae",
               "Pseudocheiridae",
               "Vombatidae")
  
# Create new detection-nondetection data matrix in the new order
y.new = data.sfMsPGOcc$y[sp.ordered, ,]
# Create a new data array
data.ordered = data.sfMsPGOcc 
# Change the data to the new ordered data
data.ordered$y = y.new
str(data.ordered)


# [optional] Specify initial value and prior distributions ---------------------------
# If not specify, the model will use default value, but model is sensitive to initial, 
# if convergence is not achieved, run a single chain with moderate number of iterations until 
# the traceplots looks settling, then extract the estimated mean values for the factor loading matrix 
# and supply these as initial values.

# Pair-wise distance between all sites
dist.data <- dist(data.ordered$coords)
# Exponential correlation model
cov.model <- "exponential"
# Specify all other initial values identical to lfMsPGOcc() from before
# Number of species
N <- nrow(data.ordered$y)
# Initiate all lambda initial values to 0.
lambda.inits <- matrix(0, N, n.factors)
# Set diagonal elements to 1
diag(lambda.inits) <- 1
# Set lower triangular elements to random values from a standard normal dist
lambda.inits[lower.tri(lambda.inits)] <- rnorm(sum(lower.tri(lambda.inits)))
# Check it out
lambda.inits


# Create list of initial values.
inits <- list(alpha.comm = 0, #Community level detection coefficients
              beta.comm = 0, #community level occurrence coefficients
              beta = 0, #species-level occurrence coefficients
              alpha = 0, #specie-level detection coefficients 
              tau.sq.beta = 1, #community-level occurrence variance parameters
              tau.sq.alpha = 1, #community-level detection variance parameters
              lambda = lambda.inits, #species-specific factor loading
              phi = 3/mean(dist.data), #spatial range parameter
              z = apply(data.ordered$y, c(1, 2), max, na.rm = TRUE)) #latent occurrence variables for all species


# Adaptive MCMC Sampler ---------------------------------------------------
accept.rate = 0.43 # leave as default

# Total MCMC samples in each chain = n.batch*batch.length
n.batch = n_batch
batch.length = 25 #default
n.burn = 3000 #default
n.thin = n_thin
n.chains = 4

#   MCMC samples in each of 4 chains 

# Setting Priors ----------------------------------------------------------
# Assume unifor priors for the spatial decay parameter phi 
# We recommend determining the
# bounds of the uniform distribution by computing the smallest distance between sites and the largest distance
# between sites in the observed data set
min.dist <- min(dist.data)
max.dist <- max(dist.data)

priors <- list(beta.comm.normal = list(mean = 0, var = 2.72),
               alpha.comm.normal = list(mean = 0, var = 2.72),
               tau.sq.beta.ig = list(a = 0.1, b = 0.1),
               tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
               phi.unif = list(3 / max.dist, 3 / min.dist))

# Tunning parameters
# we also need to specify an initial value for the tuning parameters for the spatial decay and
# smoothness parameters (if applicable). These values are supplied as input in the form of a list with tags
# phi and nu. The initial tuning value can be any value greater than 0, but we recommend starting the value
# out around 0.5. After some initial runs of the model, if you notice the final acceptance rate of a parameter
# is much larger or smaller than the target acceptance rate (accept.rate), you can then change the initial
# tuning value to get closer to the target rate. Here we set the initial tuning value for phi to 1 after some
# initial exploratory runs of the model.
tuning <- list(phi = 0.5)
# if by the end the value is well below the target acceptance rate, recommend rerunning the model with a larger initial tuning parameter (higher than the final reported value in the displayed output of model progress). Opposite vice versa.

# parallelize
n.omp.threads <- 32 # within-chain parallelization
verbose <- TRUE
n.report <- 50 # Report progress at every 50th batch.




###########################################################################
# Model fitting -----------------------------------------------------------
# Specify formula
occ.formula = ~ scale(bio5) + I(scale(bio5)^2) + scale(bio12) + I(scale(bio12)^2) + 
  scale(roadLength) + scale(pd_mean) + scale(der) + I(scale(der)^2)

det.formula = ~ scale(effort) #assuume no variation in detection probability



# Run Model 
start_time <- Sys.time()
print(paste("Start Time:", start_time))

out.sfMsPGOcc <- sfMsPGOcc(occ.formula = occ.formula,
                           det.formula = det.formula,
                           data = data.ordered,
                           inits = inits,
                           n.batch = n.batch,
                           batch.length = batch.length,
                           accept.rate = 0.43,
                           priors = priors,
                           n.factors = n.factors,
                           cov.model = cov.model,
                           tuning = tuning,
                           n.omp.threads = n.omp.threads,
                           verbose = TRUE,
                           NNGP = TRUE,
                           n.neighbors = 15,
                           n.report = n.report,
                           n.burn = n.burn,
                           n.thin = n.thin,
                           n.chains = n.chains)
                           #k.fold = 2,
                           #k.fold.threads = n.omp.threads,
                           #k.fold.seed = 123

save(out.sfMsPGOcc, 
     file = paste0("models/model_sfMsPGOcc_2011-2040_gfdl-esm4_ssp126_nthin", n_thin, "_nbatch", n_batch, ".RData"))

end_time <- Sys.time()
print(paste("End Time:", end_time))
print(paste("Duration:", end_time - start_time))