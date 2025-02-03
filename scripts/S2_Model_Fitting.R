library(spOccupancy)
setwd(".")



# Load the data -----------------------------------------------------------
load("models/data.sfMsPGOcc.RData")
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
# careful consideration of the ordering of species can lead to (1) increased interpretability
# of the factors and factor loadings; (2) faster model convergence; and (3) improved mixing.

# Current species ordering
sp.names
# [1] "Black Rat"                 "Brown Hare"                "Bush Rat"                 
# [4] "Cat"                       "Common Brushtail Possum"   "Common Ringtail Possum"   
# [7] "Fox"                       "Long-nosed Bandicoot"      "Mainland Dusky Antechinus"
# [10] "Mountain Brushtail Possum" "Rabbit"                    "rat"                      
# [13] "Short-eared Possum"        "Southern Greater Glider"   "Spotted-tailed Quoll"     
# [16] "Squirrel Glider"           "Sugar Glider"              "Swamp Rat"                
# [19] "Yellow-footed Antechinus" 

# Reorder species
sp.ordered = c( "Common Brushtail Possum",
                "Southern Greater Glider", #expected to have different occurrence pattern with the possum     
                "Bare-nosed Wombat",
                "Black-striped Wallaby",      
                "Black Rat",
                "Brown Hare",                 
                "Brush-tailed Rock-wallaby", 
                "Bush Rat",                   
                "Cat",                       
                "Common Ringtail Possum",    
                "Eastern Grey Kangaroo",      
                "Fox",                       
                "Long-nosed Bandicoot",       
                "Mainland Dusky Antechinus",  
                "Mountain Brushtail Possum",  
                "Northern long-nosed potoroo",
                "Rabbit",                     
                "rat",                    
                "Red-legged Pademelon",       
                "Red-necked Pademelon",       
                "Short-eared Possum",         
                "Spotted-tailed Quoll",       
                "Squirrel Glider",            
                "Sugar Glider",               
                "Swamp Rat",                  
                "Yellow-footed Antechinus") 
  
  
# Create new detection-nondetection data matrix in the new order
y.new = data.sfMsPGOcc$y[sp.ordered, ,]
# Create a new data array
data.ordered = data.sfMsPGOcc 
# Change the data to the new ordered data
data.ordered$y = y.new
str(data.ordered)


# [optional] Specify initial value and prior distributions ---------------------------
# If not specify, the model will use default value, but model is sensitive to initial, 
# if convergence is not achieved, fun a shingle chain with moderate number of iterations until 
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
inits <- list(alpha.comm = 0,
              beta.comm = 0,
              beta = 0,
              alpha = 0,
              tau.sq.beta = 1,
              tau.sq.alpha = 1,
              lambda = lambda.inits,
              phi = 3/mean(dist.data),
              z = apply(data.ordered$y, c(1, 2), max, na.rm = TRUE))


# Adaptive MCMC Sampler ---------------------------------------------------
accept.rate = 0.43 # leave as it is 

# Total MCMC samples in each chain = n.batch*batch.length
n.batch = 200
batch.length  = 25 # keep this and change n.batch to increase MCMC sampels

n.burn = 3000
n.thin = 2
n.chains = 1


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
tuning <- list(phi = 1)


# parallelize
n.omp.threads <- 1
verbose <- TRUE
n.report <- 50 # Report progress at every 50th batch.




###########################################################################
# Model fitting -----------------------------------------------------------
# Specify formula
occ.formula = ~ scale(bio5) + I(scale(bio5)^2) + scale(bio12) + I(scale(bio12)^2) + 
  scale(roadLength) + scale(pd_mean) + scale(der) + I(scale(der)^2)

det.formula = ~ scale(effort) #assuume no variation in detection probability

# Run Model 
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
                           n.neighbors = 5,
                           n.report = n.report,
                           n.burn = n.burn,
                           n.thin = n.thin,
                           n.chains = n.chains)

names(out.sfMsPGOcc)
summary(out.sfMsPGOcc)
