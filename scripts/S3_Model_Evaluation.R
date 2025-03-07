library(spOccupancy)


# Load model --------------------------------------------------------------
load("models/model_sfMsPGOcc_2011-2040_gfdl-esm4_ssp126_nthin100_nbatch12120.RData")



# inspect the Rhat and ESS values
summary(out.sfMsPGOcc)

# Beta occurrence covariates ------------------------------------------------------------------
# visualize community occurrence covariates convergence
plot(out.sfMsPGOcc, "beta.comm", density =FALSE)

# [Alternative] convergence of individual species
# plot(out.sfMsPGOcc, "beta", density =FALSE)


# Lamba latent loading ----------------------------------------------------
summary(out.sfMsPGOcc$lambda.samples)

# inspect latent factor loading
out.sfMsPGOcc$ESS$lambda #ESS
# Phalangeridae-1      Dasyuridae-1         Canidae-1 
# 0.000000          6.121878          5.062315 
# Felidae-1       Leporidae-1    Macropodidae-1 
# 8.729760          4.512669          7.849480 
# Muridae-1     Peramelidae-1      Petauridae-1 
# 4.002509          4.040657         28.171606 
# Potoroidae-1 Pseudocheiridae-1      Vombatidae-1 
# 7.396990         33.825440          2.579985 s
# Phalangeridae-2      Dasyuridae-2         Canidae-2 
# 0.000000          0.000000          5.378794
# Felidae-2       Leporidae-2    Macropodidae-2 
# 62.325358          5.515756          4.870279 
# Muridae-2     Peramelidae-2      Petauridae-2 
# 3.975408          6.067882         61.970080 
# Potoroidae-2 Pseudocheiridae-2      Vombatidae-2 
# 5.809854         23.895240          5.342259 
out.sfMsPGOcc$rhat$lambda.lower.tri #Rhat
# [1] 12.192902 13.818766  7.474647 10.476819  8.445917
# [6] 13.091519 17.318131  2.410477 11.987074  2.908469
# [11] 20.174667 11.760799  1.479879 11.517220 10.461248
# [16] 13.738965 16.902488  1.389088 11.269297  4.536947
# [21] 10.581946

# Generally when fitting spatial occupancy models, the detection parameters will 
# mix and converge faster than the occupancy parameters (assuming an adequate 
# amount of replication to separate occupancy from detection, which may not be 
# true with a very small number of replicates). In particular, the spatial decay 
# parameters and the occupancy intercepts may show slow mixing and/or convergence, 
# which I discuss more in depth below as to what to do in this case.

plot(out.sfMsPGOcc, 'lambda', density = FALSE) #Visual
# nothing converged atm





# Theta spatial decay parameter -------------------------------------------
# Spatial decay parameter
plot(out.sfMsPGOcc, 'theta', density = FALSE)



summary(out.sfMsPGOcc)

summary(out.sfMsPGOcc$lambda.samples)

# Takes a few seconds to run. 
ppc.occ.out <- ppcOcc(out.sfMsPGOcc, 'freeman-tukey', group = 2)

# Species-specific intercepts
plot(out.sfMsPGOcc$beta.samples[, 1:7], density = FALSE)

plot(out.sfMsPGOcc, "theta", density = FALSE)
# Calculate Bayesian p-values
summary(ppc.occ.out)
