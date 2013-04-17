library(MASS)

set.seed(1)

# Regular PCA, with tolerance set to return five PCs
prcomp(Boston, tol = 0.36, scale. = TRUE)

# Sparse PCA with different cardinalities per component. The number of components
# is derived from the length of vector k.
nsprcomp(Boston, k = c(13,7,5,5,5), scale. = TRUE)  

# Non-negative sparse PCA with subset removal as the deflation method. Note that
# the feature space is exhausted after four components due to the subset
# removal, and that the parameter k only specifies an upper bound on the 
# cardinality of the PAs.
nsprcomp(Boston, k = c(6,5,5,5,5), nneg = TRUE, deflation = "remove",
         scale. = TRUE)  
