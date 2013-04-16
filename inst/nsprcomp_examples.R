library(MASS)

set.seed(1)

# Regular PCA, with tolerance set to return five PCs
prcomp(Boston, tol = 0.36, scale. = TRUE)

# Sparse PCA with different cardinalities per component. The number of components
# is derived from the length of vector k.
nsprcomp(Boston, k = c(13,7,5,5,5), scale. = TRUE)  

# Non-negative sparse PCA with subset removal as the deflation method
# (note that the parameter k only specifies an upper bound on the PA cardinality)
nsprcomp(Boston, k = c(6,5,5,5,5), nneg = TRUE, deflation = "remove",
         scale. = TRUE)  
