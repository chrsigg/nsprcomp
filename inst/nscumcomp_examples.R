library(MASS)

# Sparse cumulative PCA with three components and 25 non-zero loadings
# of the rotation matrix. The high orthonormality penalty is necessary to
# avoid co-linear principal axes.
nscumcomp(Boston, ncomp = 3, k = 15, gamma = 1e6)  

# Non-negative sparse cumulative PCA with unequal weighting of samples
nscumcomp(Boston, ncomp = 5, k = 15, nneg = TRUE, omega = 1:nrow(Boston), 
          gamma = 1e8)  