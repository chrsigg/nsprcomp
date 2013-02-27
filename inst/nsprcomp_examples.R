library(MASS)

# sparse PCA with four components and five features per component
nsprcomp(Boston, ncomp = 4, k = 5)  

# sparse PCA with different cardinalities per component
nsprcomp(Boston, k = c(1,1,2,2))  

# non-negative sparse PCA using subset removal as the deflation method
nsprcomp(Boston, k = c(1,1,2,2,3), nneg = TRUE, deflation = "remove")  
