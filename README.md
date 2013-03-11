nsprcomp
========

An R package for non-negative and sparse principal component analysis.

_Principal component analysis_ (PCA) provides a lower dimensional
approximation of high-dimensional data, where the reconstruction error
(measured using Euclidean distance) is minimal. PCA is therefore
useful to succinctly describe a data set which has an "intrinsic"
dimensionality that is low, even though the number of features might
be high. A classical example are images of faces: even low-resolution
images quickly exceed 10e4 pixels, yet about 50 principal components
(PCs) are sufficient to achieve a good approximation of a face image
(Kirby and Sirovich, 1990).

Although PCA often provides a good approximation with few PCs, each
component is usually a linear combination of all features of the
data. Enforcing sparsity of the principal axes (PA) facilitates
identification of the relevant influence factors and is therefore an
unsupervised feature selection method. In applications where a fixed
penalty is associated with each included dimension (e.g. transaction
costs in finance), a small loss in variance for a large reduction in
cardinality can lead to an overall better solution. Enforcing
non-negativity of PAs renders PCA applicable to domains where only
positive influence of features is deemed appropriate, i.e. the total
variance is explained additively by each feature. Non-negative
solutions often show some degree of sparsity already, but a
combination of both constraints enables precise control over the
cardinality of the PAs.

This package implements a non-negative sparse PCA algorithm which is
rooted in _expectation-maximization_ (EM) for a probabilistic
generative model of PCA (Sigg and Buhmann, 2008). In each EM
iteration, a soft thresholding operator is applied to enforce sparsity
of the PA, and projection to the non-negative orthant is applied to
enforce non-negativity.  Once the EM procedure has converged and the
support of the PA has been identified, the non-zero coefficients are
recomputed to maximize the variance of the PC. Finally, because EM is
a local optimizer, random restarts are employed to (hopefully) avoid bad
local minima.

The same approach to sparse PCA (without the non-negativity constraint) 
was later also proposed by Journee et al. (2010), although the derivation 
of the algorithm is different.

Because constrained PAs no longer correspond to true eigenvectors of
the covariance matrix, this package implements three different matrix
deflation techniques to compute more than a single PC. _Orthogonal
projection deflation_ projects the data matrix onto the
orthocomplement space spanned by the principal axes.  _Schur
complement deflation_ projects the data matrix onto the
orthocomplement space spanned by the principal components. These two
deflation methods are presented in Mackey (2009).  Finally, subset
removal simply removes all columns of the data matrix which are
associated with non-zero loadings.
  
The algorithm implemented in this package is suitable for large and
high-dimensional data sets, because it entirely avoids computing the
covariance matrix. It therefore also can handle the case where the
number of features exceeds the number of observations.

References
-------------------------

Kirby, M., & Sirovich, L. (1990). Application of the Karhunen-Loeve procedure for the characterization of human faces. _Pattern Analysis and Machine Intelligence, IEEE Transactions on_, 12(1), 103-108.

Sigg, C. D., & Buhmann, J. M. (2008). Expectation-maximization for sparse and non-negative PCA. In _Proceedings of the 25th international conference on Machine learning_ (pp. 960-967).

Journee, M., Nesterov, Y., Richtarik, P., & Sepulchre, R. (2010). Generalized power method for sparse principal component analysis. _The Journal of Machine Learning Research_, 11, 517-553.

Mackey, L. (2009). Deflation methods for sparse pca. _Advances in neural information processing systems_, 21, 1017-1024.