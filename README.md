nsprcomp
========

An R package for non-negative and sparse principal component analysis.

_Principal component analysis_ (PCA) provides a lower dimensional
approximation of high-dimensional data, where the reconstruction error
(measured using Euclidean distance) is minimal. PCA is therefore
useful to succinctly describe a data set which has an "intrinsic"
dimensionality that is low, even though the number of features might
be high. A classical example is the set of face images: even low-resolution
images quickly exceed 1e4 pixels, yet about 50 principal components
(PCs) are sufficient to achieve a good approximation of a face image
(Kirby and Sirovich, 1990).

Although PCA often provides a good approximation with few PCs, each
component is usually a linear combination of all features of the
data. Enforcing sparsity of the principal axes (PA) facilitates
identification of the relevant features and is therefore an
unsupervised feature selection method. In applications where a fixed
penalty is associated with each included feature (e.g. transaction
costs in finance), a small loss in variance of the PC for a large
reduction in cardinality of the PA can lead to an overall more
desirable solution. Enforcing non-negativity of the PAs
renders PCA applicable to domains where only positive influence of
features is deemed appropriate, i.e. the total variance is explained
additively. Non-negative solutions often show some degree of sparsity
already, but a combination of both constraints enables precise control
over the cardinality of the PAs.

This package implements two non-negative and/or sparse PCA algorithms
which are rooted in _expectation-maximization_ (EM) for a
probabilistic generative model of PCA (Sigg and Buhmann, 2008). A
small example from the domain of portfolio optimization which
demonstrates the usage is provided in [this blog
entry](http://sigg-iten.ch/learningbits/2013/05/27/nsprcomp-is-on-cran/), and
a comparison to the `arrayspc` algorithm from the elasticnet package can be found
[here](http://sigg-iten.ch/learningbits/2013/09/15/non-negative-sparse-pca-comparison/).

nsprcomp Algorithm
-------------------------

In each EM iteration of the `nsprcomp`algorithm, a soft thresholding
operator is applied to enforce sparsity of the PA, and projection to
the non-negative orthant is applied to enforce non-negativity.  Once
the EM procedure has converged and the support of the PA has been
identified, the non-zero coefficients are recomputed to maximize the
variance of the PC. Finally, because EM is a local optimizer, random
restarts are employed to (hopefully) avoid bad local minima. The same
iterative procedure (without the non-negativity constraint) was later also
proposed by Journee et al. (2010), although the motivation and
derivation is different.

The `nsprcomp` algorithm computes one PC after the other. Because
constrained PAs no longer correspond to true eigenvectors of the
covariance matrix and are usually not pairwise orthogonal, special
attention needs to be paid when computing more than a single PC. The
algorithm implements the _generalized deflation method_ proposed by
Mackey (2009) to maximize the additional variance of each
component. Given a basis of the space spanned by the previous PAs, the
variance of the PC is maximized after projecting the current PA to the
ortho-complement space of the basis. This procedure maximizes the
additional variance not explained by previous components, and is
identical to standard PCA if no sparsity or non-negativity constraints
are enforced on the PAs.
  
The `nsprcomp` algorithm is suitable for large and high-dimensional
data sets, because it entirely avoids computing the covariance
matrix. It is therefore especially suited to the case where the number of
features exceeds the number of observations.

nscumcomp Algorithm
-------------------------

The `nscumcomp` algorithm jointly computes all PCs, such that the
cumulative variance of all components is maximized. The computation is
again based on EM iterations, but here the maximization step has to be
carried out numerically. Non-negativity of the PAs is achieved by
enforcing a zero lower bound in the L-BFGS-B algorithm used for
performing the maximization step. Sparsity of the PAs is achieved by a
subsequent soft-thresholding of the pseudo-rotation matrix (the matrix
which has the PAs as its columns). In contrast to sequential PCA
algorithms such as `nsprcomp`, only the total number `k` of non-zero
loadings is specified by the user, and the non-zero loadings are
distributed over all PAs by the algorithm.

The algorithm penalizes the divergence of the pseudo-rotation matrix
from ortho-normality using a Lagrange parameter `gamma`. A positive
`gamma` enforces a minimum angle between the PAs, and is sometimes
necessary to avoid PAs collapsing onto each other during the EM
procedure.

The `nscumcomp` algorithm also scales to large and high-dimensional
data sets, as long as the number of components is not too large
(i.e. inverting an nc by nc matrix can still be done, where nc is the
number of components). But due to the numerical maximization of the
M-step it is computationally more involved than the `nsprcomp`
algorithm.

The `nscumcomp` algorithm is currently unpublished.

References
-------------------------

Kirby, M., & Sirovich, L. (1990). Application of the Karhunen-Loeve procedure for the characterization of human faces. _Pattern Analysis and Machine Intelligence, IEEE Transactions on_, 12(1), 103-108.

Sigg, C. D., & Buhmann, J. M. (2008). Expectation-maximization for sparse and non-negative PCA. In _Proceedings of the 25th international conference on Machine learning_ (pp. 960-967).

Journee, M., Nesterov, Y., Richtarik, P., & Sepulchre, R. (2010). Generalized power method for sparse principal component analysis. _The Journal of Machine Learning Research_, 11, 517-553.

Mackey, L. (2009). Deflation methods for sparse pca. _Advances in neural information processing systems_, 21, 1017-1024.