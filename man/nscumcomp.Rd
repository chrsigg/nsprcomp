\name{nscumcomp}
\alias{nscumcomp}
\alias{nscumcomp.default}
\alias{nscumcomp.formula}
\title{Non-Negative Sparse Cumulative PCA}
\usage{
  nscumcomp(x, ...)

  \method{nscumcomp}{default} (x, retx = TRUE,
    ncomp = min(dim(x)), omega = rep(1, nrow(x)),
    k = d * ncomp, nneg = FALSE, gamma = 0, center = TRUE,
    scale. = FALSE, nrestart = 5, em_tol = 0.001,
    em_maxiter = 20, verbosity = 0, ...)

  \method{nscumcomp}{formula} (formula, data = NULL,
    subset, na.action, ...)
}
\arguments{
  \item{x}{a numeric matrix or data frame which provides
  the data for the analysis.}

  \item{...}{arguments passed to or from other methods.}

  \item{retx}{a logical value indicating whether the
  principal components, i.e. \code{x} projected into the
  principal subspace, should be returned.}

  \item{ncomp}{the number of principal components (PCs) to
  be computed. The default is to compute a full basis for
  \code{x}.}

  \item{omega}{a vector with as many entries as there are
  data samples, to perform weighted PCA (analogous to
  weighted least-squares regression).  The default is an
  equal weighting of all samples.}

  \item{k}{an upper bound on the total number of non-zero
  loadings of the pseudo-rotation matrix. \code{k} is
  increased if necessary to ensure at least one non-zero
  coefficient per principal axis.}

  \item{nneg}{a logical value indicating whether the
  loadings should be non-negative, i.e. the PAs should be
  constrained to the non-negative orthant.}

  \item{gamma}{a non-negative penalty on the divergence
  from orthonormality of the pseudo-rotation matrix. The
  default is not to penalize, but a positive value is
  sometimes necessary to avoid PAs collapsing onto each
  other.}

  \item{center}{a logical value indicating whether the
  empirical mean of \code{x} should be subtracted.
  Alternatively, a vector of length equal the number of
  columns of \code{x} can be supplied.  The value is passed
  to \code{\link{scale}}.}

  \item{scale.}{a logical value indicating whether the
  columns of \code{x} should be scaled to have unit
  variance before the analysis takes place. The default is
  \code{FALSE} for consistency with \code{prcomp}.
  Alternatively, a vector of length equal the number of
  columns of \code{x} can be supplied.  The value is passed
  to \code{\link{scale}}.}

  \item{nrestart}{the number of random restarts for
  computing the pseudo-rotation matrix via
  expectation-maximization (EM) iterations. The solution
  achieving the minimum of the objective function over all
  random restarts is kept. A value greater than one can
  help to avoid bad local minima.}

  \item{em_tol}{If the relative change of the objective is
  less than \code{em_tol} between iterations, the EM
  procedure is asssumed to have converged to a local
  optimum.}

  \item{em_maxiter}{the maximum number of EM iterations to
  be performed. The EM procedure is terminated if either
  the \code{em_tol} or the \code{em_maxiter} criterion is
  satisfied.}

  \item{verbosity}{an integer specifying the verbosity
  level. Greater values result in more output, the default
  is to be quiet.}

  \item{formula}{a formula with no response variable,
  referring only to numeric variables.}

  \item{data}{an optional data frame (or similar: see
  \code{\link{model.frame}}) containing the variables in
  the formula \code{formula}.  By default the variables are
  taken from \code{environment(formula)}.}

  \item{subset}{an optional vector used to select rows
  (observations) of the data matrix \code{x}.}

  \item{na.action}{a function which indicates what should
  happen when the data contain \code{NA}s.  The default is
  set by the \code{na.action} setting of
  \code{\link{options}}, and is \code{\link{na.fail}} if
  that is unset. The \sQuote{factory-fresh} default is
  \code{\link{na.omit}}.}
}
\value{
  \code{nscumcomp} returns a list with class
  \code{(nsprcomp, prcomp)} containing the following
  elements: \item{sdev}{the additional standard deviation
  explained by each component, see \code{\link{asdev}}.}
  \item{rotation}{the matrix of non-negative and/or sparse
  loadings, containing the principal axes as columns.}
  \item{x}{the scores matrix \eqn{\mathbf{XW}}{X*W}
  containing the principal components as columns (after
  centering and scaling if requested). For the formula
  method, \code{\link{napredict}} is applied to handle the
  treatment of values omitted by the \code{na.action}.}
  \item{center, scale}{the centering and scaling used, or
  \code{FALSE}.}

  The components are returned in order of decreasing
  variance for convenience.
}
\description{
  Performs a constrained PCA-like analysis on the given
  data matrix, where non-negativity and/or sparsity
  constraints are enforced on the principal axes (PAs). In
  contrast to regular PCA, which greedily maximises the
  variance of each principal component (PC),
  \code{nscumcomp} jointly optimizes the components such
  that the cumulative variance of all PCs is maximal.
}
\details{
  \code{nscumcomp} computes all PCs jointly using
  expectation-maximization (EM) iterations. The M-step is
  equivalent to minimizing the objective function

  \deqn{\left\Vert
  \mathbf{X}-\mathbf{Z}\mathbf{W}^{\top}\right\Vert
  _{F}^{2}+\gamma\left\Vert
  \mathbf{W}^{\top}\mathbf{W}-\mathbf{I}\right\Vert
  _{F}^{2}}{norm(X - Z*t(W), "F")^2 + gamma*norm(t(W)*W -
  I, "F")^2}

  w.r.t. the pseudo-rotation matrix \eqn{\mathbf{W}}{W},
  where
  \eqn{\mathbf{Z}=\mathbf{X}\mathbf{W}\left(\mathbf{W}^\top\mathbf{W}\right)^{-1}}{Z=X*W*solve(t(W)*W)}
  is the scores matrix modified to account for the
  non-orthogonality of \eqn{\mathbf{W}}{W},
  \eqn{\mathbf{I}}{I} is the identity matrix and
  \code{gamma} is the Lagrange parameter associated with
  the ortho-normality penalty on \eqn{\mathbf{W}}{W}.
  Non-negativity of the loadings is achieved by enforcing a
  zero lower bound in the L-BFGS-B algorithm used for the
  minimization of the objective, and sparsity is achieved
  by a subsequent soft thresholding of \eqn{\mathbf{W}}{W}.
}
\examples{
library(MASS)
set.seed(1)

# Regular PCA, with tolerance set to return five PCs
prcomp(Boston, tol = 0.35, scale. = TRUE)

# Sparse cumulative PCA with five components and 20 non-zero loadings. 
# The orthonormality penalty is set to a value which avoids co-linear principal 
# axes. Note that the non-zero loadings are not distributed uniformly over 
# the components.
nscumcomp(Boston, ncomp = 5, k = 20, gamma = 1e3, scale. = TRUE)  

# Non-negative sparse cumulative PCA
nscumcomp(Boston, ncomp = 5, nneg = TRUE, k = 20, gamma = 1e3, scale. = TRUE)  
}
\seealso{
  \code{\link{asdev}}, \code{\link{peav}},
  \code{\link{nsprcomp}}, \code{\link{scale}}
}

