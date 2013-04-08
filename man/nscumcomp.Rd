\name{nscumcomp}
\alias{nscumcomp}
\alias{nscumcomp.default}
\alias{nscumcomp.formula}
\title{Non-Negative Sparse Cumulative PCA}
\usage{
  nscumcomp(x, ...)

  \method{nscumcomp}{default} (x, retx = TRUE, ncomp,
    omega = rep(1, nrow(x)), k = length(x), nneg = FALSE,
    gamma, center = TRUE, scale. = FALSE, nrestart = 5,
    em.tol = 0.001, ...)

  \method{nscumcomp}{formula} (formula, data = NULL,
    subset, na.action, ...)
}
\arguments{
  \item{x}{a numeric matrix or data frame which provides
  the data for the analysis.}

  \item{...}{arguments passed to or from other methods.}

  \item{retx}{a logical value indicating whether the
  principal components (i.e. \code{x} projected into the
  principal subspace) should be returned.}

  \item{ncomp}{an integer indicating the number of
  principal components to be computed.}

  \item{omega}{a vector with as many entries as there are
  data samples, to perform weighted PCA (analogous to
  weighted least-squares). The default is an equal
  weighting of all samples.}

  \item{k}{an integer specifying an upper bound on the
  number of non-zero loadings of the rotation matrix.}

  \item{nneg}{a logical value indicating whether the
  principal axes should be constrained to the non-negative
  orthant.}

  \item{center}{a logical value indicating whether the
  empirical mean of \code{x} should be subtracted.
  Alternately, a vector of length equal the number of
  columns of \code{x} can be supplied.  The value is passed
  to \code{\link{scale}}.}

  \item{scale.}{a logical value indicating whether the
  columns of \code{x} should be scaled to have unit
  variance before the analysis takes place. The default is
  \code{FALSE} for consistency with \code{nsprcomp}, but in
  general scaling is advisable.  Alternatively, a vector of
  length equal the number of columns of \code{x} can be
  supplied.  The value is passed to \code{\link{scale}}.}

  \item{nrestart}{an integer indicating the number of
  random restarts for computing the principal component via
  EM iterations. The solution achieving maximum cumulative
  variance over all random restarts is kept. A value
  greater than one can help to avoid bad local optima.}

  \item{em.tol}{a lower bound on the minimum relative
  change of standard deviation, used as the stopping
  criterion for the EM iterations. If the relative change
  of every PC magnitude changes less than \code{em.tol}
  between iterations, the EM procedure is asssumed to have
  converged to a local optimum.}

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
  elements: \item{sdev}{the standard deviations of the
  principal components.} \item{rotation}{the matrix of
  non-negative and/or sparse variable loadings, i.e. a
  matrix whose columns contain the principal axes.}
  \item{x}{if \code{retx} is \code{TRUE} the principal
  components, i.e. the data projected into the principal
  subspace (after centering and scaling if requested). For
  the formula method, \code{\link{napredict}()} is applied
  to handle the treatment of values omitted by the
  \code{na.action}.} \item{center, scale}{the centering and
  scaling used, or \code{FALSE}.}
}
\description{
  Performs a constrained PCA-like analysis on the given
  data matrix, where non-negativity and/or sparsity
  constraints are enforced on the principal axes (PAs). In
  contrast to regular PCA, where the algorithm sequentially
  optimizes the PAs such that the variance of each
  principal component (PC) is maximal, \code{nscumcomp}
  jointly optimizes the PAs such that the cumulative
  variance of all PCs is maximal. Furthermore, only
  quasi-orthogonality is enforced between PAs, which is
  especially useful if the PAs are constrained to lie in
  the non-negative orthant.
}
\details{
  The results are returned as an object of class
  \code{nsprcomp}, which inherits from \code{prcomp}.

  \code{nscumcomp} computes all PCs jointly using
  expectation-maximization (EM) iterations. The
  maximization step employs the L-BFGS-B algorithm, where
  non-negativity of the loadings is achieved by enforcing a
  zero lower bound. Sparsity is achieved by a subsequent
  soft thresholding of the rotation matrix.
}
\examples{
library(MASS)

# Sparse cumulative PCA with three components and 25 non-zero loadings
# of the rotation matrix. The high orthonormality penalty is necessary to
# avoid co-linear principal axes.
nscumcomp(Boston, ncomp = 3, k = 15, gamma = 1e6)  

# Non-negative sparse cumulative PCA with unequal weighting of samples
nscumcomp(Boston, ncomp = 5, k = 15, nneg = TRUE, omega = 1:nrow(Boston), 
          gamma = 1e8)  
}
\seealso{
  \code{\link{nsprcomp}}, \code{\link{prcomp}},
  \code{\link{scale}}
}

