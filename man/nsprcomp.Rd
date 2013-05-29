\name{nsprcomp}
\alias{nsprcomp}
\alias{nsprcomp.default}
\alias{nsprcomp.formula}
\title{Non-Negative Sparse PCA}
\usage{
  nsprcomp(x, ...)

  \method{nsprcomp}{default} (x, retx = TRUE, ncomp = NULL,
    omega = rep(1, nrow(x)), k = ncol(x), nneg = FALSE,
    deflation = "ortho", center = TRUE, scale. = FALSE,
    tol = NULL, nrestart = 5, em.tol = 0.001, rety = FALSE,
    verbosity = 0, ...)

  \method{nsprcomp}{formula} (formula, data = NULL, subset,
    na.action, ...)
}
\arguments{
  \item{x}{a numeric matrix or data frame which provides
  the data for the principal component analysis.}

  \item{...}{arguments passed to or from other methods.}

  \item{retx}{a logical value indicating whether the
  principal components (i.e. \code{x} projected into the
  principal subspace) should be returned.}

  \item{ncomp}{\code{NULL} or an integer indicating the
  number of principal components to be computed. With the
  default setting, PCs are computed until \code{x} is fully
  deflated, or if \code{tol} is specified, until the PC
  magnitude drops below the relative tolerance threshold.
  Alternatively, \code{ncomp} can be specified implicitly
  if \code{k} is given as a vector.}

  \item{omega}{a vector with as many entries as there are
  data samples, to perform weighted PCA (analogous to
  weighted least-squares). The default is an equal
  weighting of all samples.}

  \item{k}{either a scalar or a vector of length
  \code{ncomp}, specifying the upper bounds on the
  cardinalities of the principal axes.}

  \item{nneg}{a logical value indicating whether the
  principal axes should be constrained to the non-negative
  orthant.}

  \item{deflation}{a character string which is either
  \code{"ortho"}, \code{"Schur"} or \code{"remove"},
  indicating the deflation method to be used when computing
  more than a single principal component (see the details
  section).}

  \item{center}{a logical value indicating whether the
  empirical mean of \code{x} should be subtracted.
  Alternately, a vector of length equal the number of
  columns of \code{x} can be supplied.  The value is passed
  to \code{\link{scale}}.}

  \item{scale.}{a logical value indicating whether the
  columns of \code{x} should be scaled to have unit
  variance before the analysis takes place. The default is
  \code{FALSE} for consistency with S, but in general
  scaling is advisable.  Alternatively, a vector of length
  equal the number of columns of \code{x} can be supplied.
  The value is passed to \code{\link{scale}}.}

  \item{tol}{a value indicating the magnitude below which
  components should be omitted. Components are omitted if
  their standard deviations are less than or equal to
  \code{tol} times the standard deviation of the first
  component.  With the default \code{NULL} setting, no
  components are omitted.  With \code{tol = 0} or \code{tol
  = sqrt(.Machine$double.eps)}, essentially constant
  components are omitted.}

  \item{nrestart}{an integer indicating the number of
  random restarts for computing the principal component via
  expectation-maximization (EM) iterations. The solution
  achieving maximum standard deviation over all random
  restarts is kept. A value greater than one can help to
  avoid bad local optima.}

  \item{em.tol}{a lower bound on the minimum relative
  change of standard deviation, used as the stopping
  criterion for the EM iterations. If the relative change
  of PC magnitude changes less than \code{em.tol} between
  iterations, the EM procedure is asssumed to have
  converged to a local optimum.}

  \item{rety}{a logical value indicating whether the
  deflated data matrix should be returned.}

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
  \code{nsprcomp} returns a list with class
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
  scaling used, or \code{FALSE}.} \item{y}{if \code{rety}
  is \code{TRUE} the deflated data matrix, for which all
  principal axes lie in its null space.}
}
\description{
  Performs a constrained principal component analysis,
  where non-negativity and/or sparsity constraints are
  enforced on the principal axes. The results are returned
  as an object of class \code{nsprcomp}, which inherits
  from \code{prcomp}.
}
\details{
  \code{nsprcomp} computes a principal component (PC) using
  expectation-maximization iterations, where non-negativity
  of the loadings is achieved by projection into the
  non-negative orthant, and sparsity is achieved by soft
  thresholding. Further PCs are computed by deflating the
  data matrix and computing the next PC, and so on.

  Because constrained principal axes (PAs) no longer
  correspond to true eigenvectors of the covariance matrix,
  the package implements three different matrix deflation
  techniques to compute more than a single PC. Orthogonal
  projection deflation (\code{"ortho"}) projects the data
  matrix onto the orthocomplement space spanned by the
  principal axes. Schur complement deflation
  (\code{"Schur"}) projects the data matrix onto the
  orthocomplement space spanned by the principal
  components. Finally, subset removal (\code{"remove"})
  simply removes all columns of the data matrix which are
  associated with non-zero loadings.

  See the references for further details.
}
\note{
  Deflating the data matrix accumulates numerical errors
  over successive PCs.
}
\examples{
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
}
\references{
  Sigg, C. D. and Buhmann, J. M. (2008)
  Expectation-Maximization for Sparse and Non-Negative PCA.
  In \emph{Proceedings of the 25th International Conference
  on Machine Learning} (pp. 960--967).

  Mackey, L. (2009) Deflation Methods for Sparse PCA. In
  \emph{Advances in Neural Information Processing Systems}
  (pp. 1017--1024).
}
\seealso{
  \code{\link{prcomp}}, \code{\link{scale}}
}

