\name{asdev}
\alias{asdev}
\title{Additional Explained Standard Deviation}
\usage{
  asdev(X, W, center = TRUE, scale. = FALSE)
}
\arguments{
  \item{X}{a numeric data matrix with the observations as
  rows}

  \item{W}{a numeric data matrix with the principal axes as
  columns}

  \item{center}{a logical value indicating whether the
  empirical mean of \code{X} should be subtracted.
  Alternatively, a vector of length equal the number of
  columns of \code{X} can be supplied.  The value is passed
  to \code{\link{scale}}.}

  \item{scale.}{a logical value indicating whether the
  columns of \code{X} should be scaled to have unit
  variance before the analysis takes place. The default is
  \code{FALSE} for consistency with \code{prcomp}.
  Alternatively, a vector of length equal the number of
  columns of \code{X} can be supplied.  The value is passed
  to \code{\link{scale}}.}
}
\description{
  \code{asdev} computes the additional standard deviation
  explained by each principal component, taking into
  account the possible non-orthogonality of
  \eqn{\mathbf{W}}{W}.
}
\details{
  The additional standard deviation of a component is
  measured after projecting the corresponding principal
  axis to the ortho-complement space spanned by the
  previous principal axes. This procedure ensures that the
  variance explained by non-orthogonal principal axes is
  not counted multiple times. If the principal axes are
  pairwise orthogonal (e.g. computed using standard PCA),
  the additional standard deviations are identical to the
  standard deviations of the columns of the scores matrix
  \eqn{\mathbf{XW}}{X*W}.
}
\references{
  Mackey, L. (2009) Deflation Methods for Sparse PCA. In
  \emph{Advances in Neural Information Processing Systems}
  (pp. 1017--1024).
}

