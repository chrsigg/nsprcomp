\name{peav}
\alias{peav}
\title{Percentage Explained Additional Variance}
\usage{
  peav(X, W, center = TRUE, scale. = FALSE)
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
  columns of \code{x} should be scaled to have unit
  variance before the analysis takes place. The default is
  \code{FALSE} for consistency with \code{prcomp}.
  Alternatively, a vector of length equal the number of
  columns of \code{X} can be supplied.  The value is passed
  to \code{\link{scale}}.}
}
\description{
  \code{peav} computes the percentage explained additional
  variance of each principal component, taking into account
  the possible non-orthogonality of \eqn{\mathbf{W}}{W}.
}
\details{
  The explained additional variance is computed using
  \code{\link{asdev}} and divided by the total variance of
  the data to obtain percentages. \code{sum(peav(X, W))} is
  one if \eqn{\mathbf{W}}{W} is an ortho-normal basis, e.g.
  the rotation matrix of a standard PCA.

  \code{peav} is useful to compare the solutions of various
  constrained PCA methods w.r.t. standard PCA.
}
\note{
  The method produces different results than the
  "percentage explained variance" (\code{pev}) computed by
  the \code{spca} function from the \code{elasticnet}
  package.
}
\seealso{
  \code{\link{asdev}}, \code{\link{scale}}
}

