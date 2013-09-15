#  Copyright 2013 Christian Sigg
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

#' Percentage Explained Additional Variance
#' 
#' \code{peav} computes the percentage explained additional variance  
#' of each principal component, taking into account 
#' the possible non-orthogonality of \eqn{\mathbf{W}}{W}.
#' 
#' The explained additional variance is computed using \code{\link{asdev}} and divided
#' by the total variance of the data to obtain percentages. \code{sum(peav(X, W))}
#' is  one if \eqn{\mathbf{W}}{W} is an ortho-normal basis, e.g. the
#' rotation matrix of a standard PCA.
#' 
#' \code{peav} is useful to compare the solutions of various constrained PCA
#' methods w.r.t. standard PCA.
#' 
#' @export
#' @param X a numeric data matrix with the observations as rows
#' @param W a numeric data matrix with the principal axes as columns
#' @param center a logical value indicating whether the empirical mean of \code{X}
#'   should be subtracted. Alternatively, a vector of
#'   length equal the number of columns of \code{X} can be supplied.
#'   The value is passed to \code{\link{scale}}.
#' @param scale. a logical value indicating whether the columns of \code{x} should
#'   be scaled to have unit variance before the analysis takes
#'   place. The default is \code{FALSE} for consistency with \code{prcomp}.
#'   Alternatively, a vector of length
#'   equal the number of columns of \code{X} can be supplied.  The
#'   value is passed to \code{\link{scale}}.
#'   
#' @note The method produces different results than the "percentage explained 
#'   variance" (\code{pev}) 
#'   computed by the \code{spca} function from the \code{elasticnet} package.
#'   
#' @seealso \code{\link{asdev}},  \code{\link{scale}}
peav <- function(X, W, center = TRUE, scale. = FALSE) {
    variance <- asdev(X, W, center, scale.)^2
    total_variance <- sum(scale(X, center, scale.)^2)/(nrow(X)-1)
    return(variance/total_variance)
}