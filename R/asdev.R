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

#' Additional Explained Standard Deviation
#' 
#' \code{asdev} computes the additional
#' standard deviation explained by each principal component, taking into account 
#' the possible non-orthogonality of \eqn{\mathbf{W}}{W}.
#' 
#' The additional standard deviation of a component is measured after projecting the 
#' corresponding principal axis to the ortho-complement space spanned by the 
#' previous principal axes. This procedure ensures that the variance explained
#' by non-orthogonal principal axes is not counted multiple times. If the principal
#' axes are pairwise orthogonal (e.g. computed using standard PCA), the 
#' additional standard deviations are
#' identical to the standard deviations of the columns of the scores matrix
#' \eqn{\mathbf{XW}}{X*W}.
#' 
#' @references Mackey, L. (2009) Deflation Methods for Sparse PCA. In 
#'   \emph{Advances in Neural Information Processing Systems} (pp. 1017--1024).
#'
#' @export
#' @param X a numeric data matrix with the observations as rows
#' @param W a numeric data matrix with the principal axes as columns
#' @param center a logical value indicating whether the empirical mean of \code{X}
#'   should be subtracted. Alternatively, a vector of
#'   length equal the number of columns of \code{X} can be supplied.
#'   The value is passed to \code{\link{scale}}.
#' @param scale. a logical value indicating whether the columns of \code{X} should
#'   be scaled to have unit variance before the analysis takes
#'   place. The default is \code{FALSE} for consistency with \code{prcomp}.
#'   Alternatively, a vector of length
#'   equal the number of columns of \code{X} can be supplied.  The
#'   value is passed to \code{\link{scale}}.
asdev <- function(X, W, center = TRUE, scale. = FALSE) {
    X <- scale(X, center, scale.)
    nc <- ncol(W)
    sdev <- numeric(nc)
    Q <- qr.Q(qr(W))
    Xp <- X
    for (cc in seq_len(nc)) {
        sdev[cc] <- sd(as.vector(Xp%*%W[ , cc]))  # explicit casting to avoid warning in old R versions
        Xp <- Xp - Xp%*%Q[ , cc]%*%t(Q[ , cc])   
    }
    return(sdev)
}