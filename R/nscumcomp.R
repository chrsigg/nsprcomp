#  Copyright 2013 Christian Sigg
#  Copyright 1995-2012 The R Core Team
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

#' Non-Negative Sparse Cumulative PCA
#' 
#' Performs a constrained PCA-like analysis on the given data matrix,
#' where non-negativity and/or sparsity constraints are enforced on the principal axes (PAs).
#' In contrast to regular PCA, where the algorithm sequentially optimizes the PAs
#' such that the variance of each principal component (PC) is maximal, \code{nscumcomp} 
#' jointly optimizes the PAs such that the cumulative variance of all PCs is maximal.
#' Furthermore, only quasi-orthogonality is enforced between PAs, which is especially
#' useful if the PAs are constrained to lie in the non-negative orthant.
#' 
#' \code{nscumcomp} computes all PCs jointly using expectation-maximization (EM)
#' iterations. The maximization step is equivalent to minimizing the objective function
#' 
#' \deqn{\left\Vert \mathbf{X}-\mathbf{Y}\mathbf{W}^{\top}\right\Vert _{F}^{2}+\gamma\left\Vert \mathbf{W}^{\top}\mathbf{W}-\mathbf{I}\right\Vert _{F}^{2}}{norm(X - Y*W^t, "F")^2 + gamma*norm(W^tW - I, "F")^2}
#'  
#' w.r.t. the rotation matrix W, where \code{gamma} is the Lagrange parameter 
#' associated with the ortho-normality 
#' penalty on W. Non-negativity of the loadings is achieved by enforcing a zero lower 
#' bound in the L-BFGS-B algorithm used for the minimization of the objective, 
#' and sparsity is achieved by a subsequent soft 
#' thresholding of the rotation matrix.
#' 
#' @export nscumcomp
#' @param x a numeric matrix or data frame which provides the data 
#'   for the analysis.
#' @param ... arguments passed to or from other methods.
nscumcomp <- function (x, ...) UseMethod("nscumcomp") 

#' @method nscumcomp default
#' @S3method nscumcomp default
#' @rdname nscumcomp
#' @param retx a logical value indicating whether the principal components
#' (i.e. \code{x} projected into the principal subspace) should be returned.
#' @param ncomp an integer indicating the number of principal 
#'   components to be computed. 
#' @param omega a vector with as many entries as there are data samples, to
#'   perform weighted PCA (analogous to weighted least-squares). The default is an 
#'   equal weighting of all samples.
#' @param k an integer specifying an upper bound on the number of non-zero
#'   loadings of the rotation matrix.
#' @param nneg a logical value indicating whether the principal axes should be
#'   constrained to the non-negative orthant.
#' @param gamma a positive number indicating the penalty on the divergence from 
#'   orthonormality of the rotation matrix. A too small value for \code{gamma}
#'   results in co-linear PAs and produces an error.
#' @param center a logical value indicating whether the empirical mean of \code{x}
#'   should be subtracted. Alternately, a vector of
#'   length equal the number of columns of \code{x} can be supplied.
#'   The value is passed to \code{\link{scale}}.
#' @param scale. a logical value indicating whether the columns of \code{x} should
#'   be scaled to have unit variance before the analysis takes
#'   place. The default is \code{FALSE} for consistency with \code{nsprcomp}, but
#'   in general scaling is advisable.  Alternatively, a vector of length
#'   equal the number of columns of \code{x} can be supplied.  The
#'   value is passed to \code{\link{scale}}.
#' @param nrestart an integer indicating the number of random restarts for computing
#'   the principal component via EM iterations. The solution 
#'   achieving maximum cumulative variance over all random restarts is kept. A 
#'   value greater than one can help to avoid bad local optima.
#' @param em.tol a lower bound on the minimum relative change of standard deviation, 
#'   used as the stopping criterion for the EM iterations. If the relative change
#'   of every PC magnitude changes less than \code{em.tol} between iterations, 
#'   the EM procedure is asssumed to have converged to a local optimum.
#' @param verbosity an integer specifying the verbosity level. Greater values
#'   result in more output, the default is to be quiet.
#' 
#' @return \code{nscumcomp} returns a list with class \code{(nsprcomp, prcomp)}
#' containing the following elements:
#' \item{sdev}{the standard deviations of the principal components.}
#' \item{rotation}{the matrix of non-negative and/or sparse variable loadings, 
#'   i.e. a matrix whose columns contain the principal axes.}
#' \item{x}{if \code{retx} is \code{TRUE} the principal components, i.e. the 
#'   data projected into the principal subspace (after centering and scaling 
#'   if requested). For the formula method, \code{\link{napredict}()} is applied 
#'   to handle the treatment of values omitted by the \code{na.action}.}
#' \item{center, scale}{the centering and scaling used, or \code{FALSE}.}
#' 
#' The components are returned in order of decreasing variance for convenience.
#'   
#' @seealso \code{\link{nsprcomp}}, \code{\link{prcomp}}, \code{\link{scale}}
#' 
#' @example inst/nscumcomp_examples.R
nscumcomp.default <-
    function(x, retx = TRUE, ncomp, 
             omega = rep(1, nrow(x)),
             k = length(x), nneg = FALSE, gamma,
             center = TRUE, scale. = FALSE,
             nrestart = 5, em.tol = 1e-2, 
             verbosity = 0, ...) 
{        
    if (missing(ncomp))
        stop("'ncomp' needs to be specified.")
    if (missing(gamma))
        stop("'gamma' needs to be specified.")
    
    d <- ncol(x); n <- nrow(x)
    if (k > d*ncomp) {
        k = d*ncomp
    }
    
    x <- as.matrix(x)
    x <- scale(x, center = center, scale = scale.)
    cen <- attr(x, "scaled:center")
    sc <- attr(x, "scaled:scale")
    if(any(sc == 0))
        stop("cannot rescale a constant/zero column to unit variance")
    
    cumsdev.opt <- 0
    cumsdev <- numeric(nrestart)
    for (rr in seq(nrestart)) {        
        W <- emcumca(x, omega, ncomp, k, nneg, gamma, em.tol, verbosity)
        
        # keep solution with maximum cumulative variance
        x.pc <- x%*%W
        cumsdev[rr] <- sum(apply(x.pc, 2, sd))
        if (cumsdev[rr] > cumsdev.opt) {
            x.pc.opt <- x.pc
            cumsdev.opt <- cumsdev[rr]
            W.opt <- W
        }
        if (verbosity > 0) {
            print(paste("maximal cum. variance is ", format(cumsdev.opt, digits = 4),
                        " at random restart ", rr-1, sep = ""))
        }
    }
    
    # sort components according to variance
    sdev = apply(x.pc.opt, 2, sd)
    if (length(sdev) > 1) {
        srt = sort(sdev, decreasing = TRUE, index.return = TRUE)
        sdev = srt$x
        W.opt <- W.opt[, srt$ix]
        x.pc <- x.pc[, srt$ix]    
    }
    
    dimnames(W.opt) <-
        list(colnames(x), paste0("PC", seq_len(ncol(W.opt))))
    
    r <- list(sdev = sdev, rotation = W.opt,
              center = if(is.null(cen)) FALSE else cen,
              scale = if(is.null(sc)) FALSE else sc)
    if (retx) r$x <- x.pc
    class(r) <- c("nsprcomp", "prcomp")
    return(r)
}

emcumca <- function(X, omega, ncomp, k, nneg, gamma, em.tol, verbosity = 0) {
    d <- ncol(X); n <- nrow(X)
    
    W <- matrix(rnorm(d*ncomp), d)
    if (nneg) {
        W <- abs(W)
    }
    W <- W/t(matrix(rep(sqrt(colSums(W^2)), d), ncomp))
    
    sdev.old <- 0
    repeat {   
        Y <- tryCatch(X%*%W%*%solve(t(W)%*%W), 
                 error = function(e) {
                     stop("Co-linear principal axes, try increasing the orthonormality penalty 'gamma'.")
                 }
        )
        
        sdev <- apply(Y, 2, sd)
        if (all(abs(sdev - sdev.old)/sdev < em.tol)) {
            break
        }
        sdev.old <- sdev
        
        if (verbosity > 1) {
            print(paste("cum. variance: ", format(sum(sdev), digits = 4),
                        " - ortho penalty: ", format(sum((t(W)%*%W - diag(ncomp))^2), digits = 4),
                        sep = ""
                        )
                  )
        }
        
        # objective function
        if (isTRUE(all.equal(omega, rep(1, nrow(X))))) {
            fn <- function(W) {
                dim(W) <- c(d, ncomp)
                return( sum((X-Y%*%t(W))^2) + gamma*sum((t(W)%*%W - diag(ncomp))^2) )
            }
            
        } else {
            fn <- function(W) {
                dim(W) <- c(d, ncomp)
                return( omega%*%rowSums((X-Y%*%t(W))^2) + gamma*sum((t(W)%*%W - diag(ncomp))^2) )
            }
            
        }
        
        # gradient
        if (isTRUE(all.equal(omega, rep(1, nrow(X))))) {
            tYY <- t(Y)%*%Y
            tXY <- t(X)%*%Y
            gr <- function(W) {
                dim(W) <- c(d, ncomp)
                return( 2*(W%*%tYY - tXY) + 4*gamma*(W%*%t(W)%*%W - W) )
            }
        } else {
            OY <- matrix(rep(omega, ncomp), n)*Y
            tYOY <- t(Y)%*%OY
            tXOY <- t(X)%*%OY
            gr <- function(W) {
                dim(W) <- c(d, ncomp)
                return( 2*(W%*%tYOY - tXOY) + 4*gamma*(W%*%t(W)%*%W - W) )
            }
        }
        
        control <- list()
        if (verbosity > 2) {
            control$trace <- 3    
        } else {
            control$trace <- 0
        }
        W.star <- optim(W, fn, gr, method = "L-BFGS-B",
                        lower = if (nneg) 0 else -Inf, control = control)$par
        
        if (k < d*ncomp) {
            
            # soft thresholding
            srt <- sort(abs(W.star), decreasing = TRUE)
            W <- pmax(abs(W.star) - srt[k+1], 0)*sign(W.star)      
        } else {
            W <- W.star
        }
        
        norms <- sqrt(colSums(W^2))
        if (any(norms == 0))
            stop("Principal axis is the zero vector, try increasing 'k' or decreasing 'ncomp'.")
        W <- W/t(matrix(rep(norms, d), ncomp))
    }
    return(W)
}

#' @method nscumcomp formula
#' @S3method nscumcomp formula
#' @rdname nscumcomp
#' @param formula a formula with no response variable, referring only to numeric variables.
#' @param data an optional data frame (or similar: see
#'   \code{\link{model.frame}}) containing the variables in the
#'   formula \code{formula}.  By default the variables are taken from
#'   \code{environment(formula)}.
#' @param subset an optional vector used to select rows (observations) of the
#'   data matrix \code{x}.
#' @param na.action a function which indicates what should happen
#'   when the data contain \code{NA}s.  The default is set by
#'   the \code{na.action} setting of \code{\link{options}}, and is
#'   \code{\link{na.fail}} if that is unset. The \sQuote{factory-fresh}
#'   default is \code{\link{na.omit}}.
nscumcomp.formula <- function (formula, data = NULL, subset, na.action, ...) {
    mt <- terms(formula, data = data)
    if (attr(mt, "response") > 0L)
        stop("response not allowed in formula")
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    mf$... <- NULL
    mf[[1L]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    ## this is not a `standard' model-fitting function,
    ## so no need to consider contrasts or levels
    if (stats:::.check_vars_numeric(mf))
        stop("PCA applies only to numerical variables")
    na.act <- attr(mf, "na.action")
    mt <- attr(mf, "terms")
    attr(mt, "intercept") <- 0L
    x <- model.matrix(mt, mf)
    res <- nscumcomp.default(x, ...)
    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl[[1L]] <- as.name("nscumcomp")
    res$call <- cl
    if (!is.null(na.act)) {
        res$na.action <- na.act
        if (!is.null(sc <- res$x))
            res$x <- napredict(na.act, sc)
    }
    res
}
