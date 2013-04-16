#  Copyright 2012, 2013 Christian Sigg
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

#' Non-Negative Sparse PCA
#' 
#' Performs a constrained principal component analysis,
#' where non-negativity and/or sparsity constraints are enforced on the principal axes.
#' The results are returned as an object of class \code{nsprcomp}, which inherits from
#' \code{prcomp}.
#' 
#' \code{nsprcomp} computes a principal component (PC) using expectation-maximization
#' iterations, where non-negativity of the loadings is achieved by projection
#' into the non-negative orthant, and sparsity is achieved by soft thresholding.
#' Further PCs are computed by deflating the data matrix and computing the next PC,
#' and so on.
#' 
#' Because constrained principal axes (PAs) no longer correspond to true eigenvectors of the
#' covariance matrix, the package implements three different matrix deflation 
#' techniques to compute more than a single PC. Orthogonal projection deflation 
#' (\code{"ortho"}) projects the data 
#' matrix onto the orthocomplement space spanned by the principal axes. Schur 
#' complement deflation (\code{"Schur"}) projects the data matrix onto the 
#' orthocomplement space spanned by the principal components. Finally, subset 
#' removal (\code{"remove"}) simply removes all columns of the data matrix which
#' are associated with non-zero loadings.
#' 
#' See the references for further details. 
#' 
#' @export nsprcomp
#' @param x a numeric matrix or data frame which provides the data 
#'   for the principal component analysis.
#' @param ... arguments passed to or from other methods.
nsprcomp <- function (x, ...) UseMethod("nsprcomp") 

#' @method nsprcomp default
#' @S3method nsprcomp default
#' @rdname nsprcomp
#' @param retx a logical value indicating whether the principal components
#' (i.e. \code{x} projected into the principal subspace) should be returned.
#' @param ncomp \code{NULL} or an integer indicating the number of principal 
#'   components to be computed. With the default setting, PCs are 
#'   computed until \code{x} is fully deflated, or if \code{tol} is specified,
#'   until the PC magnitude drops below the relative tolerance threshold.
#'   Alternatively, \code{ncomp} can be specified implicitly if \code{k} is given
#'   as a vector.
#' @param omega a vector with as many entries as there are data samples, to
#'   perform weighted PCA (analogous to weighted least-squares). The default is an 
#'   equal weighting of all samples.
#' @param k either a scalar or a vector of length \code{ncomp}, specifying the
#'   upper bounds on the cardinalities of the principal axes.
#' @param nneg a logical value indicating whether the principal axes should be
#'   constrained to the non-negative orthant.
#' @param deflation a character string which is either \code{"ortho"}, \code{"Schur"} 
#'   or \code{"remove"}, indicating the deflation method to be used when computing 
#'   more than a single principal component (see the details section).
#' @param center a logical value indicating whether the empirical mean of \code{x}
#'   should be subtracted. Alternately, a vector of
#'   length equal the number of columns of \code{x} can be supplied.
#'   The value is passed to \code{\link{scale}}.
#' @param scale. a logical value indicating whether the columns of \code{x} should
#'   be scaled to have unit variance before the analysis takes
#'   place. The default is \code{FALSE} for consistency with S, but
#'   in general scaling is advisable.  Alternatively, a vector of length
#'   equal the number of columns of \code{x} can be supplied.  The
#'   value is passed to \code{\link{scale}}.
#' @param tol a value indicating the magnitude below which components
#'   should be omitted. Components are omitted if their
#'   standard deviations are less than or equal to \code{tol} times the
#'   standard deviation of the first component.
#'   With the default \code{NULL} setting, no components
#'   are omitted.  With \code{tol = 0} or \code{tol = sqrt(.Machine$double.eps)}, 
#'   essentially constant components are omitted.
#' @param nrestart an integer indicating the number of random restarts for computing
#'   the principal component via expectation-maximization (EM) iterations. The solution 
#'   achieving maximum standard deviation over all random restarts is kept. A 
#'   value greater than one can help to avoid bad local optima.
#' @param em.tol a lower bound on the minimum relative change of standard deviation, 
#'   used as the stopping criterion for the EM iterations. If the relative change
#'   of PC magnitude changes less than \code{em.tol} between iterations, 
#'   the EM procedure is asssumed to have converged to a local optimum.
#' @param rety a logical value indicating whether the deflated data matrix
#'   should be returned.
#' 
#' @return \code{nsprcomp} returns a list with class \code{(nsprcomp, prcomp)}
#' containing the following elements:
#' \item{sdev}{the standard deviations of the principal components.}
#' \item{rotation}{the matrix of non-negative and/or sparse variable loadings, 
#'   i.e. a matrix whose columns contain the principal axes.}
#' \item{x}{if \code{retx} is \code{TRUE} the principal components, i.e. the 
#'   data projected into the principal subspace (after centering and scaling 
#'   if requested). For the formula method, \code{\link{napredict}()} is applied 
#'   to handle the treatment of values omitted by the \code{na.action}.}
#' \item{center, scale}{the centering and scaling used, or \code{FALSE}.}
#' \item{y}{if \code{rety} is \code{TRUE} the deflated data matrix, for which all
#'   principal axes lie in its null space.}
#'   
#' @note Deflating the data matrix accumulates numerical errors over successive
#' PCs.
#'   
#' @references Sigg, C. D. and Buhmann, J. M. (2008) Expectation-Maximization 
#'   for Sparse and Non-Negative PCA. In \emph{Proceedings of the 25th International 
#'   Conference on Machine Learning} (pp. 960--967).
#' @references Mackey, L. (2009) Deflation Methods for Sparse PCA. In 
#'   \emph{Advances in Neural Information Processing Systems} (pp. 1017--1024).
#'   
#' @seealso \code{\link{prcomp}}, \code{\link{scale}}
#' 
#' @example inst/nsprcomp_examples.R
nsprcomp.default <-
    function(x, retx = TRUE, ncomp = NULL, omega = rep(1, nrow(x)),
             k = ncol(x), nneg = FALSE, deflation = "ortho",
             center = TRUE, scale. = FALSE, tol = NULL,
             nrestart = 5, em.tol = 1e-3, rety = FALSE, ...) 
{        
    d <- ncol(x); n <- nrow(x)
    if (is.null(ncomp)) {
        nc <- min(d,n)
    } else {
        nc <- ncomp
    }
    
    if (length(k) == 1) {
        k <- rep(k, nc)
    } else if (!is.null(ncomp) && length(k) != ncomp) {
        stop("length of k must agree with 'ncomp'")
    } else {
        nc <- length(k)
    }
    
    x <- as.matrix(x)
    x <- scale(x, center = center, scale = scale.)
    cen <- attr(x, "scaled:center")
    sc <- attr(x, "scaled:scale")
    if(any(sc == 0))
        stop("cannot rescale a constant/zero column to unit variance")
    
    X.pc <- matrix(0, n, nc)  # data matrix projected into the PC subspace
    sdev <- rep(0, nc)  # standard deviations of PCs
    W <- matrix(0, d, nc)  # principal axes as columns
    if (deflation == "ortho") {
        Q <- matrix(0, d, nc)  # orthogonalized principal axes as columns
    }
    Y <- x  # data matrix, to be deflated if nc > 1
    for (cc in seq(nc)) {
        sdev.opt <- -Inf
        for (rr in seq(nrestart)) {
            w <- empca(Y, omega, k[cc], nneg, em.tol)
                   
            # variational renormalization
            idx <- abs(w) > 0
            w <- rep(0, d)
            w.sub <- empca(as.matrix(Y[ ,idx]), omega, sum(idx), nneg, em.tol)
            w[idx] <- w.sub                
            
            # keep solution with maximum variance
            x.pc <- x%*%w
            sdev.cur <- sd(as.vector(x.pc))
            if (sdev.cur > sdev.opt) {
                x.pc.opt <- x.pc
                sdev.opt <- sdev.cur
                w.opt <- w
            }
        }
        X.pc[ ,cc] <- x.pc.opt
        sdev[cc] <- sdev.opt
        w <- w.opt  
        W[ ,cc] <- w
        
        if (deflation == "ortho") {
            if (cc > 1) {
                q <- w - Q[ ,1:(cc-1)]%*%(t(Q[ ,1:(cc-1)])%*%w) 
            } else {
                q <- w
            }
            q <- q/as.vector(sqrt(t(q)%*%q))
            Y <- Y - Y%*%q%*%t(q)   
            Q[ ,cc] <- q
        } 
        else if (deflation == "Schur") {
            Yw <- Y%*%w
            Y <- Y - Yw%*%t(Yw)%*%Y/as.vector(t(Yw)%*%Yw)
        }
        else if (deflation == "remove") {
            Y[ ,abs(w) > 0] <- 0
        } else {
            stop("deflation parameter must be one of 'ortho', 'Schur' or 'remove'")
        }
        
        if (!is.null(tol) && sdev[cc] < tol*sdev[1]) {
            X.pc <- X.pc[ ,1:(cc-1)]
            sdev <- sdev[1:(cc-1)]
            W <- W[ ,1:(cc-1)]
            break
        } else if (cc < nc && all(abs(Y) < 1e-14)) {  # data matrix fully deflated
            if (!is.null(ncomp)) {
                warning("data matrix is fully deflated, less than 'ncomp' components could be computed")            
            }
            X.pc <- X.pc[ ,1:cc]
            sdev <- sdev[1:cc]
            W <- W[ ,1:cc]
            break
        }
    }
    
    W <- as.matrix(W)
    dimnames(W) <-
        list(colnames(x), paste0("PC", seq_len(ncol(W))))
    
    r <- list(sdev = sdev, rotation = W,
              center = if(is.null(cen)) FALSE else cen,
              scale = if(is.null(sc)) FALSE else sc)
    if (retx) r$x <- X.pc
    if (rety) r$y <- Y
    class(r) <- c("nsprcomp", "prcomp")
    return(r)
}

empca <- function(X, omega, k, nneg, em.tol) {
    d <- ncol(X); n <- nrow(X)
    
    w <- rnorm(d); 
    if (nneg) {
        w <- abs(w)
    }
    w <- w/sqrt(t(w)%*%w)     
    
    sdev.old <- 0
    repeat {   
        y <- X%*%w
        
        sdev <- sd(as.vector(y))
        if (abs(sdev - sdev.old)/sdev < em.tol) {
            break
        }
        sdev.old <- sdev

        w.star <- as.vector(t(X)%*%(omega*y)/as.vector(t(y)%*%(omega*y)))
        if (nneg) {
            w.star[w.star < 0] <- 0
        }
        if (k<d) {
            # soft thresholding
            s <- sort(abs(w.star), decreasing = TRUE)
            w <- pmax(abs(w.star) - s[k+1], 0)*sign(w.star)  
        } else {
            w <- w.star
        }

        w <- w/sqrt(t(w)%*%w)
    }
    return(w)
}

#' @method nsprcomp formula
#' @S3method nsprcomp formula
#' @rdname nsprcomp
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
nsprcomp.formula <- function (formula, data = NULL, subset, na.action, ...) {
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
    res <- nsprcomp.default(x, ...)
    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl[[1L]] <- as.name("nsprcomp")
    res$call <- cl
    if (!is.null(na.act)) {
        res$na.action <- na.act
        if (!is.null(sc <- res$x))
            res$x <- napredict(na.act, sc)
    }
    res
}
