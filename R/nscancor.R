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

#' @export
nscancor <- function (x, y, xcenter = TRUE, ycenter = TRUE, 
                      xscale = FALSE, yscale = FALSE, ncomp = NULL,
                      xpredict = function(Y, x) {
                          return(ginv(Y)%*%x)
                      }, 
                      ypredict = function(X, y) {
                          return(ginv(X)%*%y)
                      }, 
                      tol = NULL, nrestart = 10, em.tol = 1e-3, em.maxiter = 20,
                      verbosity = 0) {
    
    n <- nrow(x)
    dx <- ncol(x)
    dy <- ncol(y)
    
    if (is.null(ncomp)) {
        nc <- min(dx, dy, n)
    } else {
        nc <- ncomp
    }
    
    X <- as.matrix(x)
    X <- scale(X, center = xcenter, scale = xscale)
    xcen <- attr(X, "scaled:center")
    xsc <- attr(X, "scaled:scale")
    Y <- as.matrix(y)
    Y <- scale(Y, center = ycenter, scale = yscale)
    ycen <- attr(Y, "scaled:center")
    ysc <- attr(Y, "scaled:scale")
    if(any(xsc == 0) || any(ysc == 0))
        stop("cannot rescale a constant/zero column to unit variance")
    
    corr <- rep(0, nc)  # additional explained correlation
    W <- matrix(0, dx, nc)  # projection vectors for X
    V <- matrix(0, dy, nc)  # projection vectors for Y
    Qx <- matrix(0, dx, nc)  # orthonormal basis spanned by the canonical variables X%*%W
    Qy <- matrix(0, dy, nc)  # orthonormal basis spanned by the canonical variables Y%*%V
    Xp <- X  # X projected to the orthocomplement space spanned by Qx
    Yp <- Y  # Y projected to the orthocomplement space spanned by Qy
    
    for (cc in seq(nc)) {
        
        corr_opt <- -Inf
        for (rr in seq(nrestart)) {
            
            res <- emcca(X, Xp, Y, Yp, xpredict, ypredict, em.tol, em.maxiter)               
            
            # keep solution with maximum additional explained variance
            if (res$corr > corr_opt) {
                corr_opt <- res$corr
                w_opt <- res$w
                v_opt <- res$v
            }
            if (verbosity > 0) {
                print(paste("component ", cc, ": ",
                            "maximal correlation is ", format(corr_opt, digits = 4),
                            " at random restart ", rr-1, sep = ""))
            }
        }
        corr[cc] <- corr_opt
        w <- w_opt  
        v <- v_opt
        W[ ,cc] <- w
        V[ ,cc] <- v
        
        # update Qx and Qy
        XtXw <- t(X)%*%(X%*%w)
        YtYv <- t(Y)%*%(Y%*%v)
        if (cc > 1) {
            qx <- XtXw - Qx[ , 1:(cc-1)]%*%(t(Qx[ , 1:(cc-1)])%*%XtXw) 
            qy <- YtYv - Qy[ , 1:(cc-1)]%*%(t(Qy[ , 1:(cc-1)])%*%YtYv) 
        } else {
            qx <- XtXw
            qy <- YtYv
        }
        qx <- qx/normv(qx)
        qy <- qy/normv(qy)
        Xp <- Xp - Xp%*%qx%*%t(qx)  # deflate data matrices
        Yp <- Yp - Yp%*%qy%*%t(qy)  
        Qx[ , cc] <- qx
        Qy[ , cc] <- qy
        
        # current additionally explained correlation is below tol threshold
        if (!is.null(tol) && corr[cc] < tol*corr[1]) {
            corr <- corr[1:(cc-1)]
            W <- W[ , 1:(cc-1), drop=FALSE]
            V <- V[ , 1:(cc-1), drop=FALSE]
            break
            
        # at least one data matrix is fully deflated
        } else if (cc < nc && (all(abs(Xp) < 1e-14) || all(abs(Yp) < 1e-14))) {  
            if (!is.null(ncomp)) {
                warning("data matrix is fully deflated, less than 'ncomp' components could be computed")            
            }
            corr <- corr[1:cc]
            W <- W[ , 1:cc, drop=FALSE]
            V <- V[ , 1:cc, drop=FALSE]
            break
        }
    }
    
    rownames(W) <- colnames(X)
    rownames(V) <- colnames(Y)
    return(list(cor = corr, xcoef = W, ycoef = V,
                xcenter = if(is.null(xcen)) rep.int(0, dx) else xcen,  # return value follows cancor interface
                xscale = if(is.null(xsc)) FALSE else xsc,
                ycenter = if(is.null(ycen)) rep.int(0, dy) else ycen,
                yscale = if(is.null(ysc)) FALSE else ysc
    ))
}

emcca <- function(X, Xp, Y, Yp, xpredict, ypredict, em.tol, em.maxiter) {
    
    n <- nrow(X)
    dx <- ncol(X)
    dy <- ncol(Y)
    
    # initialize w and v as random unit vectors, non-negative if the prediction
    # functions return non-negative vectors
    w <- rnorm(dx);
    v <- rnorm(dy);
    if (all(ypredict(Xp, Yp[ , 1]) >= 0)) {
        w <- abs(w)
    }
    if (all(xpredict(Yp, Xp[ , 1]) >= 0)) {
        v <- abs(v)
    }
    w <- w/normv(Xp%*%w)     
    v <- v/normv(Yp%*%v)     
    
    corr_old <- -Inf
    ii <- 0
    while(ii < em.maxiter) {   
        corr <- cor(Xp%*%w, Yp%*%v)
        if (abs(corr - corr_old)/corr < em.tol) {
            break
        }
        corr_old <- corr
        
        w <- ypredict(Xp, Yp%*%v)
        if (all(w == 0))
            stop("w collapsed to the zero vector, try relaxing the constraints")
        w <- w/normv(Xp%*%w)
        
        v <- xpredict(Yp, Xp%*%w)
        if (all(v == 0))
            stop("v collapsed to the zero vector, try relaxing the constraints")
        v <- v/normv(Yp%*%v)
        
        ii <- ii + 1
    }
    if (ii == em.maxiter)
        warning("maximum number of EM iterations reached before convergence")
    
    w <- w/normv(X%*%w)     
    v <- v/normv(Y%*%v) 
    return(list(corr=corr, w=w, v=v))
}