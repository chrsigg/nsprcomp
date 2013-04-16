#  Copyright 2012, 2013 Christian Sigg
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

context("nsprcomp.spca")

test_that("cardinality", {
    set.seed(1)
    X <- matrix(rnorm(5*5), 5)
    
    nspc.model <- nsprcomp(X, k = 1)
    card <- colSums(abs(nspc.model$rotation) > 0)
    expect_true(all(card == 1))
    
    nspc.model <- nsprcomp(X, k = 4)
    card <- colSums(abs(nspc.model$rotation) > 0)
    expect_true(all(card <= 4))
    
    nspc.model <- nsprcomp(X, k = 1:5)
    card <- colSums(abs(nspc.model$rotation) > 0)
    expect_true(all(card <= 1:5))
    
    expect_error(nsprcomp(X, ncomp = 3, k = 1:2))
})

test_that("deflation", {
    set.seed(1)
    k <- 4
    d <- 20
    n <- 100
    X = matrix(runif(n*d), n)
    
    nspc <- nsprcomp(X, k = k, rety = TRUE, deflation = "ortho")
    W <- nspc$rotation
    Y <- nspc$y
    for (cc in seq(length(nspc$sdev))) {
        expect_true(sum(abs(Y%*%W[ ,cc])) < 1e-10)
    }
    
    nspc <- nsprcomp(X, k = k, rety = TRUE, deflation = "Schur")
    W <- nspc$rotation
    Y <- nspc$y
    for (cc in seq(length(nspc$sdev))) {
        expect_true(sum(abs(Y%*%W[ ,cc])) < 1e-10)
    }
    
    nspc <- nsprcomp(X, k = k, rety = TRUE, deflation = "remove")
    W <- nspc$rotation
    Y <- nspc$y
    for (cc in seq(length(nspc$sdev))) {
        expect_true(sum(abs(Y%*%W[ ,cc])) < 1e-10)
    }
})

test_that("weighted sparse PCA approximation error", {
    set.seed(1)
    X <- scale(matrix(runif(5*5), 5))
    nspc <- nsprcomp(X, omega = c(1,1,1,1,5), ncomp = 2, k = 3)
    X_hat <- nspc$x%*%t(nspc$rotation)
    
    nrm <- rowSums((X - X_hat)^2)
    expect_true(which(nrm == min(nrm)) == 5)
})

