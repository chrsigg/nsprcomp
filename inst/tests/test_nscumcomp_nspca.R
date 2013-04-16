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

context("nscumcomp.nspca")

test_that("cardinality", {
    set.seed(1)
    X <- matrix(rnorm(20*10), 20)
    
    nscc <- nscumcomp(X, ncomp = 1, gamma = 1, k = 5, nneg = TRUE)
    card <- sum(abs(nscc$rotation) > 0)
    expect_equal(card, 5)
    
    nscc <- nscumcomp(X, ncomp = 5, gamma = 100, k = 10, nneg = TRUE)
    card <- sum(abs(nscc$rotation) > 0)
    expect_equal(card, 10)
})

test_that("non-negativity", {
    set.seed(1)
    X <- matrix(rnorm(20*10), 20)
    
    nscc <- nscumcomp(X, ncomp = 5, gamma = 1e2, k = 10, nneg = TRUE)
    expect_true(all(nscc$rotation >= 0))
})

test_that("weighted non-negative sparse PCA approximation error", {
    set.seed(1)
    X <- scale(matrix(runif(5*5), 5))
    
    nscc <- nscumcomp(X, omega = c(1,10,1,1,1), ncomp = 1, k = 3, gamma = 1,
                      nneg = TRUE)
    W <- nscc$rotation
    X_hat <- X%*%W%*%solve(t(W)%*%W)%*%t(W)
    nrm <- rowSums((X - X_hat)^2)
    expect_true(which(nrm == min(nrm)) == 2)
    
    nscc <- nscumcomp(X, omega = c(1,1,1,1,5), ncomp = 2, k = 6, gamma = 1,
                      nneg = TRUE)
    W <- nscc$rotation
    X_hat <- X%*%W%*%solve(t(W)%*%W)%*%t(W)
    nrm <- rowSums((X - X_hat)^2)
    expect_true(which(nrm == min(nrm)) == 5)
})