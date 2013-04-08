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

context("nscumcomp.pca")

test_that("Equivalence of first PC", {
    set.seed(1)
    X <- matrix(rnorm(20*10), 20)
    cc <- nscumcomp(X, ncomp = 1, gamma = 1, em.tol = 1e-10)
    pc <- prcomp(X)
    
    rot_nrm <- norm(abs(cc$rotation) - abs(pc$rotation[ ,1]), "F")
    expect_true(rot_nrm < 1e-3)
    x_nrm <- norm(abs(cc$x) - abs(pc$x[ ,1]), "F")
    expect_true(x_nrm < 1e-3)
    sdev_nrm <- sqrt(sum((cc$sdev - pc$sdev[1])^2))
    expect_true(sdev_nrm < 1e-3)
})

test_that("weighted PCA equivalence to nsprcomp", {
    set.seed(1)
    X <- matrix(runif(20*10), 20)
    
    nscc <- nscumcomp(X, ncomp = 1, gamma = 1, omega = 1:20, em.tol = 1e-5)
    nspc <- nsprcomp(X, ncomp = 1, omega = 1:20, em.tol = 1e-5)
    
    w1 <- nscc$rotation[ ,1]
    w2 <- nspc$rotation[ ,1]
    expect_true(sum(abs(w1)-abs(w2)) < 1e-3)
})
