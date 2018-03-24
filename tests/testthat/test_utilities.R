context("Checking utility functions")

test_that("Diag() works correctly", {

    x1 <- Diag(c("a", "b"))
    x2 <- Diag(c(1, 2, 3, 4))
    x3 <- Diag(c("a", 10, "c"))

    expect_identical(x1, matrix(c("a", "0",
                                  "0", "b"), ncol=2))
    
    expect_identical(x2, matrix(c(1, 0, 0, 0,
                                  0, 2, 0, 0,
                                  0, 0, 3, 0,
                                  0, 0, 0, 4), ncol=4))

    expect_identical(x3, matrix(c("a", "0", "0",
                                  "0", "10", "0",
                                  "0", "0", "c"), ncol=3))

})

test_that("is.pd() works correctly", {

    x1 <- diag(1,4)
    x2 <- matrix(c(1,2,2,1), ncol=2)
    x3 <- diag(1,4)
    x3[1,2] <- x3[2,1] <- NA
    
    expect_true(is.pd(x1))
    expect_false(is.pd(x2))    
    expect_identical(is.pd(x3), NA)
    expect_identical(is.pd(list(x1, x2, x3)), c(TRUE, FALSE, NA)) 
})

test_that("vec2symMat() works correctly", {

    x1 <- vec2symMat(1:10)
    x2 <- vec2symMat(1:10, byrow=TRUE)
    x3 <- vec2symMat(1:10, diag=FALSE)
    x4 <- vec2symMat(1:10, diag=FALSE, byrow=TRUE)

    expect_true(isSymmetric(x1))
    expect_true(isSymmetric(x2))
    expect_true(isSymmetric(x3))
    expect_true(isSymmetric(x4))
    
    expect_identical(x1, matrix(c(1,2,3,4,
                                  2,5,6,7,
                                  3,6,8,9,
                                  4,7,9,10), ncol=4))
    expect_identical(x2, matrix(c(1,2,4,7,
                                  2,3,5,8,
                                  4,5,6,9,
                                  7,8,9,10), ncol=4))
    expect_identical(x3, matrix(c(1,1,2,3,4,
                                  1,1,5,6,7,
                                  2,5,1,8,9,
                                  3,6,8,1,10,
                                  4,7,9,10,1), ncol=5))
    expect_identical(x4, matrix(c(1,1,2,4,7,
                                  1,1,3,5,8,
                                  2,3,1,6,9,
                                  4,5,6,1,10,
                                  7,8,9,10,1), ncol=5))
    
})

test_that("bdiagMat() works correctly", {
    x1 <- bdiagMat( list(matrix(1:4,nrow=2,ncol=2),
                         matrix(5:6,nrow=1,ncol=2)) )

    x2 <- bdiagMat( list(matrix(letters[1:4],nrow=2,ncol=2),
                         matrix(letters[5:6],nrow=1,ncol=2)) )

    expect_identical(x1, matrix(c(1, 3, 0, 0,
                                  2, 4, 0, 0,
                                  0, 0, 5, 6), ncol=4, nrow=3,
                                byrow=TRUE))
    expect_identical(x2, matrix(c("a", "c", "0", "0",
                                  "b", "d", "0", "0",
                                  "0", "0", "e", "f"),
                                ncol=4, nrow=3, byrow=TRUE))
})

test_that("list2matrix() works correctly", {
    x1 <- matrix(c(1,0.5,0.4,0.5,1,0.2,0.4,0.2,1), ncol=3)
    x2 <- matrix(c(1,0.4,NA,0.4,1,NA,NA,NA,NA), ncol=3)  

    expect_identical(list2matrix(list(x1, x2), diag=FALSE),
                     matrix(c(0.5, 0.4, 0.2,
                              0.4, NA, NA),
                            byrow=TRUE, nrow=2, ncol=3,
                            dimnames=list(NULL, c("x2_x1", "x3_x1", "x3_x2"))))

    expect_identical(list2matrix(list(x1, x2), diag=TRUE),
                     matrix(c(1, 0.5, 0.4, 1,  0.2, 1,
                              1, 0.4, NA, 1, NA, NA),
                            byrow=TRUE, nrow=2, ncol=6,
                            dimnames=list(NULL, c("x1_x1", "x2_x1", "x3_x1",
                                                  "x2_x2", "x3_x2", "x3_x3"))))    

    dimnames(x1) <- list( c("x","y","z"), c("x","y","z") )
    dimnames(x2) <- list( c("x","y","z"), c("x","y","z") )

    expect_identical(list2matrix(list(x1, x2), diag=FALSE),
                     matrix(c(0.5, 0.4, 0.2,
                              0.4, NA, NA),
                            byrow=TRUE, nrow=2, ncol=3,
                            dimnames=list(NULL, c("y_x", "z_x", "z_y"))))
    
    expect_identical(list2matrix(list(x1, x2), diag=TRUE),
                     matrix(c(1, 0.5, 0.4, 1,  0.2, 1,
                              1, 0.4, NA, 1, NA, NA),
                            byrow=TRUE, nrow=2, ncol=6,
                            dimnames=list(NULL, c("x_x", "y_x", "z_x",
                                                  "y_y", "z_y", "z_z"))))  

    x3 <- matrix(c(1,0.5,0.5,1), ncol=2)
    x4 <- matrix(c(1,0.4,0.4,1), ncol=2)  

    expect_identical(list2matrix(list(x3, x4), diag=FALSE),
                     matrix(c(0.5,
                              0.4),
                            byrow=TRUE, nrow=2, ncol=1,
                            dimnames=list(NULL, c("x2_x1"))))

    expect_identical(list2matrix(list(x3, x4), diag=TRUE),
                     matrix(c(1, 0.5, 1,
                              1, 0.4, 1),
                            byrow=TRUE, nrow=2, ncol=3,
                            dimnames=list(NULL, c("x1_x1", "x2_x1", "x2_x2"))))    

    dimnames(x3) <- list( c("x","y"), c("x","y") )
    dimnames(x4) <- list( c("x","y"), c("x","y") )

    expect_identical(list2matrix(list(x3, x4), diag=FALSE),
                     matrix(c(0.5,
                              0.4),
                            byrow=TRUE, nrow=2, ncol=1,
                            dimnames=list(NULL, c("y_x"))))

    expect_identical(list2matrix(list(x3, x4), diag=TRUE),
                     matrix(c(1, 0.5, 1,
                              1, 0.4, 1),
                            byrow=TRUE, nrow=2, ncol=3,
                            dimnames=list(NULL, c("x_x", "y_x", "y_y"))))
})
