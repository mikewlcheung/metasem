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

test_that("smdMTS() works correctly", {
    ## Means
    m <- c(5,NA,7,9,NA)
    ## Sample variances
    v <- c(10,0,11,12,0)
    ## Sample sizes
    n <- c(50,0,52,53,0)

    index <- !is.na(m)
    ## index.y: index on comparisons against the first group
    index.y <- index[-1]    

    ## Comparing against the first group
    x1 <- smdMTS(m=m, v=v, n=n, homogeneity="variance", bias.adjust=TRUE,
                 all.comparisons=FALSE, list.output=TRUE, lavaan.output=FALSE)

    x2 <- smdMTS(m=m[index], v=v[index], n=n[index], homogeneity="variance",
                 bias.adjust=TRUE, all.comparisons=FALSE,
                 list.output=TRUE, lavaan.output=FALSE)    

    ## Check NA in y
    expect_identical(!index.y, unname(is.na(x1$y))) 
    ## Check NA in V
    expect_identical(TRUE, all(is.na(x1$V[!index.y, !index.y])))  
    ## Check the content in y
    expect_identical(unname(x1$y[!is.na(x1$y)]), unname(x2$y))
    ## Check the content in V
    expect_identical(unname(x1$V[!is.na(x1$y), !is.na(x1$y)]), unname(x2$V))

    ## Conducting all comparisons
    x3 <- suppressWarnings( smdMTS(m=m, v=v, n=n, homogeneity="none",
                                   bias.adjust=FALSE, all.comparisons=TRUE,
                                   list.output=TRUE, lavaan.output=FALSE) )

    x4 <- suppressWarnings( smdMTS(m=m[index], v=v[index], n=n[index],
                                   homogeneity="none", bias.adjust=FALSE,
                                   all.comparisons=TRUE, list.output=TRUE,
                                   lavaan.output=FALSE) )
    ## index for y
    k <- length(index)
    index.y <- rep(NA, k*(k-1)/2)
    p <- 1
    for (i in 1:(k-1)) {
        for (j in (i+1):k) {
            index.y[p] <- index[i]&index[j]
            p <- p+1
        }
    }

    ## Check NA in y
    expect_identical(!index.y, unname(is.na(x3$y)))
    ## Check NA in y
    expect_identical(TRUE, all(is.na(x3$V[!index.y, !index.y])))
    ## Check the content in y
    expect_identical(unname(x3$y[!is.na(x3$y)]), unname(x4$y))
    ## Check the content in V
    expect_identical(unname(x3$V[!is.na(x3$y), !is.na(x3$y)]), unname(x4$V))
})

test_that("smdMES() works correctly", {
    ## Sample means of the first group
    m1 <- c(4, NA, 5)
    ## Sample means of the second group
    m2 <- c(5, NA, 6)

    index <- !is.na(m1)
    
    ## Sample covariance matrices
    V1 <- V2 <- matrix(NA, ncol=3, nrow=3)
    V1[index, index] <- c(3,2,2,3)
    V2[index, index] <- c(3.5,2.1,2.1,3.5)

    ## Sample size in Group 1
    n1 <- 20

    ## Sample size in Group 2
    n2 <- 25
    
    ## Assuming homogeneity of covariance matrix
    x1 <- smdMES(m=m1, m2=m2, V1=V1, V2=V2, n1=n1, n2=n2, homogeneity="covariance",
                 bias.adjust=TRUE, list.output=TRUE, lavaan.output=FALSE)
    x2 <- smdMES(m=m1[index], m2=m2[index], V1=V1[index, index],
                 V2=V2[index, index], n1=n1, n2=n2, homogeneity="covariance",
                 bias.adjust=TRUE, list.output=TRUE, lavaan.output=FALSE)
    
    ## Check NA in y
    expect_identical(!index, unname(is.na(x1$y))) 
    ## Check NA in V
    expect_identical(TRUE, all(is.na(x1$V[!index, !index])))    
    ## Check the content in y
    expect_identical(unname(x1$y[!is.na(x1$y)]), unname(x2$y))
    ## Check the content in V
    expect_identical(unname(x1$V[!is.na(x1$y), !is.na(x1$y)]), unname(x2$V))

    ## Without assuming homogeneity of covariance matrix
    x3 <- smdMES(m=m1, m2=m2, V1=V1, V2=V2, n1=n1, n2=n2, homogeneity="none",
                 bias.adjust=FALSE, list.output=TRUE, lavaan.output=FALSE)
    x4 <- smdMES(m=m1[index], m2=m2[index], V1=V1[index, index],
                 V2=V2[index, index], n1=n1, n2=n2, homogeneity="none",
                 bias.adjust=FALSE, list.output=TRUE, lavaan.output=FALSE)
    
    ## Check NA in y
    expect_identical(!index, unname(is.na(x3$y))) 
    ## Check NA in V
    expect_identical(TRUE, all(is.na(x3$V[!index, !index])))    
    ## Check the content in y
    expect_identical(unname(x3$y[!is.na(x3$y)]), unname(x4$y))
    ## Check the content in V
    expect_identical(unname(x3$V[!is.na(x3$y), !is.na(x3$y)]), unname(x4$V))

    
})
