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

test_that("as.mxMatrix() works correctly", {

    x1 <- matrix(c(1, "2*a", "3@b", 4), ncol=2, nrow=2)
    x1.labels <- c(NA, "a", "b", NA)
    x1.values <- 1:4
    x1.free <- c(FALSE, TRUE, FALSE, FALSE)
    x2 <- mxMatrix(type="Full", nrow=2, ncol=2,
                   free=x1.free, values=x1.values,
                   labels=x1.labels, name="x1")
 
    expect_identical(x2, as.mxMatrix(x1))
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

test_that("lavaan2RAM() works correctly", {
    ## Multiple regression with 2 groups
    model1 <- "y ~ 1 + c(b1, b2)*x1 + c(b3, b4)*x2
               fn1 := b1 + b2
               b3 == b4"
    model2 <- list("1"="y ~ 1 + b1*x1 + b3*x2
                        fn1 := b1 + b2
                        b3 == b4",
                   "2"="y ~ 1 + b2*x1 + b4*x2")       
    RAM1 <- lavaan2RAM(model1, ngroups=2)
    RAM2 <- lapply(model2, lavaan2RAM)
    names(RAM1) <- c("1", "2")
    expect_identical(RAM1, RAM2)

    ## CFA with 2 groups
    model3 <- "f =~ c(a, a)*x1 + c(b1, b2)*x2 + c(c1, c2)*x3 + c(d1, d2)*x4"
    model4 <- list("1"="f =~ a*x1 + b1*x2 + c1*x3 + d1*x4",
                   "2"="f =~ a*x1 + b2*x2 + c2*x3 + d2*x4")
    RAM3 <- lavaan2RAM(model3, ngroups=2)
    RAM4 <- lapply(model4, lavaan2RAM)
    names(RAM3) <- c("1", "2")
    expect_identical(RAM3, RAM4)

    ## Single group multiple regression
    model5 <- "y ~ 1 + b1*x1 + b2*x2"
    RAM5a <- lavaan2RAM(model5)
    ## RAM5b: hard-coded
    RAM5b <- list(A = structure(c("0", "0", "0", "0.1*b1", "0", "0", "0.1*b2", 
                                  "0", "0"), .Dim = c(3L, 3L),
                                .Dimnames = list(c("y", "x1", "x2"),
                                                 c("y", "x1", "x2"))),
                  S = structure(c("0.5*yWITHy", "0", "0", 
                                  "0", "0.5*x1WITHx1", "0*x1WITHx2", "0",
                                  "0*x1WITHx2", "0.5*x2WITHx2"), .Dim = c(3L, 3L),
                                .Dimnames = list(c("y", "x1", "x2"),
                                                 c("y", "x1", "x2"))),
                  F = structure(c(1, 0, 0, 0, 1, 0, 0, 0, 1), .Dim = c(3L, 3L),
                                .Dimnames = list(c("y", "x1", "x2"),
                                                 c("y", "x1", "x2"))), 
                  M = structure(c("0*ymean", "0", "0"), .Dim = c(1L, 3L),
                                .Dimnames = list("1", c("y", "x1", "x2"))))
    expect_identical(RAM5a, RAM5b)

 })

test_that("as.symMatrix() works correctly", {
    A1 <- matrix(c(1:3, "a", "*b", "6*c", 7:9), ncol=3, nrow=3)
    A2 <- matrix(c(1:3, "a", "b", "c", 7:9), ncol=3, nrow=3)    
    A3 <- as.symMatrix(A1)
    expect_identical(A2, A3)

    B1 <- diag(4)
    B2 <- Diag(rep("1", 4))
    B3 <- as.symMatrix(B1)
    expect_identical(B2, B3)
    
    model <- "y ~ b*m + c*x
              m ~ a*x
              x ~~ 1*x
              m ~~ Errm*m
              y ~~ Erry*y
              x ~ meanx*1
              m ~ interceptm*1
              y ~ intercepty*1"
    RAM1 <- lavaan2RAM(model, obs.variables =c("y", "m", "x"))
    RAM2 <- RAM1
    RAM2$A[1, 2] <- "b"
    RAM2$A[1, 3] <- "c"
    RAM2$A[2, 3] <- "a"
    RAM2$S[1, 1] <- "Erry"
    RAM2$S[2, 2] <- "Errm"
    RAM2$M[1, 1] <- "intercepty"
    RAM2$M[1, 2] <- "interceptm"
    RAM2$M[1, 3] <- "meanx"
    RAM2$F[] <- as.character(RAM2$F)
    RAM3 <- as.symMatrix(RAM1)
    expect_identical(RAM2, RAM3)
})


context("Checking functions calculating effect sizes")

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

context("Checking OSMASEM functions")

test_that("Cor2DataFrame() works correctly", {

    ## No moderators
    my.df1 <- Cor2DataFrame(Nohe15A1$data, Nohe15A1$n)
    my.df2 <- Cor2DataFrame(Nohe15A1, append.vars=FALSE)
    expect_equal(my.df1, my.df2, tolerance = .001)
    
    ## Append additional variables
    my.df1$data <- data.frame(my.df1$data,
                              RelW1=Nohe15A1$RelW1,
                              RelW2=Nohe15A1$RelW2,
                              RelS1=Nohe15A1$RelS1,
                              RelS2=Nohe15A1$RelS2,
                              FemalePer=Nohe15A1$FemalePer,
                              Publication=Nohe15A1$Publication,
                              Lag=Nohe15A1$Lag,
                              Country=Nohe15A1$Country,
                              check.names=FALSE)
    my.df2 <- Cor2DataFrame(Nohe15A1, append.vars=TRUE)
    expect_equal(my.df1, my.df2, tolerance = .001)  
})


test_that("checkRAM() works correctly", {
    ## Checking A
    
    ## OK
    A1 <- matrix(c("0", "0", "0",
                   "1*a", "0", "0",
                   "1*b", "1*c", "0"),
                 nrow=3, ncol=3, byrow=TRUE)
    expect_silent(checkRAM(Amatrix=A1))
    expect_silent(checkRAM(Amatrix=as.mxMatrix(A1)))

    ## Diagonals are not zero
    A2 <- matrix(c("0", "0", "0",
                   "1*a", "1", "0",
                   "1*b", "1*c", "0"),
                 nrow=3, ncol=3, byrow=TRUE)
    expect_warning(checkRAM(Amatrix=A2))
    expect_warning(checkRAM(Amatrix=as.mxMatrix(A2)))
    
    A3 <- matrix(c("0", "0", "0",
                   "1*a", "0*d", "0",
                   "1*b", "1*c", "0"),
                 nrow=3, ncol=3, byrow=TRUE)
    expect_warning(checkRAM(Amatrix=A3))
    expect_warning(checkRAM(Amatrix=as.mxMatrix(A3)))

    ## Non-recursive model
    A4 <- matrix(c("0", "0*d", "0",
                   "1*a", "0", "0",
                   "1*b", "1*c", "0"),
                 nrow=3, ncol=3, byrow=TRUE)
    expect_warning(checkRAM(Amatrix=A4))
    expect_warning(checkRAM(Amatrix=as.mxMatrix(A4)))

    ## Checking S
    
    ## OK
    S1 <- matrix(c("1", "0", "0",
                   "0", "0*a", "0*b",
                   "0", "0*b", "0*c"),
                 nrow=3, ncol=3, byrow=TRUE)
    expect_silent(checkRAM(Smatrix=S1, cor.analysis=TRUE))
    expect_silent(checkRAM(Smatrix=as.mxMatrix(S1), cor.analysis=TRUE)) 
    expect_silent(checkRAM(Smatrix=S1, cor.analysis=FALSE))
    expect_silent(checkRAM(Smatrix=as.mxMatrix(S1), cor.analysis=FALSE))
    
    ## Not symmetric in labels
    S2 <- matrix(c("1", "0", "0",
                   "0", "0*a", "0*b1",
                   "0", "0*b2", "0*c"),
                 nrow=3, ncol=3, byrow=TRUE)
    expect_warning(checkRAM(Smatrix=S2, cor.analysis=TRUE))
    expect_warning(checkRAM(Smatrix=as.mxMatrix(S2), cor.analysis=TRUE))
    expect_warning(checkRAM(Smatrix=S2, cor.analysis=FALSE))
    expect_warning(checkRAM(Smatrix=as.mxMatrix(S2), cor.analysis=FALSE))
    
    ## Not symmetric in values
    S3 <- matrix(c("1", "0", "0",
                   "1", "0*a", "0*b",
                   "0", "0*b", "0*c"),
                 nrow=3, ncol=3, byrow=TRUE)
    expect_warning(checkRAM(Smatrix=S3, cor.analysis=TRUE))
    expect_warning(checkRAM(Smatrix=as.mxMatrix(S3), cor.analysis=TRUE))
    expect_warning(checkRAM(Smatrix=S3, cor.analysis=FALSE))
    expect_warning(checkRAM(Smatrix=as.mxMatrix(S3), cor.analysis=FALSE))    

    ## Not symmetric in free parameters
    S4 <- matrix(c("1", "0", "0",
                   "1*d", "0*a", "0*b",
                   "0", "0*b", "0*c"),
                 nrow=3, ncol=3, byrow=TRUE)
    expect_warning(checkRAM(Smatrix=S4, cor.analysis=TRUE))
    expect_warning(checkRAM(Smatrix=as.mxMatrix(S4), cor.analysis=TRUE))
    expect_warning(checkRAM(Smatrix=S4, cor.analysis=FALSE))
    expect_warning(checkRAM(Smatrix=as.mxMatrix(S4), cor.analysis=FALSE))    
 
    ## Checking both A and S
    ## OK
    expect_silent(checkRAM(A=A1, S=S1, cor.analysis=TRUE))

    ## Variance of the IV is a free parameter
    S5 <- matrix(c("1*Err_IV", "0", "0",
                   "0", "0*a", "0*b",
                   "0", "0*b", "0*c"),
                 nrow=3, ncol=3, byrow=TRUE)
    expect_warning(checkRAM(Amatrix=A1, Smatrix=S5, cor.analysis=TRUE))
    ## OK when S is for a covariance structure
    expect_silent(checkRAM(Amatrix=A1, Smatrix=S5, cor.analysis=FALSE))
    
    ## Variance of the IV is not fixed at 1
    S6 <- matrix(c("0", "0", "0",
                   "0", "0*a", "0*b",
                   "0", "0*b", "0*c"),
                 nrow=3, ncol=3, byrow=TRUE)
    expect_warning(checkRAM(Amatrix=A1, Smatrix=S6, cor.analysis=TRUE))
    ## OK when S is for a covariance structure (fewer checking)
    expect_silent(checkRAM(Amatrix=A1, Smatrix=S6, cor.analysis=FALSE))     
})

test_that("create.Tau2() works correctly", {
  ## Symmetric variance component  
  T0 <- create.Tau2(no.var=6, RE.type="Symm", Transform="expLog", 
                    RE.startvalues=0.01)
  vecTau0 <- create.mxMatrix(paste0(log(0.01), "*Tau1_", seq(6)),
                             ncol=1, nrow=6, name="vecTau1")
  Cor0 <- create.mxMatrix(vechs(outer(seq(6), seq(6),
                                      function(x,y) paste0("0*Cor_", x, "_", y))),
                          type="Stand", ncol=6, nrow=6,
                          lbound=-0.99, ubound=0.99, name="Cor")
  expect_identical(T0$vecTau1, vecTau0)
  expect_identical(T0$Cor, Cor0)
  
  ## Diagonal variance component  
  T1 <- create.Tau2(no.var=6, RE.type="Diag", Transform="expLog", 
                    RE.startvalues=0.01)
  vecTau1 <- create.mxMatrix(paste0(log(0.01), "*Tau1_", seq(6)),
                             ncol=1, nrow=6, name="vecTau1")
  Cor1 <- as.mxMatrix(diag(6), name="Cor")
  expect_identical(T1$vecTau1, vecTau1)
  expect_identical(T1$Cor, Cor1)

  ## Zero variance component  
  T2 <- create.Tau2(no.var=6, RE.type="Zero", Transform="expLog", 
                    RE.startvalues=0.01)
  vecTau2 <- create.mxMatrix(rep(log(0),6), type="Full", ncol=1,
                             nrow=6, name="vecTau1")
  Cor2 <- as.mxMatrix(diag(6), name="Cor")
  expect_identical(T2$vecTau1, vecTau2)
  expect_identical(T2$Cor, Cor2)  

  ## User specified diagonal matrix   
  RE.User <- diag(c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE))
  T3 <- create.Tau2(no.var=6, RE.type="User", 
                    Transform="expLog", 
                    RE.User=RE.User, 
                    RE.startvalues=0.01)
  vecTau3 <- paste0(log(0.01), "*Tau1_", seq(6))
  ## Fixed a bug that the values should be log(0) rather than 0 when they are fixed parameters.  
  vecTau3[diag(RE.User)==FALSE] <- log(0)
  vecTau3 <- create.mxMatrix(vecTau3, ncol=1, nrow=6, name="vecTau1")
  Cor3 <- outer(seq(6), seq(6),
                function(x,y) paste0("0*Cor_", x, "_", y))
  Cor3[RE.User==FALSE] <- 0
  Cor3 <- create.mxMatrix(vechs(Cor3), type="Stand", ncol=6, nrow=6,
                          lbound=-0.99, ubound=0.99, name="Cor")
  expect_identical(T3$vecTau1, vecTau3)
  expect_identical(T3$Cor, Cor3)   
  
  ## User specified symmetric matrix   
  RE.User <- diag(c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE))
  RE.User[2,1] <- RE.User[1,2] <- TRUE
  T4 <- create.Tau2(no.var=6, RE.type="User", 
                    Transform="expLog", 
                    RE.User=RE.User, 
                    RE.startvalues=0.01)
  vecTau4 <- paste0(log(0.01), "*Tau1_", seq(6))
  vecTau4[diag(RE.User)==FALSE] <- log(0)
  vecTau4 <- create.mxMatrix(vecTau4, ncol=1, nrow=6, name="vecTau1")
  Cor4 <- outer(seq(6), seq(6),
                function(x,y) paste0("0*Cor_", x, "_", y))
  Cor4[RE.User==FALSE] <- 0
  Cor4 <- create.mxMatrix(vechs(Cor4), type="Stand", ncol=6, nrow=6,
                          lbound=-0.99, ubound=0.99, name="Cor")
  expect_identical(T4$vecTau1, vecTau4)
  expect_identical(T4$Cor, Cor4)  
  
  ## User specified symmetric matrix with errors
  RE.User <- diag(c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE))
  ## Okay
  expect_silent( create.Tau2(no.var=6, RE.type="User", 
                             Transform="expLog", 
                             RE.User=RE.User, 
                             RE.startvalues=0.01) )  
  ## Asymmetric
  RE.User[3,1] <- TRUE
  expect_error( create.Tau2(no.var=6, RE.type="User", 
                           Transform="expLog", 
                           RE.User=RE.User, 
                           RE.startvalues=0.01) )
  ## Estimating covariance but variances are fixed
  RE.User[1,3] <- TRUE 
  expect_error( create.Tau2(no.var=6, RE.type="User", 
                            Transform="expLog", 
                            RE.User=RE.User, 
                            RE.startvalues=0.01) )  
})

context("Checking metaFIML functions")

test_that("metaFIML() works correctly", {

    ## Univariate meta-analysis without AV
    fit1a <- metaFIML(y=r, v=r_v, x=JP_alpha, data=Jaramillo05)

    m1 <- "fy =~ 1*r
           r ~~ data.r_v*r
           fx =~ 1*JP_alpha
           JP_alpha ~~ 0*JP_alpha
           fy ~ Slope1_1*fx
           fy ~~ Tau2_1_1*fy
           fx ~~ CovX1_X1*fx
           fx ~ MeanX1*1
           fy ~ Intercept1*1"

    RAM1 <- lavaan2RAM(m1, obs.variables = c("r", "JP_alpha"), std.lv=FALSE)
    fit1b <- sem(RAM=RAM1, data=Jaramillo05)

    coef1a <- coef(fit1a)
    names1 <- names(coef1a)
    coef1b <- coef(fit1b)[names1] 
    
    ## Equal coefficients within the tolerance
    tolerance <- 1e-3
    expect_equal(coef1a, coef1b, tolerance=tolerance)
    expect_equal(vcov(fit1a), vcov(fit1b)[names1, names1], tolerance=tolerance)
    expect_equal(fit1a$mx.fit$output$Minus2LogLikelihood,
                 fit1b$mx.fit$output$Minus2LogLikelihood)

    ## Univariate meta-analysis with AV
    fit2a <- metaFIML(y=r, v=r_v, x=JP_alpha, av=IDV, data=Jaramillo05)

    m2 <- "fy =~ 1*r
           r ~~ data.r_v*r
           fx =~ 1*JP_alpha
           JP_alpha ~~ 0*JP_alpha
           fy ~ Slope1_1*fx
           fy ~~ Tau2_1_1*fy
           fx ~~ CovX1_X1*fx
           fx ~ MeanX1*1
           fy ~ Intercept1*1

           fz =~ 1*IDV
           IDV ~~ 0*IDV
           fz ~ MeanX2*1
           fz ~~ CovX2_X2*fz + start(818)*fz
           fx ~~ CovX2_X1*fz
           fy ~~ CovX2_Y1*fz"

    RAM2 <- lavaan2RAM(m2, obs.variables = c("r", "JP_alpha", "IDV"), std.lv=FALSE)
    fit2b <- sem(RAM=RAM2, data=Jaramillo05)

    coef2a <- coef(fit2a)
    names2 <- names(coef2a)
    coef2b <- coef(fit2b)[names2] 
    
    ## Equal coefficients within the tolerance
    expect_equal(coef2a, coef2b, tolerance=tolerance)
    ## Remove CovX2_X2 in comparisons as it is too big
    v_fit2a <- vcov(fit2a)[-4, -4]
    v_fit2b <- vcov(fit2b)[names2, names2][-4, -4]
    expect_equal(v_fit2a, v_fit2b, tolerance=tolerance)
    expect_equal(fit2a$mx.fit$output$Minus2LogLikelihood,
                 fit2b$mx.fit$output$Minus2LogLikelihood)

    ## Multivariate meta-analysis without AV
    wvs94a$gnp <- scale(wvs94a$gnp)
    fit3a <- metaFIML(y=cbind(lifesat, lifecon),
                      v=cbind(lifesat_var, inter_cov, lifecon_var),
                      x=gnp, data=wvs94a)

    m3 <- "fy1 =~ 1*lifesat
           lifesat ~~ data.lifesat_var*lifesat
           fy2 =~ 1*lifecon
           lifecon ~~ data.lifecon_var*lifecon
           lifesat ~~ data.inter_cov*lifecon

           fx =~ 1*gnp
           gnp ~~ 0*gnp
           fy1 ~ Slope1_1*fx
           fy2 ~ Slope2_1*fx

           fy1 ~~ Tau2_1_1*fy1
           fy2 ~~ Tau2_2_2*fy2
           fy1 ~~ Tau2_2_1*fy2
           fx ~~ CovX1_X1*fx
           fx ~ MeanX1*1
           fy1 ~ Intercept1*1
           fy2 ~ Intercept2*1"

    RAM3 <- lavaan2RAM(m3, obs.variables = c("lifesat", "lifecon", "gnp"), std.lv=FALSE)
    fit3b <- sem(RAM=RAM3, data=wvs94a)

    coef3a <- coef(fit3a)
    names3 <- names(coef3a)
    coef3b <- coef(fit3b)[names3] 
    
    ## Equal coefficients within the tolerance
    expect_equal(coef3a, coef3b, tolerance=tolerance)
    expect_equal(vcov(fit3a), vcov(fit3b)[names3, names3], tolerance=tolerance)
    expect_equal(fit3a$mx.fit$output$Minus2LogLikelihood,
                 fit3b$mx.fit$output$Minus2LogLikelihood)    
})

test_that("Handling NA in diagonals in tssem1FEM() correctly", {

    var.names <- paste0("x", 1:4) 
    ## All correlations of a variables are NA but the diagonal is 1
    C1 <- matrix(.5, ncol=4, nrow=4)
    diag(C1) <- 1
    C2 <- matrix(.5, ncol=4, nrow=4)
    C2[2, ] <- C2[, 2] <- NA
    diag(C2) <- 1
    C3 <- matrix(.5, ncol=4, nrow=4)
    C3[1, ] <- C3[, 1] <- NA
    diag(C3) <- 1
    dimnames(C1) <- dimnames(C2) <- dimnames(C3) <- list(var.names, var.names)

    C2.NA <- C2
    C2.NA[2,2] <- NA
    C3.NA <- C3
    C3.NA[1,1] <- NA
    
    fit <- tssem1(Cov=list(C1, C2,C3), n=c(50, 50, 50), method="FEM")
    expect_identical(list(C1, C2.NA, C3.NA), fit$data)

    ## Not all correlations are NA. Thus, they cannot be corrected.
    C2[2,3] <- C2[3,2] <- .5
    C3[1,2] <- C3[2,1] <- .5
    expect_error(tssem1(Cov=list(C1, C2,C3), n=c(50, 50, 50), method="FEM"))
})

test_that("Testing new asyCov() correctly", {

    set.seed(123456)
    
    ## Lower tolerance
    tolerance <- 1e-3
    
    new  <- asyCov(x=Becker92$data, n=Becker92$n, acov="individual")
    row.names(new) <- NULL
    old <- asyCovOld(x=Becker92$data, n=Becker92$n, acov="individual")
    expect_equal(new, old, tolerance=tolerance)

    new  <- asyCov(x=Becker92$data, n=Becker92$n, acov="weighted")
    row.names(new) <- NULL
    old <- asyCovOld(x=Becker92$data, n=Becker92$n, acov="weighted")
    expect_equal(new, old, tolerance=tolerance)

    new  <- asyCov(x=Becker92$data, n=Becker92$n, acov="unweighted")
    row.names(new) <- NULL
    old <- asyCovOld(x=Becker92$data, n=Becker92$n, acov="unweighted")
    expect_equal(new, old, tolerance=tolerance)

    new  <- asyCov(x=Becker92$data, n=Becker92$n, acov="individual", as.matrix=FALSE)
    old <- asyCovOld(x=Becker92$data, n=Becker92$n, acov="individual", as.matrix=FALSE)
    expect_equal(new, old, tolerance=tolerance)    

    new  <- asyCov(x=Cheung09$data, n=Cheung09$n, acov="individual")
    row.names(new) <- NULL
    old <- asyCovOld(x=Cheung09$data, n=Cheung09$n, acov="individual")
    expect_equal(new, old, tolerance=tolerance)

    new  <- asyCov(x=Cheung09$data, n=Cheung09$n, acov="weighted")
    row.names(new) <- NULL
    old <- asyCovOld(x=Cheung09$data, n=Cheung09$n, acov="weighted")
    expect_equal(new, old, tolerance=tolerance)

    new  <- asyCov(x=Cheung09$data, n=Cheung09$n, acov="unweighted")
    row.names(new) <- NULL
    old <- asyCovOld(x=Cheung09$data, n=Cheung09$n, acov="unweighted")
    expect_equal(new, old, tolerance=tolerance)

    new  <- asyCov(x=Cheung09$data, n=Cheung09$n, acov="individual", as.matrix=FALSE)
    old <- asyCovOld(x=Cheung09$data, n=Cheung09$n, acov="individual", as.matrix=FALSE)
    expect_equal(new, old, tolerance=tolerance)

    ## Lower tolerance of cor.analysis=F
    tolerance <- 1e-3
    new  <- asyCov(x=Becker92$data, n=Becker92$n, acov="individual", cor.analysis=FALSE)
    row.names(new) <- NULL
    old <- asyCovOld(x=Becker92$data, n=Becker92$n, acov="individual", cor.analysis=FALSE)
    expect_equal(new, old, tolerance=tolerance)

    new  <- asyCov(x=Becker92$data, n=Becker92$n, acov="weighted", cor.analysis=FALSE)
    row.names(new) <- NULL
    old <- asyCovOld(x=Becker92$data, n=Becker92$n, acov="weighted", cor.analysis=FALSE)
    expect_equal(new, old, tolerance=tolerance)

    new  <- asyCov(x=Becker92$data, n=Becker92$n, acov="unweighted", cor.analysis=FALSE)
    row.names(new) <- NULL
    old <- asyCovOld(x=Becker92$data, n=Becker92$n, acov="unweighted", cor.analysis=FALSE)
    expect_equal(new, old, tolerance=tolerance)

    ## Not equal
    ## new  <- asyCov(x=Becker92$data, n=Becker92$n, acov="individual", as.matrix=FALSE, cor.analysis=FALSE)
    ## old <- asyCovOld(x=Becker92$data, n=Becker92$n, acov="individual", as.matrix=FALSE, cor.analysis=FALSE)
    ## expect_equal(new, old, tolerance=tolerance)   
    
})

context("Checking meta function")
test_that("meta() observed statistics is correct", {

    fit <- summary(meta(r, r_v, data=Jaramillo05))
    expect_equal(fit$obsStat, 61)
})
