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
    ## Parameters not explicitly shared with c(l,l) must default to being
    ## FREE and DISTINCT per group (lavaan's own default), so a model
    ## written out separately per group must repeat that with explicit,
    ## per-group-suffixed labels to be equivalent to the combined form.
    ## x1/x2's own means are also free-by-default (lavaan2RAM()'s
    ## documented convention: observed variables' means are free unless
    ## explicitly fixed), so the hand-written per-group text needs its own
    ## explicit, per-group-suffixed mean lines for x1/x2 too, matching the
    ## auto-generated "0*x1mean_g1" label lavaan2RAM(model1, ngroups=2)
    ## produces for them.
    model2 <- list("1"="y ~ ymean_g1*1 + b1*x1 + b3*x2
                        y ~~ yWITHy_g1*y
                        x1 ~~ x1WITHx1_g1*x1 + x1WITHx2_g1*x2
                        x2 ~~ x2WITHx2_g1*x2
                        x1 ~ x1mean_g1*1
                        x2 ~ x2mean_g1*1
                        fn1 := b1 + b2
                        b3 == b4",
                   "2"="y ~ ymean_g2*1 + b2*x1 + b4*x2
                        y ~~ yWITHy_g2*y
                        x1 ~~ x1WITHx1_g2*x1 + x1WITHx2_g2*x2
                        x2 ~~ x2WITHx2_g2*x2
                        x1 ~ x1mean_g2*1
                        x2 ~ x2mean_g2*1")
    RAM1 <- lavaan2RAM(model1, ngroups=2)
    RAM2 <- lapply(model2, lavaan2RAM)
    ## RAM1 is tagged "RAM_multigroup" by lavaan2RAM(ngroups=2); a manually
    ## assembled list of per-group RAM objects (as RAM2 is here) is not
    ## tagged automatically and must be tagged by hand to be treated the
    ## same way by sem()'s class-based dispatch.
    class(RAM2) <- c("RAM_multigroup", class(RAM2))
    names(RAM1) <- c("1", "2")
    expect_identical(RAM1, RAM2)
    expect_s3_class(RAM1, "RAM_multigroup")

    ## CFA with 2 groups
    model3 <- "f =~ c(a, a)*x1 + c(b1, b2)*x2 + c(c1, c2)*x3 + c(d1, d2)*x4"
    model4 <- list("1"="f =~ a*x1 + b1*x2 + c1*x3 + d1*x4
                        x1 ~~ x1WITHx1_g1*x1
                        x2 ~~ x2WITHx2_g1*x2
                        x3 ~~ x3WITHx3_g1*x3
                        x4 ~~ x4WITHx4_g1*x4
                        x1 ~ x1mean_g1*1
                        x2 ~ x2mean_g1*1
                        x3 ~ x3mean_g1*1
                        x4 ~ x4mean_g1*1",
                   "2"="f =~ a*x1 + b2*x2 + c2*x3 + d2*x4
                        x1 ~~ x1WITHx1_g2*x1
                        x2 ~~ x2WITHx2_g2*x2
                        x3 ~~ x3WITHx3_g2*x3
                        x4 ~~ x4WITHx4_g2*x4
                        x1 ~ x1mean_g2*1
                        x2 ~ x2mean_g2*1
                        x3 ~ x3mean_g2*1
                        x4 ~ x4mean_g2*1")
    RAM3 <- lavaan2RAM(model3, ngroups=2)
    RAM4 <- lapply(model4, lavaan2RAM)
    class(RAM4) <- c("RAM_multigroup", class(RAM4))
    names(RAM3) <- c("1", "2")
    expect_identical(RAM3, RAM4)
    expect_s3_class(RAM3, "RAM_multigroup")

    ## Auto-generated labels must be distinct per group by default (the bug
    ## this fix addresses): only user-requested equalities (c(b,b)*x) are
    ## shared across groups.
    modelShared <- "y ~ c(a1, a2)*1 + c(b, b)*x"
    RAMshared <- lavaan2RAM(modelShared, ngroups=2)
    ## User-shared slope stays identical across groups ...
    expect_identical(RAMshared[["1"]]$A["y","x"], RAMshared[["2"]]$A["y","x"])
    ## ... but the auto-generated residual variance of y does not.
    expect_false(identical(RAMshared[["1"]]$S["y","y"], RAMshared[["2"]]$S["y","y"]))
    ## Single-group output is unaffected (no suffix, no new class).
    expect_false(inherits(lavaan2RAM("y ~ x"), "RAM_multigroup"))

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
                  ## x1/x2's own means are free by default too (not just
                  ## y's), like every other observed variable.
                  M = structure(c("0*ymean", "0*x1mean", "0*x2mean"), .Dim = c(1L, 3L),
                                .Dimnames = list("1", c("y", "x1", "x2"))))
    expect_identical(RAM5a, RAM5b)

 })

test_that("lavaan2RAM() frees an observed predictor's own mean by default, without breaking the single-indicator-marker identification pattern", {
    ## Regression test for a high-severity bug (external review): observed
    ## variables' means are documented as free by default, but lavaanify's
    ## own auto-generated, still-fixed-at-0 "~1" row for an unlabelled
    ## predictor used to silently override that default with a hard 0.
    ## Invisible whenever the predictor's population mean happens to be
    ## ~0 (as in most of this file's other tests, which is why this went
    ## unnoticed for so long), but for a genuinely off-centre predictor
    ## the matching "~~" (variance) parameter then absorbs the UNCENTERED
    ## second moment (mean^2 + variance) instead of the actual variance,
    ## and -2LL diverges sharply from an equivalent direct lavaan fit.
    set.seed(31)
    nA <- 300; nB <- 300
    xA <- rnorm(nA, 10, 1); yA <- 2 + 0.5*xA + rnorm(nA)
    xB <- rnorm(nB, -8, 1); yB <- -1 + 0.3*xB + rnorm(nB)
    df <- data.frame(y=c(yA,yB), x=c(xA,xB), g=rep(c("A","B"), c(nA,nB)))

    model <- "y ~ c(b1,b2)*x
              y ~ c(a1,a2)*1
              y ~~ c(e1,e2)*y
              x ~~ c(vx1,vx2)*x"
    RAM <- lavaan2RAM(model, ngroups=2)
    fit <- sem(RAM=RAM, data=df, group="g")
    fit.lav <- lavaan::sem(model, data=df, group="g", fixed.x=FALSE)

    cf <- coef(fit); lav.cf <- coef(fit.lav)
    expect_equal(unname(cf["vx1"]), unname(lav.cf["vx1"]), tolerance=1e-3)
    expect_equal(unname(cf["vx2"]), unname(lav.cf["vx2"]), tolerance=1e-3)
    expect_equal(fit$mx.fit@output$Minus2LogLikelihood,
                 -2*as.numeric(lavaan::fitMeasures(fit.lav, "logl")), tolerance=1e-2)
    ## The bug specifically produced an inflated, mean^2-contaminated
    ## "variance": guard against a future regression collapsing back to
    ## roughly mean(x)^2 (~100 and ~64 here) rather than the true ~1.
    expect_true(unname(cf["vx1"]) < 5)
    expect_true(unname(cf["vx2"]) < 5)

    ## The single-indicator-marker pattern used throughout this package's
    ## own meta-analysis functions ("f =~ 1*yi", with the free mean
    ## modelled entirely through "f ~ mu*1", never through yi's own "~1")
    ## must stay well-identified: freeing yi's own mean too, on top of
    ## f's, would make the two literally redundant (their SUM, not each
    ## individually, is identified) -- confirmed this is a real risk, not
    ## hypothetical: real lavaan::sem() has the exact same non-
    ## identification for this exact pattern (it emits its own "model may
    ## not be identified" warning), so lavaan2RAM() deliberately does NOT
    ## follow lavaan's literal default here.
    model.meta <- "f =~ 1*yi
                   f ~ mu*1
                   f ~~ tau2*f
                   yi ~~ data.vi*yi"
    RAM.meta <- lavaan2RAM(model.meta, obs.variables="yi", std.lv=FALSE)
    expect_identical(unname(RAM.meta$M["1","yi"]), "0")   ## still fixed, not free
    fit.meta <- sem(RAM=RAM.meta, data=Hox02, intervals="LB")
    expect_true(fit.meta$mx.fit@output$status$code %in% c(0, 1))
    expect_equal(unname(coef(fit.meta)["mu"]), 0.579, tolerance=1e-2)
})

test_that("an ordinary label merely resembling 'data.' (e.g. 'database') is not misclassified as a definition variable", {
    ## Regression test for a high-severity bug (external review):
    ## as.mxMatrix()/create.mxMatrix() used grep("data.", labels) to spot
    ## definition-variable labels (e.g. "data.vi", resolved per-row from a
    ## data column) and force them non-free -- but "." in a regex means
    ## "any character", so an ordinary label like "database" (data + b +
    ## ase) ALSO matched and was silently fixed instead of estimated. Only
    ## a literal "data." PREFIX should trigger this. Checked against a
    ## direct lavaan fit, in single-group, multiple-group, and two-level
    ## sem() fits (all three route through the same buggy match).
    set.seed(41)
    n <- 200
    x <- rnorm(n); y <- 2*x + rnorm(n)
    df <- data.frame(y=y, x=x)
    model <- "y ~ database*x
              y ~~ e*y"
    RAM <- lavaan2RAM(model)
    expect_true(RAM$A["y", "x"] %in% grep("^0(\\.[0-9]+)?\\*database$", RAM$A["y","x"], value=TRUE))
    fit <- sem(RAM=RAM, data=df)
    fit.lav <- lavaan::sem(model, data=df, fixed.x=FALSE)
    expect_true("database" %in% names(coef(fit)))
    expect_equal(unname(coef(fit)["database"]), unname(coef(fit.lav)["database"]), tolerance=1e-3)

    ## Multiple-group
    df2 <- data.frame(y=c(y, -y), x=c(x, -x), g=rep(c("g1","g2"), each=n))
    RAM2 <- lavaan2RAM("y ~ c(database1,database2)*x", ngroups=2)
    fit2 <- sem(RAM=RAM2, data=df2, group="g")
    expect_true(all(c("database1","database2") %in% names(coef(fit2))))
    expect_false(unname(coef(fit2)["database1"]) == unname(coef(fit2)["database2"]))

    ## Two-level
    nClusters <- 30; nPerCluster <- 5; N <- nClusters*nPerCluster
    clusterID <- rep(seq_len(nClusters), each=nPerCluster)
    xl <- rnorm(N); yl <- 0.5*xl + rnorm(N)
    df3 <- data.frame(clusterID=clusterID, x=xl, y=yl)
    model3 <- "level: 1
  y ~ database*x
  y ~~ residW*y
level: 2
  y ~ 1
  y ~~ residB*y"
    RAM3 <- lavaan2RAM(model3)
    fit3 <- sem(RAM=RAM3, data=df3, cluster="clusterID")
    expect_true("database" %in% names(coef(fit3)))

    ## A genuine definition variable ("data.vi") must still be correctly
    ## recognised and left non-free -- the fix must not have gone too far
    ## the other way.
    model4 <- "f =~ 1*yi
              f ~ mu*1
              f ~~ tau2*f
              yi ~~ data.vi*yi"
    RAM4 <- lavaan2RAM(model4, obs.variables="yi", std.lv=FALSE)
    fit4 <- sem(RAM=RAM4, data=Hox02)
    expect_false("data.vi" %in% names(coef(fit4)))
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

context("Checking sem() with multiple-group RAM objects")

test_that("sem() fits multiple-group regression matching lavaan", {
    set.seed(1)
    n1 <- 100; n2 <- 120
    x1 <- rnorm(n1); y1 <- 2 + 0.5*x1 + rnorm(n1, sd=0.8)
    x2 <- rnorm(n2); y2 <- -1 + 0.5*x2 + rnorm(n2, sd=1.2)
    df <- rbind(data.frame(y=y1, x=x1, g="g1"),
               data.frame(y=y2, x=x2, g="g2"))

    model <- "y ~ c(b,b)*x + c(a1,a2)*1"
    RAM <- lavaan2RAM(model, ngroups=2)
    expect_s3_class(RAM, "RAM_multigroup")

    fit <- sem(RAM=RAM, data=df, group="g")

    ## lavaan::sem() (not the lower-level lavaan::lavaan()) is the correct
    ## reference here: lavaan2RAM() documents observed variables' means as
    ## free by default, matching lavaan::sem()'s own int.ov.free=TRUE
    ## default -- raw lavaan::lavaan(..., fixed.x=FALSE) does NOT set
    ## int.ov.free=TRUE, so it instead fixes x's mean at 0 (with an
    ## explicit "automatically added intercepts are set to zero" warning),
    ## a DIFFERENT convention that happened to match an old lavaan2RAM()
    ## bug (fixed: lavaanify()'s own auto-generated, still-fixed-at-0 "~1"
    ## row for x was incorrectly allowed to overwrite lavaan2RAM()'s
    ## documented free-mean default). Confirmed by inspecting
    ## parTable()'s "free"/"user" columns directly for both lavaan() and
    ## lavaan::sem() on this exact model, not assumed.
    fit.lav <- lavaan::sem(model, data=df, group="g", fixed.x=FALSE)
    lav.m2ll <- -2*as.numeric(lavaan::fitMeasures(fit.lav, "logl"))

    expect_equal(fit$mx.fit$output$fit, lav.m2ll, tolerance=1e-3)

    cf <- coef(fit)
    ## The shared slope "b" is estimated once and equal across groups.
    expect_equal(unname(cf["b"]), unname(coef(fit.lav)["b"]), tolerance=1e-3)
    expect_equal(unname(cf["a1"]), unname(coef(fit.lav)["a1"]), tolerance=1e-3)
    expect_equal(unname(cf["a2"]), unname(coef(fit.lav)["a2"]), tolerance=1e-3)
})

test_that("sem() respects := and == constraints across groups", {
    ## b3 == b4 must equate group 1's and group 2's second slope; fn1 is a
    ## derived mxAlgebra of the (group-specific) first slopes.
    model1 <- "y ~ 1 + c(b1, b2)*x1 + c(b3, b4)*x2
               fn1 := b1 + b2
               b3 == b4"
    set.seed(2)
    n <- 150
    x1 <- rnorm(n); x2 <- rnorm(n); y <- 1 + 0.3*x1 + 0.3*x2 + rnorm(n)
    df <- data.frame(y=y, x1=x1, x2=x2, g=sample(c("g1","g2"), n, replace=TRUE))

    RAM1 <- lavaan2RAM(model1, ngroups=2)
    fit1 <- sem(RAM=RAM1, data=df, group="g")

    cf <- coef(fit1)
    expect_equal(unname(cf["b3"]), unname(cf["b4"]))

    ## The mxalgebra fn1 = b1 + b2 is reported and numerically correct.
    s1 <- summary(fit1)
    est <- s1$mxalgebras["fn1", "estimate"]
    expect_equal(unname(est), unname(cf["b1"] + cf["b2"]), tolerance=1e-6)

    ## fn1's SE/CI: mxSE()'s delta method (var(b1)+var(b2)+2*cov(b1,b2),
    ## via vcov(), since b1/b2 have independent per-group likelihoods here
    ## but this formula does not assume that) must agree with the reported
    ## Wald CI, since the CI is built as estimate +/- 1.96*SE internally --
    ## cross-checked against an INDEPENDENT computation (plain vcov()
    ## algebra, not mxSE()) rather than trusting mxSE()'s own output.
    V <- vcov(fit1)
    se.manual <- sqrt(V["b1","b1"] + V["b2","b2"] + 2*V["b1","b2"])
    se.implied <- (s1$mxalgebras["fn1","ubound"] - est)/1.96
    expect_equal(unname(se.implied), unname(se.manual), tolerance=1e-6)

    ## Likelihood-based (LB) CI for the same defined parameter: must exist,
    ## straddle the point estimate, and be centred close to it. Indexed as
    ## [row, col] throughout (not a single-row [row, ] slice), since the
    ## LB-intervals mxalgebras result is a data.frame (unlike the z-
    ## intervals path's matrix) -- a single-row [row, ] slice would stay a
    ## 1-row data.frame instead of dropping to a plain named vector.
    fit1.lb <- sem(RAM=RAM1, data=df, group="g", intervals.type="LB")
    mxalg.lb <- summary(fit1.lb)$mxalgebras
    expect_true(mxalg.lb["fn1","lbound"] < mxalg.lb["fn1","estimate"])
    expect_true(mxalg.lb["fn1","estimate"] < mxalg.lb["fn1","ubound"])
    expect_equal(unname(mxalg.lb["fn1","estimate"]), unname(est), tolerance=1e-3)
})

test_that("sem()'s := SE/CI is consistent for a nonlinear (indirect-effect) formula, in both single- and multiple-group fits, and respects robust=TRUE", {
    ## a*b (a mediation indirect effect) exercises mxSE()'s delta-method
    ## Jacobian on a genuinely nonlinear formula, not just a sum -- and,
    ## since a and b live in the SAME (sub)model with real (nonzero,
    ## correlated) sampling covariance, this is a much stronger check than
    ## the additive fn1 case above.
    set.seed(21)
    n <- 300
    m <- rnorm(n); x <- 0.5*m + rnorm(n); y <- 0.6*m + rnorm(n)

    ## -- Single-group --
    df1 <- data.frame(x=x, m=m, y=y)
    model.sg <- "m ~ a*x
                 y ~ b*m + cp*x
                 ind := a*b"
    fit.sg <- sem(RAM=lavaan2RAM(model.sg), data=df1)
    s.sg <- summary(fit.sg)
    a <- coef(fit.sg)["a"]; b <- coef(fit.sg)["b"]
    g <- c(b, a)   # d(a*b)/d(a,b) = (b, a)
    V <- vcov(fit.sg)[c("a","b"), c("a","b")]
    se.manual <- sqrt(as.numeric(t(g) %*% V %*% g))
    se.implied <- (s.sg$mxalgebras["ind","ubound"] - s.sg$mxalgebras["ind","estimate"])/1.96
    expect_equal(unname(se.implied), unname(se.manual), tolerance=1e-6)

    ## -- Multiple-group (2 groups, group-specific a/b/indirect effect) --
    df2 <- data.frame(x=x, m=m, y=y, g=sample(c("g1","g2"), n, replace=TRUE))
    model.mg <- "m ~ c(a1,a2)*x
                y ~ c(b1,b2)*m + c(cp1,cp2)*x
                ind1 := a1*b1
                ind2 := a2*b2"
    fit.mg <- sem(RAM=lavaan2RAM(model.mg, ngroups=2), data=df2, group="g")
    s.mg <- summary(fit.mg)
    a1 <- coef(fit.mg)["a1"]; b1 <- coef(fit.mg)["b1"]
    g <- c(b1, a1)
    V <- vcov(fit.mg)[c("a1","b1"), c("a1","b1")]
    se.manual <- sqrt(as.numeric(t(g) %*% V %*% g))
    se.implied <- (s.mg$mxalgebras["ind1","ubound"] - s.mg$mxalgebras["ind1","estimate"])/1.96
    expect_equal(unname(se.implied), unname(se.manual), tolerance=1e-6)

    ## robust=TRUE regression test: mxSE()'s own default silently ignores
    ## "robust" (it always falls back to the non-robust vcov(model) unless
    ## an explicit "cov=" is supplied), so a defined parameter's CI used to
    ## come out IDENTICAL whether robust=TRUE or FALSE, even though the
    ## ordinary coefficients' SEs correctly switched to robust ones --
    ## fixed by passing imxRobustSE(..., details=TRUE)$cov as mxSE()'s
    ## "cov" when robust=TRUE. Checked two ways: (a) the := CI must now
    ## actually differ between robust=TRUE/FALSE, and (b) it must match an
    ## independent (non-mxSE) delta-method computation using that same
    ## robust covariance matrix.
    s.mg.robust <- summary(fit.mg, robust=TRUE)
    expect_false(isTRUE(all.equal(s.mg.robust$mxalgebras["ind1", ],
                                  s.mg$mxalgebras["ind1", ])))
    robust.cov <- suppressMessages(imxRobustSE(fit.mg$mx.fit, details=TRUE))$cov
    Vr <- robust.cov[c("a1","b1"), c("a1","b1")]
    se.manual.robust <- sqrt(as.numeric(t(g) %*% Vr %*% g))
    se.implied.robust <- (s.mg.robust$mxalgebras["ind1","ubound"] -
                          s.mg.robust$mxalgebras["ind1","estimate"])/1.96
    expect_equal(unname(se.implied.robust), unname(se.manual.robust), tolerance=1e-6)
})

test_that("sem() errors clearly on unsupported multiple-group input", {
    model <- "y ~ c(b,b)*x"
    RAM <- lavaan2RAM(model, ngroups=2)

    set.seed(3)
    df <- data.frame(y=rnorm(20), x=rnorm(20), g=rep(c("g1","g2"), 10),
                     g3=rep(c("g1","g2","g3"), length.out=20))

    ## Missing 'group'
    expect_error(sem(RAM=RAM, data=df))
    ## 'group' not a column in data
    expect_error(sem(RAM=RAM, data=df, group="nosuchcolumn"))
    ## Number of groups implied by data does not match RAM
    expect_error(sem(RAM=RAM, data=df, group="g3"))
    ## replace.constraints=TRUE IS now supported for multiple-group RAM
    ## (see the dedicated tests below) -- with no "==" constraints in this
    ## model at all, it is simply a no-op and must NOT error.
    expect_error(sem(RAM=RAM, data=df, group="g", replace.constraints=TRUE), NA)
    ## Cov/numObs-based input not yet supported for multi-group RAM
    expect_error(sem(RAM=RAM, Cov=diag(2), numObs=100))

    ## Regression test: 'cluster' (a two-level-only argument) used to be
    ## silently ignored for a multi-group RAM instead of rejected -- a
    ## caller passing both 'group' and 'cluster' could be misled into
    ## thinking both dimensions were being modeled, when only 'group' was.
    expect_error(sem(RAM=RAM, data=df, group="g", cluster="g3"), "cluster")

    ## Regression test: an unused factor level in the group column used to
    ## silently pass split()'s group-count check (it creates an entry, with
    ## 0 rows, for every level) and produce a submodel with no data at all
    ## -- mxRun() did not even error on this, it returned a catastrophic-
    ## failure status with every parameter frozen at its starting value,
    ## easy to mistake for a converged (if odd) fit. Fixed by switching to
    ## unique(data[[group]]) (first-appearance order), which by
    ## construction never includes a level that doesn't actually occur --
    ## this now surfaces as an (equally clear) group-count mismatch error
    ## instead, one level up from where it used to be caught.
    df.empty.level <- data.frame(y=rnorm(20), x=rnorm(20),
                                 g=factor(rep("g1", 20), levels=c("g1","g2")))
    expect_error(sem(RAM=RAM, data=df.empty.level, group="g"))

    ## Regression test: rows with a missing group value are silently
    ## dropped by split() (a defensible default -- unlike the cluster case
    ## below, this does not create a wrong "phantom" group) but should
    ## still warn, not disappear without a trace.
    df.na.group <- data.frame(y=rnorm(20), x=rnorm(20),
                              g=rep(c("g1","g2"), 10))
    df.na.group$g[1:2] <- NA
    expect_warning(sem(RAM=RAM, data=df.na.group, group="g"), "missing")
})

test_that("sem() sanitizes group values that collide with OpenMx-reserved or internal names", {
    ## Regression test (external review): the sanitizer only replaced
    ## illegal CHARACTERS, so a group value that was already a legal name
    ## string could still collide with (a) mxModel()'s own globally
    ## reserved names ("data", "fitfunction", "expectation", "compute" --
    ## confirmed directly: mxModel("data") itself errors with "is illegal
    ## because it is a reserved name", regardless of context), (b) this
    ## call's own container model.name ("sem" by default -- a model
    ## cannot contain a child entity sharing its own name), or (c) the
    ## fixed internal matrix/algebra names .build.group.mxmodel() gives
    ## each group submodel ("Amatrix" etc. -- same self-containment
    ## problem once a submodel named "Amatrix" tries to hold a matrix
    ## also named "Amatrix"). All four surfaced as a cryptic OpenMx error
    ## rather than a submodel name collision the existing make.unique()
    ## machinery already knows how to resolve.
    set.seed(42)
    n <- 100
    x1 <- rnorm(n); y1 <- 2*x1 + rnorm(n)
    x2 <- rnorm(n); y2 <- -1*x2 + rnorm(n)
    RAM <- lavaan2RAM("y ~ c(b1,b2)*x", ngroups=2)

    for (reserved.value in c("sem", "data", "fitfunction", "Amatrix")) {
        df <- data.frame(y=c(y1,y2), x=c(x1,x2),
                         g=rep(c(reserved.value, "other"), each=n))
        fit <- expect_no_error(sem(RAM=RAM, data=df, group="g"))
        ## Not just "doesn't error" -- the reserved value's own submodel
        ## must still be locatable via group.map and give the right
        ## (first-appearing group's, here "b1" ~ slope 2) estimate.
        sub.name <- fit$group.map[[reserved.value]]
        expect_true(sub.name %in% names(fit$mx.fit@submodels))
        expect_equal(unname(coef(fit)["b1"]), 2, tolerance=0.5)
    }
})

test_that("sem() sanitizes group values that collide with a genuine free-parameter LABEL, with and without replace.constraints", {
    ## Regression test (external review): the "reserved" set covered
    ## OpenMx-/package-internal names (see the test above) but never a
    ## group's OWN free-parameter labels -- so a group value equal to a
    ## shared slope's label (e.g. "b" via lavaan's c(b,b) syntax) collided
    ## with that label at the OpenMx level, since a submodel is a NAMED
    ## ENTITY and OpenMx rejects a name used as both a named entity and a
    ## free parameter within it ("In model 'b' the following are used
    ## both as named entities and free parameters: 'b'"). Reproduced with
    ## AND without replace.constraints=TRUE (a broader multi-group bug,
    ## not specific to that feature).
    set.seed(1)
    n1 <- 100; n2 <- 100
    x1 <- rnorm(n1); y1 <- 2 + 0.5*x1 + rnorm(n1)
    x2 <- rnorm(n2); y2 <- -1 + 0.5*x2 + rnorm(n2)
    df <- rbind(data.frame(y=y1, x=x1, g="b"), data.frame(y=y2, x=x2, g="other"))
    model <- "y ~ c(b,b)*x + c(a1,a2)*1"
    RAM <- lavaan2RAM(model, ngroups=2)

    fit1 <- expect_no_error(sem(RAM=RAM, data=df, group="g"))
    expect_false(fit1$group.map[["b"]] == "b")   # sanitized away from the collision
    expect_equal(unname(coef(fit1)["b"]), 0.5, tolerance=0.2)

    fit2 <- expect_no_error(sem(RAM=RAM, data=df, group="g", replace.constraints=TRUE))
    expect_false(fit2$group.map[["b"]] == "b")
})

test_that("sem() sanitizes group values that collide with a label INTRODUCED by replace.constraints=TRUE", {
    ## Regression test (external review, follow-up to the test above): the
    ## first fix only reserved labels already present in the RAM matrices
    ## BEFORE substitution (para.labels). A label that only comes into
    ## existence via a replace.constraints=TRUE constraint's RHS (e.g. "a0"
    ## in "sigma1 == exp(a0)") was NOT reserved, because it used to be
    ## discovered only after group names were already finalized -- so a
    ## group literally named "a0" hit the same named-entity/free-parameter
    ## collision as an ordinary label ("In model 'a0' the following are
    ## used both as named entities and free parameters: 'a0'").
    set.seed(1)
    df <- rbind(
        data.frame(y=rnorm(120, 2, 1), g="a0"),
        data.frame(y=rnorm(120, -1, 1.3), g="g2"))
    model <- "
    y ~ c(mu1,mu2)*1
    y ~~ c(sigma1,sigma2)*y
    sigma1 == exp(a0)
    "
    fit1 <- expect_no_error(sem(RAM=lavaan2RAM(model, ngroups=2), data=df,
                                group="g", replace.constraints=TRUE))
    expect_false(fit1$group.map[["a0"]] == "a0")

    ## Same collision, but via a definition-variable-adjacent introduced
    ## label ("a1" in "exp(a0+a1*data.x)") -- "data.x" itself must NOT be
    ## reserved (it is a definition variable, never a free parameter).
    set.seed(2)
    df2 <- rbind(
        data.frame(y=rnorm(120, 2, 1), x=rnorm(120), g="a1"),
        data.frame(y=rnorm(120, -1, 1.3), x=rnorm(120), g="g2"))
    model2 <- "
    y ~ c(mu1,mu2)*1
    y ~~ c(sigma1,sigma2)*y
    sigma1 == exp(a0+a1*data.x)
    "
    fit2 <- expect_no_error(sem(RAM=lavaan2RAM(model2, ngroups=2), data=df2,
                                group="g", replace.constraints=TRUE))
    expect_false(fit2$group.map[["a1"]] == "a1")
})

test_that("sem() handles numeric group columns and 3+ groups", {
    ## Regression test: a purely-numeric group column (e.g. an integer site
    ## ID) produced submodel names like "1"/"2" via split(); OpenMx's
    ## mxModel() rejects a name that can be interpreted as a number
    ## ("The name '1' is illegal..."). Fixed by prefixing such names.
    set.seed(10)
    n <- 100
    x <- rnorm(n); y <- 0.5*x + rnorm(n)
    df.num <- data.frame(y=y, x=x, g=rep(c(1,2), n/2))
    RAM <- lavaan2RAM("y ~ c(b,b)*x", ngroups=2)
    fit.num <- sem(RAM=RAM, data=df.num, group="g")
    expect_identical(fit.num$mx.fit@output$status[[1]], 0L)
    expect_identical(sort(names(fit.num$mx.fit@submodels)), c("group_1","group_2"))

    ## 3 groups (not just the usual 2), all with a shared slope.
    set.seed(11)
    n3 <- 90
    x3 <- rnorm(n3); y3 <- 0.5*x3 + rnorm(n3)
    df3 <- data.frame(y=y3, x=x3, g=rep(c("g1","g2","g3"), n3/3))
    RAM3 <- lavaan2RAM("y ~ c(b,b,b)*x", ngroups=3)
    fit3 <- sem(RAM=RAM3, data=df3, group="g")
    expect_identical(fit3$mx.fit@output$status[[1]], 0L)
    expect_identical(names(fit3$mx.fit@submodels), c("g1","g2","g3"))
})

test_that("sem() orders groups by first appearance in the data, matching lavaan", {
    ## Regression test for a HIGH-severity bug found on external review: an
    ## earlier version of .sem.multigroup() split data via
    ## split(data, data[[group]]), which sorts group values alphabetically
    ## (or by a factor's levels(), if the column is a factor). lavaan does
    ## neither -- confirmed directly against lavaan() with deliberately
    ## reversed data order (group "b" appears before "a" in the rows) --
    ## it orders groups by first appearance in the data, full stop, even
    ## for a factor column with levels() in a different order. The old
    ## behaviour silently swapped which label's estimate landed on which
    ## group whenever the data's first-appearing group wasn't also the
    ## alphabetically-first one.
    set.seed(15)
    x1 <- rnorm(50); y1 <- 2 + 0.5*x1 + rnorm(50, sd=0.5)    # group "grpB": true int=2
    x2 <- rnorm(50); y2 <- -3 + 0.5*x2 + rnorm(50, sd=0.5)   # group "grpA": true int=-3
    df <- rbind(data.frame(y=y1, x=x1, g="grpB"),            # "grpB" appears FIRST
               data.frame(y=y2, x=x2, g="grpA"))

    RAM <- lavaan2RAM("y ~ c(slope,slope)*x + c(int_b,int_a)*1", ngroups=2)
    fit <- sem(RAM=RAM, data=df, group="g")

    ## Submodels, and thus which label got which group's data, follow
    ## first-appearance order: "grpB" (int_b) then "grpA" (int_a).
    expect_identical(names(fit$mx.fit@submodels), c("grpB","grpA"))
    cf <- coef(fit)
    expect_equal(unname(cf["int_b"]), 2, tolerance=0.5)
    expect_equal(unname(cf["int_a"]), -3, tolerance=0.5)

    ## Also true for a factor column whose levels() order disagrees with
    ## first appearance in the data.
    df.f <- df
    df.f$g <- factor(df.f$g, levels=c("grpA","grpB"))   # levels: A before B ...
    fit.f <- sem(RAM=RAM, data=df.f, group="g")
    expect_identical(names(fit.f$mx.fit@submodels), c("grpB","grpA"))  # ... but data: B first
})

test_that("sem() sanitizes illegal group names and keeps a mapping back to the original value", {
    ## Regression test for a medium-severity bug found on external review:
    ## the original numeric-only-name fix ("1" -> "group1") didn't cover
    ## other characters mxModel() rejects, e.g. "." and "-" (confirmed
    ## directly: "site.1"/"site-2" both fail with "The name '...' is
    ## illegal because it contains the '.'/'-' character").
    set.seed(16)
    n <- 80
    x <- rnorm(n); y <- 0.5*x + rnorm(n)
    df <- data.frame(y=y, x=x, g=rep(c("site.1","site-2"), n/2))
    RAM <- lavaan2RAM("y ~ c(b,b)*x", ngroups=2)
    fit <- sem(RAM=RAM, data=df, group="g")

    expect_identical(fit$mx.fit@output$status[[1]], 0L)
    sub.names <- names(fit$mx.fit@submodels)
    expect_true(all(grepl("^[A-Za-z0-9_]+$", sub.names)))  # legal mxModel names
    expect_identical(unname(fit$group.map[["site.1"]]), sub.names[1])
    expect_identical(unname(fit$group.map[["site-2"]]), sub.names[2])

    ## plot.mxsem(group=) accepts the ORIGINAL (unsanitized) value.
    skip_if_not_installed("semPlot")
    grDevices::pdf(nullfile())
    on.exit(grDevices::dev.off(), add=TRUE)
    expect_silent(plot(fit, group="site-2"))
})

test_that("plot.mxsem(group=) resolves sanitized-name collisions correctly", {
    ## Regression test for a medium-severity bug found on a follow-up
    ## external review: original groups "a.b" and "a_b" sanitize/
    ## de-duplicate to submodels "a_b" and "a_b_1". plot.mxsem() used to
    ## match against submodel names *before* group.map, so
    ## plot(fit, group="a_b") matched submodel "a_b" directly (since that
    ## string happens to already be a legal submodel name in its own
    ## right) -- silently selecting original group "a.b"'s data instead of
    ## the group literally named "a_b". Fixed by checking group.map first
    ## -- but naively doing that for *every* case also breaks the default
    ## (group=NULL) selection in this exact scenario, since sub.names[1]
    ## ("a_b") can itself coincidentally be a valid group.map key; the fix
    ## only consults group.map when the caller actually supplied "group".
    skip_if_not_installed("semPlot")
    set.seed(18)
    n <- 100
    x <- rnorm(n); y <- 0.5*x + rnorm(n)
    df <- data.frame(y=y, x=x, g=rep(c("a.b","a_b"), n/2))
    RAM <- lavaan2RAM("y ~ c(b,b)*x", ngroups=2)
    fit <- sem(RAM=RAM, data=df, group="g")

    sub.names <- names(fit$mx.fit@submodels)
    expect_identical(sub.names, c("a_b","a_b_1"))
    expect_identical(unname(fit$group.map[["a.b"]]), "a_b")
    expect_identical(unname(fit$group.map[["a_b"]]), "a_b_1")

    ## Correctness, not just "doesn't error": .resolve.plot.group() must
    ## pick the submodel matching the ORIGINAL group value, not whichever
    ## submodel happens to share that literal string.
    expect_identical(sub.names[.resolve.plot.group("a_b", sub.names, fit$group.map)],
                     "a_b_1")  # original group "a_b"
    expect_identical(sub.names[.resolve.plot.group("a.b", sub.names, fit$group.map)],
                     "a_b")    # original group "a.b"
    expect_identical(sub.names[.resolve.plot.group(NULL, sub.names, fit$group.map)],
                     "a_b")    # default: first submodel, i.e. original "a.b"

    grDevices::pdf(nullfile())
    on.exit(grDevices::dev.off(), add=TRUE)
    expect_silent(plot(fit, group="a_b"))
    expect_silent(plot(fit, group="a.b"))
    ## plot(fit) with neither "group" nor "level" set defaults to
    ## combine=TRUE (all panels), not the single-panel default resolution
    ## exercised directly above via .resolve.plot.group(); combine=FALSE
    ## is needed here to exercise that same default resolution through
    ## plot.mxsem() itself rather than only the unit test above.
    expect_silent(plot(fit))
    expect_silent(plot(fit, combine=FALSE))
})

test_that("sem() supports LB intervals, run=FALSE, and bounds for multi-group fits", {
    set.seed(12)
    n <- 100
    x <- rnorm(n); y <- 2*x + rnorm(n, sd=0.3)   # strong signal, true slope=2
    df <- data.frame(y=y, x=x, g=rep(c("g1","g2"), n/2))
    RAM <- lavaan2RAM("y ~ c(b,b)*x", ngroups=2)

    ## Likelihood-based CI
    fit.lb <- sem(RAM=RAM, data=df, group="g", intervals.type="LB")
    s.lb <- summary(fit.lb)
    expect_true(all(c("lbound","ubound") %in% colnames(s.lb$coefficients)))
    expect_true(s.lb$coefficients["b","lbound"] < s.lb$coefficients["b","ubound"])

    ## run=FALSE returns an unrun MxModel that can be run manually.
    fit.norun <- sem(RAM=RAM, data=df, group="g", run=FALSE)
    expect_s4_class(fit.norun$mx.fit, "MxModel")
    fit.manual <- mxRun(fit.norun$mx.fit, silent=TRUE)
    expect_identical(fit.manual@output$status[[1]], 0L)

    ## ubound actually constrains the optimizer, not just accepted silently.
    fit.bound <- sem(RAM=RAM, data=df, group="g", ubound=list(b=1))
    expect_equal(unname(coef(fit.bound)["b"]), 1, tolerance=1e-4)

    ## anova()/vcov() against a properly nested pair of fits.
    RAM.null <- lavaan2RAM("y ~ 0*x", ngroups=2)
    fit.null <- sem(RAM=RAM.null, data=df, group="g")
    cmp <- anova(fit.bound, fit.null)
    expect_equal(cmp$diffdf[2], 1)
    expect_equal(dim(vcov(fit.bound)), rep(length(coef(fit.bound)), 2))

    ## summary(robust=TRUE) is supported for multi-group (unlike two-level,
    ## see below): OpenMx's imxRobustSE() only needs mxFitFunctionMultigroup,
    ## which this path already uses.
    expect_silent(summary(sem(RAM=RAM, data=df, group="g"), robust=TRUE))
})

test_that("sem() supports replace.constraints=TRUE for multiple-group RAM within one group, fully-constant, and across groups", {
    ## .build.group.mxmodel() builds each group's A/S/M as ordinary named
    ## mxMatrix objects (dimnames=var.names passed directly to
    ## mxExpectationRAM()), the same construction style single-group sem()
    ## already uses for its own (already-shipped) replace.constraints
    ## branch -- so, unlike two-level, no remove-and-readd/dimname patch is
    ## needed here; the same unchanged/fully-constant/algebra branch is
    ## just applied per group.
    set.seed(1)
    n1 <- 200; n2 <- 220
    x1 <- rnorm(n1); y1 <- 2 + 0.5*x1 + rnorm(n1, sd=1)
    x2 <- rnorm(n2); y2 <- -1 + 0.5*x2 + rnorm(n2, sd=1.3)
    df <- rbind(data.frame(y=y1, x=x1, g="g1"), data.frame(y=y2, x=x2, g="g2"))

    model.plain <- "y ~ c(b,b)*x + c(a1,a2)*1
  y ~~ c(resid1,resid2)*y
  x ~~ c(varx1,varx2)*x
  x ~ c(mx1,mx2)*1"
    fit.plain <- sem(RAM=lavaan2RAM(model.plain, ngroups=2), data=df, group="g")
    m2ll.plain <- fit.plain$mx.fit@output$Minus2LogLikelihood
    cf.plain <- coef(fit.plain)

    ## -- Within one group (S matrix): resid1 == exp(a0). Pure
    ## reparameterisation, so -2LL and every OTHER estimate must match the
    ## unconstrained fit exactly, and exp(a0) must recover resid1.
    model.s <- paste(model.plain, "\nresid1 == exp(a0)")
    fit.s <- sem(RAM=lavaan2RAM(model.s, ngroups=2), data=df, group="g",
                replace.constraints=TRUE)
    expect_equal(fit.s$mx.fit@output$Minus2LogLikelihood, m2ll.plain, tolerance=1e-4)
    expect_equal(unname(exp(coef(fit.s)["a0"])), unname(cf.plain["resid1"]), tolerance=1e-4)

    ## -- A matrix, SHARED cross-group label: b == exp(logb) (b is the same
    ## free parameter in both groups via c(b,b)). Same reparameterisation
    ## logic, but the label being replaced lives in both groups at once.
    model.a <- paste(model.plain, "\nb == exp(logb)")
    fit.a <- sem(RAM=lavaan2RAM(model.a, ngroups=2), data=df, group="g",
                replace.constraints=TRUE)
    expect_equal(fit.a$mx.fit@output$Minus2LogLikelihood, m2ll.plain, tolerance=1e-4)
    expect_equal(unname(exp(coef(fit.a)["logb"])), unname(cf.plain["b"]), tolerance=1e-4)

    ## -- M matrix: a1 == exp(loga1) (group 1's own intercept).
    model.m <- paste(model.plain, "\na1 == exp(loga1)")
    fit.m <- sem(RAM=lavaan2RAM(model.m, ngroups=2), data=df, group="g",
                replace.constraints=TRUE)
    expect_equal(fit.m$mx.fit@output$Minus2LogLikelihood, m2ll.plain, tolerance=1e-4)
    expect_equal(unname(exp(coef(fit.m)["loga1"])), unname(cf.plain["a1"]), tolerance=1e-4)

    ## -- Cross-group chained constraint: resid1 == exp(a0) (group 1);
    ## resid2 == 2*resid1 (group 2), i.e. group 2's residual variance
    ## references group 1's OWN label -- no special handling needed for
    ## this since OpenMx labels are already global across a model tree.
    ## This genuinely RESTRICTS the model (resid2 is no longer free), so
    ## -2LL must be no better than (>=) the unconstrained fit, "resid1"
    ## must NOT reappear as an independent free parameter (the bug fixed
    ## earlier this session for two-level/single-group), and the
    ## constraint must hold exactly at the fitted values.
    model.cross <- paste(model.plain, "\nresid1 == exp(a0)\nresid2 == 2*resid1")
    fit.cross <- sem(RAM=lavaan2RAM(model.cross, ngroups=2), data=df, group="g",
                     replace.constraints=TRUE)
    expect_true(fit.cross$mx.fit@output$Minus2LogLikelihood >= m2ll.plain - 1e-6)
    cf.cross <- coef(fit.cross)
    expect_false("resid1" %in% names(cf.cross))
    g2.name <- fit.cross$group.map[["g2"]]
    resid2.fitted <- as.numeric(mxEval(Smatrix[1,1], fit.cross$mx.fit[[g2.name]], compute=TRUE))
    expect_equal(resid2.fitted, 2*exp(unname(cf.cross["a0"])), tolerance=1e-8)

    ## -- Fully-constant branch: fix resid1 at a literal number via "==",
    ## not a free algebra -- must be dropped as a free parameter entirely
    ## and genuinely restrict the model (worse or equal -2LL). Note: this
    ## particular model still has "varx1" free in group 1's S matrix, so
    ## this exercises the algebra branch (a single constant CELL inside an
    ## otherwise-free matrix), not the matrix-WIDE constant branch tested
    ## next -- both branches are exercised, deliberately, by these two
    ## different tests.
    model.const <- paste(model.plain, "\nresid1 == 5")
    fit.const <- sem(RAM=lavaan2RAM(model.const, ngroups=2), data=df, group="g",
                     replace.constraints=TRUE)
    expect_false("resid1" %in% names(coef(fit.const)))
    expect_true(fit.const$mx.fit@output$Minus2LogLikelihood >= m2ll.plain - 1e-6)

    ## -- Matrix-WIDE fully-constant branch with a non-trivial constant
    ## EXPRESSION (not a bare number): regression test for a high-severity
    ## bug (found on external review). When EVERY free cell in a group's
    ## A/S/M matrix gets eliminated, leaving the WHOLE matrix constant,
    ## the fully-constant branch used to hand the raw, unevaluated
    ## expression text straight to as.mxMatrix() -- which only recognises
    ## bare numbers or "value*label" strings, so a non-trivial expression
    ## like "exp(0.4)" (no "*" in it) was silently treated as an
    ## UNLABELLED FREE parameter starting at 0, not a fixed value.
    ## Confirmed by testing: without evaluating by hand first, this
    ## produced an anonymous free parameter fit to the data's own MLE
    ## instead of the constant 1.491825. A single-variable, single-group
    ## model (y is the ONLY variable, so S is 1x1) forces the WHOLE
    ## matrix, not just one cell, to go constant.
    set.seed(31)
    n1c <- 300; n2c <- 300
    dfc <- rbind(data.frame(y=rnorm(n1c, 2, 1), g="g1"),
                data.frame(y=rnorm(n2c, -1, 1.3), g="g2"))
    model.wholeconst <- "y ~ c(a1,a2)*1
  y ~~ c(sigma1,sigma2)*y
  sigma1 == exp(0.4)"
    fit.wc <- sem(RAM=lavaan2RAM(model.wholeconst, ngroups=2), data=dfc, group="g",
                 replace.constraints=TRUE)
    g1.name.wc <- fit.wc$group.map[["g1"]]
    sigma1.fitted <- as.numeric(mxEval(Smatrix[1,1], fit.wc$mx.fit[[g1.name.wc]], compute=TRUE))
    expect_equal(sigma1.fitted, exp(0.4), tolerance=1e-8)
    expect_false("sigma1" %in% names(coef(fit.wc)))

    ## -- Starting-value collision guard: regression test for a medium-
    ## severity bug (found on external review). The guard collected every
    ## group's own existing starting values via do.call(c, lapply(RAM,
    ## ...)) -- since RAM (a RAM_multigroup list) is itself named "1","2",
    ## ..., combining a NAMED outer list of named inner lists this way
    ## triggers R's automatic outer.inner name-concatenation, silently
    ## producing names like "1.sigma1" that never match anything in
    ## as.mxAlgebra()'s own startvalues lookup -- leaving the guard
    ## completely inert. "sigma2 == 2*sigma1" (sigma1 left as an
    ## ORDINARY, untouched free parameter, unlike the chained case above
    ## where sigma1 is ALSO replaced) is exactly the scenario that needs
    ## this guard: without it, this errors outright with "the free
    ## parameter 'sigma1' has been assigned multiple starting values".
    model.guard <- paste(model.plain, "\nresid2 == 2*resid1")
    fit.guard <- expect_no_error(
      expect_no_warning(sem(RAM=lavaan2RAM(model.guard, ngroups=2), data=df, group="g",
                            replace.constraints=TRUE)))
    cf.guard <- coef(fit.guard)
    g2.name.guard <- fit.guard$group.map[["g2"]]
    resid2.guard <- as.numeric(mxEval(Smatrix[1,1], fit.guard$mx.fit[[g2.name.guard]], compute=TRUE))
    expect_equal(resid2.guard, 2*unname(cf.guard["resid1"]), tolerance=1e-8)

    ## -- intervals.type="LB": documented as broken for single-group fits
    ## when combined with replace.constraints=TRUE (see sem()'s own
    ## @note); confirmed by testing that this does NOT carry over to
    ## multiple-group fits either (matching two-level's own result).
    fit.lb <- sem(RAM=lavaan2RAM(model.s, ngroups=2), data=df, group="g",
                 replace.constraints=TRUE, intervals.type="LB")
    ci <- summary(fit.lb)$coefficients
    expect_true(all(ci[, "lbound"] <= ci[, "Estimate"] + 1e-8, na.rm=TRUE))
    expect_true(all(ci[, "Estimate"] <= ci[, "ubound"] + 1e-8, na.rm=TRUE))

    ## -- replace.constraints=FALSE (the default) must be completely
    ## unaffected -- same fit as the very first call.
    fit.default <- sem(RAM=lavaan2RAM(model.plain, ngroups=2), data=df, group="g")
    expect_equal(fit.default$mx.fit@output$Minus2LogLikelihood, m2ll.plain, tolerance=1e-10)
})

test_that("sem() substitutes eliminated labels into REMAINING ':=' and unconsumed '==' or '<'/'>' entries, in all three RAM dispatch paths", {
    ## Regression test for a high-severity bug (found on external review,
    ## confirmed not unique to multi-group): once replace.constraints=TRUE
    ## eliminates a label (folding it into an algebra), a REMAINING lavaan
    ## ":=" defined parameter, or an unconsumed "=="/"<"/">" constraint,
    ## that still references that label by name used to be added
    ## UNCHANGED -- ":=" rows were never part of the "==" substitution set
    ## to begin with, in any of the three RAM dispatch paths. This failed
    ## at mxRun() with "Unknown reference '<label>' detected in the entity
    ## '<name>'" the moment the model tried to resolve it. Fixed via
    ## .substitute.remaining.mxalgebras(), applied identically in
    ## single-group sem(), .sem.twolevel(), and .sem.multigroup().

    ## -- Single-group --
    set.seed(1)
    n <- 400
    df.sg <- data.frame(y=rnorm(n, 2, sqrt(1.5)))
    model.sg <- "y ~ mu*1
  y ~~ sigma1*y
  sigma1 == exp(a0)
  foo := sigma1 + 3"
    fit.sg <- sem(RAM=lavaan2RAM(model.sg, obs.variables="y"), data=df.sg,
                 replace.constraints=TRUE)
    cf.sg <- coef(fit.sg)
    foo.sg <- summary(fit.sg)$mxalgebras["foo","estimate"]
    expect_equal(unname(foo.sg), unname(exp(cf.sg["a0"])+3), tolerance=1e-6)

    ## -- Two-level --
    set.seed(123)
    nClusters <- 60; nPerCluster <- 6
    clusterID <- rep(seq_len(nClusters), each=nPerCluster)
    u <- rep(rnorm(nClusters, 0, 0.7), each=nPerCluster)
    df.tl <- data.frame(clusterID=clusterID, y=1.5+u+rnorm(nClusters*nPerCluster))
    model.tl <- "level: 1
  y ~~ residW*y
level: 2
  y ~ 1
  y ~~ residB*y
residW == exp(a0)
foo := residW + residB"
    fit.tl <- sem(RAM=lavaan2RAM(model.tl), data=df.tl, cluster="clusterID",
                 replace.constraints=TRUE)
    cf.tl <- coef(fit.tl)
    foo.tl <- summary(fit.tl)$mxalgebras["foo","estimate"]
    expect_equal(unname(foo.tl), unname(exp(cf.tl["a0"])+cf.tl["residB"]), tolerance=1e-6)

    ## -- Multiple-group --
    set.seed(1)
    n1 <- 200; n2 <- 220
    x1 <- rnorm(n1); y1 <- 2 + 0.5*x1 + rnorm(n1, sd=1)
    x2 <- rnorm(n2); y2 <- -1 + 0.5*x2 + rnorm(n2, sd=1.3)
    df.mg <- rbind(data.frame(y=y1, x=x1, g="g1"), data.frame(y=y2, x=x2, g="g2"))
    model.mg <- "y ~ c(b,b)*x + c(a1,a2)*1
  y ~~ c(sigma1,sigma2)*y
  sigma1 == exp(a0)
  foo := sigma1 + sigma2"
    fit.mg <- sem(RAM=lavaan2RAM(model.mg, ngroups=2), data=df.mg, group="g",
                 replace.constraints=TRUE)
    cf.mg <- coef(fit.mg)
    foo.mg <- summary(fit.mg)$mxalgebras["foo","estimate"]
    expect_equal(unname(foo.mg), unname(exp(cf.mg["a0"])+cf.mg["sigma2"]), tolerance=1e-6)

    ## -- Unconsumed inequality ("<") constraint referencing an eliminated
    ## label: "<"/">" rows are NEVER in the "==" consumable set at all
    ## (only form[1]=="==" ever qualifies), so this is a genuinely
    ## different code path from the ":=" cases above.
    model.ineq <- "y ~ mu*1
  y ~~ sigma1*y
  sigma1 == exp(a0)
  mu < sigma1 + 100"
    expect_error(sem(RAM=lavaan2RAM(model.ineq, obs.variables="y"), data=df.sg,
                     replace.constraints=TRUE), NA)

    ## -- An entry that references NO eliminated label must be left
    ## completely alone (not needlessly rebuilt, and definitely not
    ## broken).
    model.unrelated <- "y ~ mu*1
  y ~~ sigma1*y
  x ~~ 1*x
  x ~ meanx*1
  sigma1 == exp(a0)
  bar := meanx + 3"
    df.unrelated <- data.frame(y=rnorm(300,2,1), x=rnorm(300))
    fit.unrel <- sem(RAM=lavaan2RAM(model.unrelated, obs.variables=c("y","x")),
                     data=df.unrelated, replace.constraints=TRUE)
    cf.unrel <- coef(fit.unrel)
    bar.unrel <- summary(fit.unrel)$mxalgebras["bar","estimate"]
    expect_equal(unname(bar.unrel), unname(cf.unrel["meanx"])+3, tolerance=1e-6)
})

test_that("sem() fits a model with NO remaining free parameters after replace.constraints=TRUE, in all three RAM dispatch paths", {
    ## Regression test (external review): mxCI() was called unconditionally
    ## on the (possibly empty) post-substitution label vector -- but
    ## mxCI(character(0)) (and mxCI(NULL)) both reject an empty reference
    ## vector outright ("'reference' argument must be a character vector")
    ## rather than treating it as "no CIs requested". Reachable whenever
    ## replace.constraints=TRUE fixes every "==" target at a literal
    ## constant, leaving nothing free at all -- a valid, if unusual, model
    ## to fit. Partly pre-existing (a fully-fixed model without any
    ## replace.constraints involved fails the same way), but
    ## replace.constraints=TRUE makes it easy to hit by accident. The fit
    ## itself, once reachable, must be numerically correct -- checked
    ## against a directly hand-computed -2LL, not just "doesn't error".
    ## OpenMx's own optimizer status $code is NA (not 0/1) when there is
    ## nothing to optimize, which is expected, not a symptom of a problem
    ## -- $status is still 0.

    ## -- Single-group --
    set.seed(1)
    df.sg <- data.frame(y=rnorm(200, 2, 1))
    model.sg <- "y ~ mu*1
  y ~~ sigma*y
  mu == 2
  sigma == exp(0.4)"
    fit.sg <- expect_no_error(
      sem(RAM=lavaan2RAM(model.sg, obs.variables="y"), data=df.sg,
         replace.constraints=TRUE))
    expect_length(coef(fit.sg), 0)
    expect_equal(fit.sg$mx.fit@output$status$status, 0)
    expect_equal(fit.sg$mx.fit@output$Minus2LogLikelihood,
                -2*sum(dnorm(df.sg$y, 2, sqrt(exp(0.4)), log=TRUE)), tolerance=1e-6)

    ## -- Multiple-group --
    df.mg <- rbind(data.frame(y=rnorm(100, 2, 1), g="g1"),
                   data.frame(y=rnorm(100, -1, 1), g="g2"))
    model.mg <- "y ~ c(a1,a2)*1
  y ~~ c(sigma1,sigma2)*y
  a1 == 2
  a2 == -1
  sigma1 == exp(0.4)
  sigma2 == exp(0.4)"
    fit.mg <- expect_no_error(
      sem(RAM=lavaan2RAM(model.mg, ngroups=2), data=df.mg, group="g",
         replace.constraints=TRUE))
    expect_length(coef(fit.mg), 0)
    manual.mg <- -2*sum(dnorm(df.mg$y[df.mg$g=="g1"], 2, sqrt(exp(0.4)), log=TRUE)) +
      -2*sum(dnorm(df.mg$y[df.mg$g=="g2"], -1, sqrt(exp(0.4)), log=TRUE))
    expect_equal(fit.mg$mx.fit@output$Minus2LogLikelihood, manual.mg, tolerance=1e-6)

    ## -- Two-level --
    nClusters <- 40; nPerCluster <- 5
    clusterID <- rep(seq_len(nClusters), each=nPerCluster)
    df.tl <- data.frame(clusterID=clusterID, y=rnorm(nClusters*nPerCluster, 2, 1))
    model.tl <- "level: 1
  y ~~ residW*y
level: 2
  y ~ 1
  y ~~ residB*y
residW == exp(0.4)
residB == exp(0.1)
ymean_between == 2"
    fit.tl <- expect_no_error(
      sem(RAM=lavaan2RAM(model.tl), data=df.tl, cluster="clusterID",
         replace.constraints=TRUE))
    expect_length(coef(fit.tl), 0)
    expect_equal(fit.tl$mx.fit@output$status$status, 0)
})

test_that("sem() supports a definition variable inside a replace.constraints=TRUE multiple-group constraint, resolved per group", {
    ## Each group already gets its own full data slice (mxData(data.list
    ## [[i]], type='raw')), unlike two-level's between level (which only
    ## ever carries the cluster ID) -- so, unlike two-level, there is no
    ## separate "genuine between-level covariate" gap to guard against
    ## here: a definition variable resolves against whichever group's own
    ## data it is used in, automatically.
    set.seed(21)
    n1 <- 400; n2 <- 400
    x1 <- rnorm(n1); a0.true <- -0.4; a1.true <- 0.7
    y1 <- 2 + rnorm(n1, 0, sd=sqrt(exp(a0.true + a1.true*x1)))
    x2 <- rnorm(n2)
    y2 <- -1 + rnorm(n2, 0, sd=sqrt(exp(a0.true + a1.true*x2)))
    df <- rbind(data.frame(y=y1, x=x1, g="g1"), data.frame(y=y2, x=x2, g="g2"))

    ## SAME a0/a1 shared across both groups' constraints (chained via
    ## identical labels), each resolved against that group's own data.
    model <- "y ~ c(mu1,mu2)*1
  y ~~ c(sigma1,sigma2)*y
  sigma1 == exp(a0+a1*data.x)
  sigma2 == exp(a0+a1*data.x)"
    fit <- sem(RAM=lavaan2RAM(model, ngroups=2), data=df, group="g",
              replace.constraints=TRUE)
    expect_identical(fit$mx.fit@output$status[[1]], 0L)
    cf <- coef(fit)
    expect_false(any(c("sigma1","sigma2") %in% names(cf)))
    expect_equal(unname(cf["a0"]), a0.true, tolerance=0.15)
    expect_equal(unname(cf["a1"]), a1.true, tolerance=0.15)
})

test_that("sem() rejects a cyclic replace.constraints=TRUE substitution across groups with a clear error", {
    set.seed(1)
    n <- 200
    df <- data.frame(y=rnorm(n), g=rep(c("g1","g2"), n/2))
    model <- "y ~ c(mu1,mu2)*1
  y ~~ c(sigma1,sigma2)*y
  sigma1 == 2*sigma2
  sigma2 == 2*sigma1"
    expect_error(sem(RAM=lavaan2RAM(model, ngroups=2), data=df, group="g",
                     replace.constraints=TRUE),
                "[Cc]yclic")
})

context("Checking sem() with two-level (within/between) RAM objects")

test_that("lavaan2RAM() parses two-level syntax into a within/between RAM", {
    model <- "level: 1
  fw =~ y1 + y2 + y3
level: 2
  fb =~ y1 + y2 + y3"
    RAM <- lavaan2RAM(model)
    expect_s3_class(RAM, "RAM_twolevel")
    expect_true(all(c("within","between") %in% names(RAM)))
    ## Within-level manifest means default to fixed at zero: the overall
    ## intercept lives at the between level (see dev/PLAN..., section 4).
    expect_true(all(RAM$within$M == "0"))
    ## Between-level "shadow" variables (y1/y2/y3, which are also level-1
    ## manifest variables) default to a free mean; the genuine level-2
    ## latent factor (fb) does not.
    expect_true(all(grepl("mean_between", RAM$between$M[1, c("y1","y2","y3")])))
    expect_identical(unname(RAM$between$M[1, "fb"]), "0")
    ## No manifest variables at the between level in this scope.
    expect_equal(nrow(RAM$between$F), 0)
})

test_that("sem() fits a two-level random-intercept regression matching lavaan", {
    set.seed(7)
    nClusters <- 60; nPerCluster <- 6; N <- nClusters*nPerCluster
    clusterID <- rep(seq_len(nClusters), each=nPerCluster)
    u <- rep(rnorm(nClusters, 0, 0.7), each=nPerCluster)
    x <- rnorm(N)
    y <- 1.5 + 0.6*x + u + rnorm(N)
    df <- data.frame(clusterID=clusterID, x=x, y=y)

    model <- "level: 1
  y ~ b1*x
  y ~~ residW*y
  x ~ meanX*1
  x ~~ varX*x
level: 2
  y ~ 1
  y ~~ residB*y"

    RAM <- lavaan2RAM(model)
    fit <- sem(RAM=RAM, data=df, cluster="clusterID")

    fit.lav <- suppressWarnings(
      lavaan::lavaan(model, data=df, cluster="clusterID", fixed.x=FALSE))
    lav.m2ll <- -2*as.numeric(lavaan::fitMeasures(fit.lav, "logl"))

    expect_equal(fit$mx.fit$output$fit, lav.m2ll, tolerance=1e-2)

    cf <- coef(fit)
    lav.cf <- coef(fit.lav)
    expect_equal(unname(cf["b1"]), unname(lav.cf["b1"]), tolerance=1e-2)
    expect_equal(unname(cf["residB"]), unname(lav.cf["residB"]), tolerance=1e-2)
})

test_that(".RAM2paths() keeps zero-prefixed definition-variable/@-tied cells", {
    ## Regression test for a bug found and fixed while implementing the
    ## two-level path: a fixed (free=FALSE) cell with a real label (a
    ## definition variable like "0*data.mod", or an "@"-tied fixed
    ## constant) but a numeric prefix of exactly 0 was being silently
    ## dropped (treated as "no path") by an earlier skip condition that
    ## only checked "!free && value==0", instead of "!free && value==0 &&
    ## no label at all". Only reachable via a hand-built RAM or explicit
    ## A.start=0/S.start=0 (lavaan2RAM()'s own defaults, A.start=0.1 and
    ## S.start=0.5, are never exactly 0), but a real correctness gap.
    A <- matrix(c("0", "0", "0*data.mod", "0"), nrow=2,
               dimnames=list(c("y","x"), c("y","x")))
    S <- matrix(c("0", "0", "0", "0"), nrow=2,
               dimnames=list(c("y","x"), c("y","x")))
    M <- matrix(c("0", "0"), nrow=1, dimnames=list(1, c("y","x")))

    built <- .RAM2paths(A, S, M)
    labels <- vapply(built$paths, function(p) {
      lb <- p@labels
      if (is.null(lb) || is.na(lb)) "" else lb
    }, character(1))
    expect_true("data.mod" %in% labels)
})

test_that("sem() fits two-level and multi-group models with definition variables", {
    ## Two-level: within-level residual variance is a per-row KNOWN
    ## definition variable (data.v), not estimated -- mirrors the package's
    ## own single-group idiom (e.g. metaFIML's "yi ~~ data.vi*yi").
    set.seed(11)
    nClusters <- 50; nPerCluster <- 5; N <- nClusters*nPerCluster
    clusterID <- rep(seq_len(nClusters), each=nPerCluster)
    u <- rep(rnorm(nClusters, 0, 0.6), each=nPerCluster)
    v <- runif(N, 0.3, 0.8)
    x <- rnorm(N)
    y <- 1 + 0.5*x + u + rnorm(N, sd=sqrt(v))
    df <- data.frame(clusterID=clusterID, x=x, y=y, v=v)

    model <- "level: 1
  y ~ b1*x
  y ~~ data.v*y
level: 2
  y ~ 1
  y ~~ residB*y"
    RAM <- lavaan2RAM(model)
    expect_identical(unname(RAM$within$S["y","y"]), "0.5*data.v")

    fit <- sem(RAM=RAM, data=df, cluster="clusterID")
    expect_identical(fit$mx.fit@output$status[[1]], 0L)
    ## The definition variable must not appear as an estimated parameter.
    expect_false(any(grepl("^data\\.v$|WITHy_within", names(coef(fit)))))

    ## Multi-group: same "data.v" convention, resolved against each
    ## group's own data.
    set.seed(12)
    n <- 200
    v2 <- runif(n, 0.3, 0.8)
    x2 <- rnorm(n)
    y2 <- 1 + 0.5*x2 + rnorm(n, sd=sqrt(v2))
    g <- rep(c("g1","g2"), n/2)
    df2 <- data.frame(y=y2, x=x2, v=v2, g=g)

    model2 <- "y ~ c(b,b)*x
              y ~~ data.v*y"
    ## "data.v" resolves against each group's own data despite sharing one
    ## textual label, so lavaanify's cross-group-equality warning here is a
    ## false positive for this specific (definition-variable) case.
    RAM2 <- suppressWarnings(lavaan2RAM(model2, ngroups=2))
    fit2 <- sem(RAM=RAM2, data=df2, group="g")
    expect_identical(fit2$mx.fit@output$status[[1]], 0L)
    expect_false(any(grepl("^data\\.v$", names(coef(fit2)))))
})

test_that("sem() respects := constraints for two-level RAM, with correct SE/CI", {
    set.seed(8)
    nClusters <- 40; nPerCluster <- 5; N <- nClusters*nPerCluster
    clusterID <- rep(seq_len(nClusters), each=nPerCluster)
    u <- rep(rnorm(nClusters, 0, 0.5), each=nPerCluster)
    x <- rnorm(N)
    y <- 1 + 0.4*x + u + rnorm(N)
    df <- data.frame(clusterID=clusterID, x=x, y=y)

    ## fn1 (additive) and fn2 (product, nonlinear -- exercises mxSE()'s
    ## delta-method Jacobian more rigorously) both combine a WITHIN-level
    ## parameter (b1) with a parameter estimated at the SAME within-level
    ## submodel (residW), i.e. genuinely correlated free parameters, not
    ## independent ones.
    model <- "level: 1
  y ~ b1*x
  y ~~ residW*y
level: 2
  y ~ 1
  y ~~ residB*y
fn1 := b1 + residW
fn2 := b1 * residW"

    RAM <- lavaan2RAM(model)
    fit <- sem(RAM=RAM, data=df, cluster="clusterID")

    cf <- coef(fit)
    s <- summary(fit)
    est1 <- s$mxalgebras["fn1", "estimate"]
    expect_equal(unname(est1), unname(cf["b1"] + cf["residW"]), tolerance=1e-6)
    est2 <- s$mxalgebras["fn2", "estimate"]
    expect_equal(unname(est2), unname(cf["b1"] * cf["residW"]), tolerance=1e-6)

    ## SE/CI for both, cross-checked against an INDEPENDENT (non-mxSE)
    ## delta-method computation built from vcov() directly -- this is the
    ## part most worth checking for a two-level fit specifically, since
    ## the relational (primaryKey/joinKey) "between" submodel is a less
    ## standard OpenMx configuration than a flat or multigroup model.
    V <- vcov(fit)[c("b1","residW"), c("b1","residW")]
    se1.manual <- sqrt(V["b1","b1"] + V["residW","residW"] + 2*V["b1","residW"])
    se1.implied <- (s$mxalgebras["fn1","ubound"] - est1)/1.96
    expect_equal(unname(se1.implied), unname(se1.manual), tolerance=1e-6)

    b1 <- cf["b1"]; residW <- cf["residW"]
    grad2 <- c(residW, b1)   # d(b1*residW)/d(b1,residW) = (residW, b1)
    se2.manual <- sqrt(as.numeric(t(grad2) %*% V %*% grad2))
    se2.implied <- (s$mxalgebras["fn2","ubound"] - est2)/1.96
    expect_equal(unname(se2.implied), unname(se2.manual), tolerance=1e-6)

    ## Likelihood-based (LB) CI for the same defined parameters.
    fit.lb <- sem(RAM=RAM, data=df, cluster="clusterID", intervals.type="LB")
    ci.lb <- summary(fit.lb)$mxalgebras
    expect_true(all(ci.lb[, "lbound"] < ci.lb[, "estimate"]))
    expect_true(all(ci.lb[, "estimate"] < ci.lb[, "ubound"]))
    expect_equal(unname(ci.lb["fn1", "estimate"]), unname(est1), tolerance=1e-3)
    expect_equal(unname(ci.lb["fn2", "estimate"]), unname(est2), tolerance=1e-3)

    ## robust=TRUE is rejected outright for two-level fits (see
    ## summary.mxsem()) -- := parameters are no exception, so there is no
    ## separate "robust SE ignored" failure mode to worry about here, only
    ## confirming the rejection still fires with a := present in the model.
    expect_error(summary(fit, robust=TRUE), "two-level")
})

test_that("sem() supports replace.constraints=TRUE for two-level RAM at the within level, the between level, and across both", {
    ## .sem.twolevel() does not build named A/S/M mxMatrix objects the way
    ## the single-group path does -- it flattens everything to mxPath()
    ## calls, because OpenMx's relational (primaryKey/joinKey) mechanism
    ## was found (by testing) to silently break when a plain mxMatrix is
    ## supplied instead. replace.constraints=TRUE works around this by
    ## removing the auto-built matrix and re-adding an as.mxAlgebra()
    ## decomposition with its dimnames set explicitly (also required, and
    ## also confirmed necessary by testing -- unlike the single-group path,
    ## which instead passes dimnames=var.names directly to
    ## mxExpectationRAM() and needs no dimnames on the matrices themselves).
    set.seed(123)
    nClusters <- 60; nPerCluster <- 6; N <- nClusters*nPerCluster
    clusterID <- rep(seq_len(nClusters), each=nPerCluster)
    u <- rep(rnorm(nClusters, 0, 0.7), each=nPerCluster)
    x <- rnorm(N)
    y <- 1.5 + 0.6*x + u + rnorm(N)
    df <- data.frame(clusterID=clusterID, x=x, y=y)

    model.plain <- "level: 1
  y ~ b1*x
  y ~~ residW*y
level: 2
  y ~ 1
  y ~~ residB*y"
    fit.plain <- sem(RAM=lavaan2RAM(model.plain), data=df, cluster="clusterID")
    m2ll.plain <- fit.plain$mx.fit@output$Minus2LogLikelihood
    cf.plain <- coef(fit.plain)

    ## -- Within-level S: residW == exp(logResidW). A pure reparameterisation
    ## (not a restriction), so -2LL and every OTHER estimate must match the
    ## unconstrained fit exactly, and exp(logResidW) must recover residW.
    model.withinS <- paste(model.plain, "\nresidW == exp(logResidW)")
    fit.ws <- sem(RAM=lavaan2RAM(model.withinS), data=df, cluster="clusterID",
                 replace.constraints=TRUE)
    expect_equal(fit.ws$mx.fit@output$Minus2LogLikelihood, m2ll.plain, tolerance=1e-4)
    expect_equal(unname(exp(coef(fit.ws)["logResidW"])), unname(cf.plain["residW"]),
                tolerance=1e-4)

    ## -- Within-level A: b1 == exp(logB1). Same reparameterisation logic.
    model.withinA <- paste(model.plain, "\nb1 == exp(logB1)")
    fit.wa <- sem(RAM=lavaan2RAM(model.withinA), data=df, cluster="clusterID",
                 replace.constraints=TRUE)
    expect_equal(fit.wa$mx.fit@output$Minus2LogLikelihood, m2ll.plain, tolerance=1e-4)
    expect_equal(unname(exp(coef(fit.wa)["logB1"])), unname(cf.plain["b1"]), tolerance=1e-4)

    ## -- Between-level M: ymean_between == exp(logMean). Same logic, this
    ## time on the mean structure of the "between" submodel.
    model.betweenM <- paste(model.plain, "\nymean_between == exp(logMean)")
    fit.bm <- sem(RAM=lavaan2RAM(model.betweenM), data=df, cluster="clusterID",
                 replace.constraints=TRUE)
    expect_equal(fit.bm$mx.fit@output$Minus2LogLikelihood, m2ll.plain, tolerance=1e-4)
    expect_equal(unname(exp(coef(fit.bm)["logMean"])), unname(cf.plain["ymean_between"]),
                tolerance=1e-4)

    ## -- Cross-level constraint: residB == 2*residW, i.e. the between-level
    ## algebra references a label that lives only as an ordinary free
    ## parameter in the WITHIN model -- no special handling is needed for
    ## this since OpenMx labels are already global across a model tree.
    ## This genuinely RESTRICTS the model (residB is no longer free), so
    ## -2LL must be no better than (i.e. >=) the unconstrained fit, and the
    ## constraint must hold exactly at the fitted values.
    model.cross <- paste(model.plain, "\nresidB == 2*residW")
    fit.cross <- sem(RAM=lavaan2RAM(model.cross), data=df, cluster="clusterID",
                     replace.constraints=TRUE)
    expect_true(fit.cross$mx.fit@output$Minus2LogLikelihood >= m2ll.plain - 1e-6)
    cf.cross <- coef(fit.cross)
    expect_false("residB" %in% names(cf.cross))
    ## residB itself is reported via mxCI as a defined quantity is not
    ## automatic here (it was fully absorbed into the S matrix, not a ":="),
    ## so recompute it directly from the fitted "S" algebra instead.
    mxEval.residB <- mxEval(S[1,1], fit.cross$mx.fit$between, compute=TRUE)
    expect_equal(as.numeric(mxEval.residB), 2*unname(cf.cross["residW"]), tolerance=1e-6)

    ## -- Fully-constant branch: fix residW at a literal number via "==",
    ## not a free algebra -- must be dropped as a free parameter entirely
    ## (not left dangling as an unresolved algebra) and genuinely restrict
    ## the model (worse or equal -2LL).
    model.const <- paste(model.plain, "\nresidW == 5")
    fit.const <- sem(RAM=lavaan2RAM(model.const), data=df, cluster="clusterID",
                     replace.constraints=TRUE)
    expect_false("residW" %in% names(coef(fit.const)))
    expect_true(fit.const$mx.fit@output$Minus2LogLikelihood >= m2ll.plain - 1e-6)

    ## -- intervals.type="LB": this exact combination (replace.constraints=
    ## TRUE + "LB") is documented as BROKEN for single-group fits (see
    ## sem()'s own @note) -- confirmed by testing that this limitation does
    ## NOT carry over to two-level fits.
    fit.lb <- sem(RAM=lavaan2RAM(model.withinS), data=df, cluster="clusterID",
                 replace.constraints=TRUE, intervals.type="LB")
    ci <- summary(fit.lb)$coefficients
    expect_true(all(ci[, "lbound"] <= ci[, "Estimate"] + 1e-8, na.rm=TRUE))
    expect_true(all(ci[, "Estimate"] <= ci[, "ubound"] + 1e-8, na.rm=TRUE))

    ## -- replace.constraints=FALSE (the default) must be completely
    ## unaffected by any of the above -- same fit as the very first call.
    fit.default <- sem(RAM=lavaan2RAM(model.plain), data=df, cluster="clusterID")
    expect_equal(fit.default$mx.fit@output$Minus2LogLikelihood, m2ll.plain, tolerance=1e-10)

    ## -- Fully-constant M replacement: regression test for a bug (found
    ## on external review) where the constant branch of
    ## .twolevel.replace.one() built its numeric mxMatrix via
    ## as.mxMatrix(), which silently SKIPS setting dimnames whenever
    ## rownames(x)[1]=="1" -- always true for an M (mean) matrix. Left
    ## unfixed, this reached mxRun() as "the M matrix ... does not
    ## contain dimnames" for any constant mean replacement, e.g. fixing
    ## the between-level intercept at a literal value.
    model.constM <- paste(model.plain, "\nymean_between == 1")
    fit.constM <- sem(RAM=lavaan2RAM(model.constM), data=df, cluster="clusterID",
                      replace.constraints=TRUE)
    expect_identical(fit.constM$mx.fit@output$status[[1]], 0L)
    expect_false("ymean_between" %in% names(coef(fit.constM)))
    expect_equal(as.numeric(mxEval(M[1,1], fit.constM$mx.fit$between, compute=TRUE)), 1)
})

test_that("sem()'s single-group fully-constant replace.constraints=TRUE branch evaluates a constant EXPRESSION, not just a bare number", {
    ## Regression test for a high-severity bug (found on external review):
    ## the fully-constant branch (all free symbols eliminated from A/S/M)
    ## handed the raw, un-evaluated expression TEXT straight to
    ## as.mxMatrix(), which only recognises bare numbers or "value*label"
    ## strings -- silently treating a non-trivial constant expression like
    ## "exp(0.4)" (no "*" in it) as an UNLABELLED FREE parameter with no
    ## starting value, not a fixed value. This bug was already fixed in
    ## .twolevel.replace.one() and .build.group.mxmodel() (multi-group)
    ## earlier this session, but the identical single-group branch (the
    ## one those two were themselves copied from) was left unfixed at the
    ## time as believed out of scope -- found still present on a
    ## follow-up external review. Fixed identically in all three A/S/M
    ## branches here.
    set.seed(1)
    n <- 400
    df.s <- data.frame(y=rnorm(n, 2, sqrt(1.5)))
    model.s <- "y ~ mu*1
              y ~~ sigma*y
              sigma == exp(0.4)"
    fit.s <- sem(RAM=lavaan2RAM(model.s, obs.variables="y"), data=df.s,
                replace.constraints=TRUE)
    sigma.fitted <- as.numeric(mxEval(Smatrix[1,1], fit.s$mx.fit, compute=TRUE))
    expect_equal(sigma.fitted, exp(0.4), tolerance=1e-8)
    expect_false("sigma" %in% names(coef(fit.s)))

    set.seed(2)
    n2 <- 400
    x <- rnorm(n2); y2 <- 0.4*x + rnorm(n2)
    df.a <- data.frame(y=y2, x=x)
    model.a <- "y ~ b*x
               x ~~ 1*x
               b == exp(0.4)"
    fit.a <- sem(RAM=lavaan2RAM(model.a, obs.variables=c("y","x")), data=df.a,
                replace.constraints=TRUE)
    b.fitted <- as.numeric(mxEval(Amatrix[1,2], fit.a$mx.fit, compute=TRUE))
    expect_equal(b.fitted, exp(0.4), tolerance=1e-8)
    expect_false("b" %in% names(coef(fit.a)))

    set.seed(3)
    n3 <- 400
    df.m <- data.frame(y=rnorm(n3, 2, 1))
    model.m <- "y ~ mu*1
               y ~~ sigma*y
               mu == exp(0.4)"
    fit.m <- sem(RAM=lavaan2RAM(model.m, obs.variables="y"), data=df.m,
                replace.constraints=TRUE)
    mu.fitted <- as.numeric(mxEval(Mmatrix[1,1], fit.m$mx.fit, compute=TRUE))
    expect_equal(mu.fitted, exp(0.4), tolerance=1e-8)
    expect_false("mu" %in% names(coef(fit.m)))
})

test_that("sem() resolves CHAINED replace.constraints=TRUE substitutions correctly, in both single-group and two-level fits", {
    ## Regression test for a high-severity bug (found on external review):
    ## the substitution loop applied each "==" constraint's OWN,
    ## un-expanded RHS text in a single pass over the RAM matrices. When
    ## one constraint's RHS referenced another constraint's OWN LHS label
    ## (a chain), the referenced label -- already eliminated from its own
    ## original location -- was reinserted verbatim and, having no
    ## remaining declaration anywhere, was silently treated by
    ## as.mxAlgebra() as a brand-new, DISCONNECTED free parameter. E.g.
    ## "residW == exp(a0); residB == 2*residW" ended up fitting a
    ## genuinely different (two-parameter) model than the one written,
    ## with no error or warning. Fixed via .resolve.chained.constraints(),
    ## which expands every constraint's RHS to a fixed point (order-
    ## independent, with cycle/self-reference detection) before any of
    ## them are applied to the matrices.

    ## -- Two-level: within residW == exp(a0); between residB == 2*residW.
    set.seed(123)
    nClusters <- 60; nPerCluster <- 6; N <- nClusters*nPerCluster
    clusterID <- rep(seq_len(nClusters), each=nPerCluster)
    u <- rep(rnorm(nClusters, 0, 0.7), each=nPerCluster)
    x <- rnorm(N)
    y <- 1.5 + 0.6*x + u + rnorm(N)
    df <- data.frame(clusterID=clusterID, x=x, y=y)

    model <- "level: 1
  y ~ b1*x
  y ~~ residW*y
level: 2
  y ~ 1
  y ~~ residB*y
residW == exp(a0)
residB == 2*residW"
    fit <- sem(RAM=lavaan2RAM(model), data=df, cluster="clusterID",
              replace.constraints=TRUE)
    expect_true(fit$mx.fit@output$status[[1]] %in% c(0,1))
    cf <- coef(fit)
    ## The eliminated label must NOT resurface as its own free parameter.
    expect_false("residW" %in% names(cf))
    ## The between-level S cell (residB) must equal 2*exp(a0) EXACTLY --
    ## not merely "close", since this is the same deterministic quantity
    ## computed two different ways (fitted algebra vs. recomputed by hand).
    residB.fitted <- as.numeric(mxEval(S[1,1], fit$mx.fit$between, compute=TRUE))
    expect_equal(residB.fitted, 2*exp(unname(cf["a0"])), tolerance=1e-8)

    ## Order must not matter: writing the chain in the OPPOSITE dependency
    ## direction gives an identical fit.
    model.rev <- "level: 1
  y ~ b1*x
  y ~~ residW*y
level: 2
  y ~ 1
  y ~~ residB*y
residB == 2*residW
residW == exp(a0)"
    fit.rev <- sem(RAM=lavaan2RAM(model.rev), data=df, cluster="clusterID",
                  replace.constraints=TRUE)
    expect_equal(fit.rev$mx.fit@output$Minus2LogLikelihood,
                fit$mx.fit@output$Minus2LogLikelihood, tolerance=1e-6)

    ## -- Single-group: two dimensionally-consistent, independently
    ## identified variance parameters chained together (sigma2_y on y's
    ## residual, sigma2_x on x's variance) -- well-behaved and directly
    ## checkable against sample variances, unlike an arbitrary constraint.
    set.seed(6)
    n4 <- 800
    y4 <- rnorm(n4, 3, sqrt(1.5))
    x4 <- rnorm(n4, 0, sqrt(2*1.5))
    df4 <- data.frame(y=y4, x=x4)
    model4 <- "y ~ mu*1
  y ~~ sigma2_y*y
  x ~~ sigma2_x*x
  x ~ meanx*1
  sigma2_y == exp(a0)
  sigma2_x == 2*sigma2_y"
    fit4 <- sem(RAM=lavaan2RAM(model4, obs.variables=c("y","x")), data=df4,
               replace.constraints=TRUE)
    expect_true(fit4$mx.fit@output$status[[1]] %in% c(0,1))
    cf4 <- coef(fit4)
    expect_false("sigma2_y" %in% names(cf4))
    sigma2x.fitted <- as.numeric(mxEval(Smatrix[2,2], fit4$mx.fit, compute=TRUE))
    expect_equal(sigma2x.fitted, 2*exp(unname(cf4["a0"])), tolerance=1e-8)
    expect_equal(exp(unname(cf4["a0"])), 1.5, tolerance=0.3)

    ## -- Long RHS (multiple moderators): regression test for a bug (found
    ## on external review) where .resolve.chained.constraints()'s final
    ## vapply(exprs, deparse, character(1)) failed outright for any
    ## resolved expression long enough that deparse() wraps it onto more
    ## than one line ("values must be length 1, but FUN(X[[1]]) result is
    ## length N") -- a perfectly valid constraint, rejected purely because
    ## of its length, before the model was even built. Fixed via
    ## deparse1(), which always collapses onto one line. Six moderators is
    ## comfortably past deparse()'s default line-wrap width.
    set.seed(9)
    n5 <- 500
    xs <- as.data.frame(matrix(rnorm(n5*6), ncol=6))
    names(xs) <- paste0("x", 1:6)
    df5 <- cbind(y=rnorm(n5, 2), xs)
    model5 <- paste0(
      "y ~ mu*1\n y ~~ sigma2*y\n",
      "sigma2 == exp(a0 + ", paste0("a", 1:6, "*data.x", 1:6, collapse=" + "), ")")
    fit5 <- sem(RAM=lavaan2RAM(model5, obs.variables="y"), data=df5,
               replace.constraints=TRUE)
    expect_true(fit5$mx.fit@output$status[[1]] %in% c(0,1))
    cf5 <- coef(fit5)
    expect_false("sigma2" %in% names(cf5))
    expect_true(all(paste0("a", 0:6) %in% names(cf5)))

    ## -- Cyclic/self-referential constraints must be rejected with a
    ## clear error, not silently mis-fitted or hang.
    model.cycle <- "y ~ mu*1
  y ~~ sigma2_y*y
  x ~~ sigma2_x*x
  x ~ meanx*1
  sigma2_y == 2*sigma2_x
  sigma2_x == 2*sigma2_y"
    expect_error(sem(RAM=lavaan2RAM(model.cycle, obs.variables=c("y","x")), data=df4,
                     replace.constraints=TRUE),
                "[Cc]yclic")

    model.selfref <- "y ~ mu*1
  y ~~ sigma2_y*y
  sigma2_y == sigma2_y + 1"
    expect_error(sem(RAM=lavaan2RAM(model.selfref, obs.variables="y"), data=df4[,"y",drop=FALSE],
                     replace.constraints=TRUE),
                "[Cc]yclic")
})

test_that("sem() supports a definition variable inside a replace.constraints=TRUE two-level constraint (within-level location-scale heterogeneity)", {
    ## The realistic use case: making the WITHIN-study heterogeneity a
    ## function of an effect-size-level moderator, i.e. a two-level
    ## location-scale meta-analysis (paralleling the single-group
    ## location-scale model, e.g. "sigma2 == exp(a0+a1*data.xi)"). "tauW"
    ## never appears as a free parameter itself: lavaanify() would
    ## otherwise reject a0/a1 as undeclared, and an ordinary (non-replaced)
    ## "==" constraint cannot introduce new free parameters this way --
    ## only replace.constraints=TRUE can.
    set.seed(7)
    nClusters <- 150; nPerCluster <- 8; N <- nClusters*nPerCluster
    clusterID <- rep(seq_len(nClusters), each=nPerCluster)
    u <- rep(rnorm(nClusters, 0, 0.5), each=nPerCluster)
    xi <- rnorm(N)
    a0.true <- -0.3; a1.true <- 0.6
    e <- rnorm(N, 0, sd=sqrt(exp(a0.true + a1.true*xi)))
    y <- 1.5 + u + e
    df <- data.frame(clusterID=clusterID, xi=xi, y=y)

    model <- "level: 1
  y ~~ residW*y
level: 2
  y ~ 1
  y ~~ residB*y
residW == exp(a0 + a1*data.xi)"
    RAM <- lavaan2RAM(model, obs.variables="y", std.lv=FALSE)
    ## "data.xi" must be flagged as a fixed definition variable (free=FALSE)
    ## in the auto-built parameter matrix, not mistaken for a free
    ## parameter to estimate -- same convention as an ordinary "0*data.xi"
    ## RAM cell, but here reached via as.mxAlgebra()'s own machinery.
    fit <- sem(RAM=RAM, data=df, cluster="clusterID", replace.constraints=TRUE)
    expect_identical(fit$mx.fit@output$status[[1]], 0L)

    cf <- coef(fit)
    expect_false("data.xi" %in% names(cf))
    expect_equal(unname(cf["a0"]), a0.true, tolerance=0.15)
    expect_equal(unname(cf["a1"]), a1.true, tolerance=0.15)

    ## Independent cross-check: plugging the TRUE per-row heterogeneity in
    ## directly (an ordinary, non-algebra definition variable) must give a
    ## -2LL close to, and no better than, the fitted (MLE) model's -2LL --
    ## confirming "data.xi" is genuinely being resolved per-row and not,
    ## say, silently defaulting to a constant.
    df$residWtrue <- exp(a0.true + a1.true*df$xi)
    model.true <- "level: 1
  y ~~ data.residWtrue*y
level: 2
  y ~ 1
  y ~~ residB*y"
    fit.true <- sem(RAM=lavaan2RAM(model.true, obs.variables="y", std.lv=FALSE),
                    data=df, cluster="clusterID")
    m2ll.fit <- fit$mx.fit@output$Minus2LogLikelihood
    m2ll.true <- fit.true$mx.fit@output$Minus2LogLikelihood
    expect_true(m2ll.fit <= m2ll.true + 1e-6)
    expect_equal(m2ll.fit, m2ll.true, tolerance=2)
})

test_that("sem() rejects a genuine between-level (cluster-level) definition variable with a clear, actionable error", {
    ## "between.data" (built inside .sem.twolevel()) only ever carries the
    ## cluster ID column -- a genuine cluster-level covariate is not yet
    ## plumbed through. Left unchecked, this would reach mxRun() as a hard
    ## R error that the mxRun-then-mxTryHard() fallback mistakes for
    ## "maybe just needs a better start", burning through 50 retries before
    ## re-raising OpenMx's own cryptic message. Must instead be caught
    ## immediately with a message that actually explains the limitation.
    set.seed(7)
    nClusters <- 100; nPerCluster <- 5; N <- nClusters*nPerCluster
    clusterID <- rep(seq_len(nClusters), each=nPerCluster)
    z <- rep(rnorm(nClusters), each=nPerCluster)
    u <- rep(rnorm(nClusters, 0, 0.5), each=nPerCluster)
    y <- 1.5 + u + rnorm(N)
    df <- data.frame(clusterID=clusterID, z=z, y=y)

    ## Case 1: written directly in the model (no replace.constraints).
    model1 <- "level: 1
  y ~~ residW*y
level: 2
  y ~ 1
  y ~~ data.z*y"
    RAM1 <- lavaan2RAM(model1, obs.variables="y", std.lv=FALSE)
    expect_error(sem(RAM=RAM1, data=df, cluster="clusterID"),
                "level-2-only.*definition variable")

    ## Case 2: introduced only via a replace.constraints substitution.
    model2 <- "level: 1
  y ~~ residW*y
level: 2
  y ~ 1
  y ~~ residB*y
residB == exp(c0 + c1*data.z)"
    RAM2 <- lavaan2RAM(model2, obs.variables="y", std.lv=FALSE)
    expect_error(sem(RAM=RAM2, data=df, cluster="clusterID", replace.constraints=TRUE),
                "level-2-only.*definition variable")

    ## Must NOT false-positive: a legitimate WITHIN-level ("level: 1")
    ## definition variable, with or without replace.constraints, is
    ## unaffected by this check.
    v <- runif(N, 0.3, 0.8)
    df3 <- data.frame(clusterID=clusterID, v=v, y=y)
    model3 <- "level: 1
  y ~~ data.v*y
level: 2
  y ~ 1
  y ~~ residB*y"
    expect_error(sem(RAM=lavaan2RAM(model3, obs.variables="y", std.lv=FALSE),
                     data=df3, cluster="clusterID"), NA)

    model4 <- "level: 1
  y ~~ residW*y
level: 2
  y ~ 1
  y ~~ residB*y
residW == exp(a0 + a1*data.z)"
    df4 <- data.frame(clusterID=clusterID, z=z, y=y)
    expect_error(sem(RAM=lavaan2RAM(model4, obs.variables="y", std.lv=FALSE),
                     data=df4, cluster="clusterID", replace.constraints=TRUE), NA)
})

test_that("sem() errors clearly on unsupported two-level input", {
    model <- "level: 1
  y ~ x
level: 2
  y ~ 1"
    RAM <- lavaan2RAM(model)

    set.seed(9)
    df <- data.frame(y=rnorm(20), x=rnorm(20), clusterID=rep(1:4, 5))

    ## Missing 'cluster'
    expect_error(sem(RAM=RAM, data=df))
    ## 'cluster' not a column in data
    expect_error(sem(RAM=RAM, data=df, cluster="nosuchcolumn"))
    ## replace.constraints=TRUE IS now supported for two-level RAM (see the
    ## dedicated tests below) -- with no "==" constraints in this model at
    ## all, it is simply a no-op and must NOT error.
    expect_error(sem(RAM=RAM, data=df, cluster="clusterID", replace.constraints=TRUE), NA)
    ## Cov/numObs-based input not yet supported for two-level RAM
    expect_error(sem(RAM=RAM, Cov=diag(2), numObs=100))
    ## Regression test: 'group' (a multi-group-only argument) used to be
    ## silently ignored for a two-level RAM instead of rejected -- a caller
    ## passing both 'cluster' and 'group' could be misled into thinking
    ## both dimensions were being modeled, when only 'cluster' was.
    expect_error(sem(RAM=RAM, data=df, cluster="clusterID", group="clusterID"),
                "group")
    ## ngroups > 1 has no effect on "level:" syntax and is rejected
    expect_error(lavaan2RAM(model, ngroups=2))

    ## Regression test: a missing cluster ID used to silently create a
    ## spurious "phantom" between-level row (unique(data[[cluster]])
    ## includes NA as if it were a real cluster), which rows with a missing
    ## cluster ID then joined to -- confirmed empirically to change fitted
    ## values, not just cosmetic. Now rejected explicitly.
    df.na.cluster <- df
    df.na.cluster$clusterID[1:2] <- NA
    expect_error(sem(RAM=RAM, data=df.na.cluster, cluster="clusterID"))
})

test_that("lavaan2RAM() rejects genuine level-2-only observed covariates", {
    ## Regression test for a HIGH-severity bug found on external review:
    ## the documentation already said genuine level-2-only observed
    ## covariates (e.g. "z" below, intended as a real cluster-level
    ## covariate with its own data column, as opposed to a level-1
    ## manifest variable's between-cluster component) are "not yet
    ## supported" -- but nothing actually enforced that. "z" was silently
    ## absorbed as an ordinary latent variable with no connection to any
    ## real data column at all: no error, no warning, just a free
    ## parameter estimated in a vacuum. Reproduced exactly as reported,
    ## then fixed by rejecting it explicitly.
    model <- "level: 1
  y ~~ y
level: 2
  y ~ z
  y ~~ y
  z ~~ z"
    expect_error(lavaan2RAM(model), "z")

    ## A genuine level-2-only latent factor indicator (never appearing at
    ## level 1 either) is the same underlying problem and is rejected too.
    model.q <- "level: 1
  y ~~ y
level: 2
  fb =~ y + q
  y ~~ y"
    expect_error(lavaan2RAM(model.q), "q")

    ## Legitimate models (every level-2 variable is either a shadow of a
    ## level-1 manifest, or a genuine latent factor) are unaffected.
    expect_error(lavaan2RAM("level: 1
  fw =~ y1+y2+y3
level: 2
  fb =~ y1+y2+y3"), NA)
})

test_that("vcov.mxsem(robust=TRUE) is guarded for two-level fits, matching summary()", {
    ## Regression test: summary.mxsem()'s two-level robust=TRUE guard was
    ## not mirrored in vcov.mxsem(), which called imxRobustSE() directly
    ## and hit the same underlying OpenMx limitation with an unguarded,
    ## cryptic error. An inconsistently-applied fix, found on external
    ## review, not a design gap.
    set.seed(17)
    nClusters <- 30; nPerCluster <- 5; N <- nClusters*nPerCluster
    clusterID <- rep(seq_len(nClusters), each=nPerCluster)
    x <- rnorm(N)
    y <- 0.5*x + rep(rnorm(nClusters,0,0.5), each=nPerCluster) + rnorm(N)
    df <- data.frame(clusterID=clusterID, x=x, y=y)
    RAM <- lavaan2RAM("level: 1\n y~b1*x\n y~~residW*y\nlevel: 2\n y~1\n y~~residB*y")
    fit <- sem(RAM=RAM, data=df, cluster="clusterID")
    expect_error(vcov(fit, robust=TRUE), "two-level")
})

test_that("sem() accepts character and numeric-double cluster ID columns", {
    ## Regression test for a medium-severity bug found on a follow-up
    ## external review: OpenMx's relational primaryKey/joinKey mechanism
    ## requires the key column to be an integer or a factor.
    ## .sem.twolevel() built between.data directly from data[[cluster]]
    ## with no type check, so a character cluster ID (e.g. "school-1") or
    ## a plain numeric/double one (as.numeric(1:10)) -- both very ordinary
    ## ways to code a cluster ID -- hit OpenMx's "primary key must be an
    ## integer or factor column" deep inside mxRun(). Confirmed on
    ## reproduction that this is worse than just an error: because it's a
    ## hard R error, sem()'s own mxRun-then-mxTryHard() fallback treats it
    ## as "maybe just needs a better start" and burns through all 50
    ## retries (several seconds) before finally re-raising it. Fixed by
    ## coercing an unsupported cluster column to a factor up front.
    set.seed(19)
    nClusters <- 20; nPerCluster <- 5; N <- nClusters*nPerCluster
    x <- rnorm(N)
    u <- rep(rnorm(nClusters, 0, 0.5), each=nPerCluster)
    y <- 0.5*x + u + rnorm(N)
    RAM <- lavaan2RAM("level: 1\n y~b1*x\n y~~residW*y\nlevel: 2\n y~1\n y~~residB*y")

    df.int <- data.frame(clusterID=rep(seq_len(nClusters), each=nPerCluster), x=x, y=y)
    df.char <- data.frame(clusterID=paste0("school-", rep(seq_len(nClusters), each=nPerCluster)), x=x, y=y)
    df.dbl <- data.frame(clusterID=as.numeric(rep(seq_len(nClusters), each=nPerCluster)), x=x, y=y)

    fit.int <- sem(RAM=RAM, data=df.int, cluster="clusterID")
    fit.char <- sem(RAM=RAM, data=df.char, cluster="clusterID")
    fit.dbl <- sem(RAM=RAM, data=df.dbl, cluster="clusterID")

    expect_identical(fit.int$mx.fit@output$status[[1]], 0L)
    expect_identical(fit.char$mx.fit@output$status[[1]], 0L)
    expect_identical(fit.dbl$mx.fit@output$status[[1]], 0L)
    expect_equal(unname(coef(fit.char)["b1"]), unname(coef(fit.int)["b1"]), tolerance=1e-6)
    expect_equal(unname(coef(fit.dbl)["b1"]), unname(coef(fit.int)["b1"]), tolerance=1e-6)
})

test_that("sem() supports LB intervals, run=FALSE, bounds, and missing data for two-level fits", {
    set.seed(13)
    nClusters <- 50; nPerCluster <- 6; N <- nClusters*nPerCluster
    clusterID <- rep(seq_len(nClusters), each=nPerCluster)
    u <- rep(rnorm(nClusters, 0, 0.6), each=nPerCluster)
    x <- rnorm(N); y <- 2*x + u + rnorm(N, sd=0.3)   # strong signal, true b1=2

    model <- "level: 1
  y ~ b1*x
  y ~~ residW*y
  x ~ meanX*1
  x ~~ varX*x
level: 2
  y ~ 1
  y ~~ residB*y"
    RAM <- lavaan2RAM(model)
    df <- data.frame(clusterID=clusterID, x=x, y=y)

    ## Likelihood-based CI
    fit.lb <- sem(RAM=RAM, data=df, cluster="clusterID", intervals.type="LB")
    s.lb <- summary(fit.lb)
    expect_true(s.lb$coefficients["b1","lbound"] < s.lb$coefficients["b1","ubound"])

    ## run=FALSE returns an unrun MxModel that can be run manually.
    fit.norun <- sem(RAM=RAM, data=df, cluster="clusterID", run=FALSE)
    expect_s4_class(fit.norun$mx.fit, "MxModel")
    fit.manual <- mxRun(fit.norun$mx.fit, silent=TRUE)
    expect_identical(fit.manual@output$status[[1]], 0L)

    ## ubound actually constrains the optimizer, not just accepted silently.
    fit.bound <- sem(RAM=RAM, data=df, cluster="clusterID", ubound=list(b1=1))
    expect_equal(unname(coef(fit.bound)["b1"]), 1, tolerance=1e-4)

    ## anova()/vcov() against a properly nested pair of fits (fixing the
    ## same free path at 0, not removing the variable from the model).
    model.null <- "level: 1
  y ~ 0*x
  y ~~ residW*y
  x ~ meanX*1
  x ~~ varX*x
level: 2
  y ~ 1
  y ~~ residB*y"
    fit.null <- sem(RAM=lavaan2RAM(model.null), data=df, cluster="clusterID")
    cmp <- anova(fit.bound, fit.null)
    expect_equal(cmp$diffdf[2], 1)
    expect_equal(dim(vcov(fit.bound)), rep(length(coef(fit.bound)), 2))

    ## robust=TRUE is rejected with a clear, two-level-specific message
    ## (OpenMx's imxRobustSE() cannot compute row-wise gradients through the
    ## relational "between" submodel -- see summary.mxsem()'s comment).
    expect_error(summary(fit.lb, robust=TRUE), "two-level")

    ## Missing individual-level data (FIML) does not error.
    df.na <- df
    df.na$y[sample(N, 20)] <- NA
    expect_error(sem(RAM=RAM, data=df.na, cluster="clusterID"), NA)
})

test_that("sem() fits a two-level structural model and handles partial level-2 declarations", {
    ## A regression between two GENUINE between-level latent factors (not
    ## just a random-intercept/CFA shadow structure). 3 indicators per
    ## factor (not 2): a 2-indicator between-level factor is only
    ## just-identified and, empirically, needs a much better-than-default
    ## starting point to converge on the first optimizer attempt (confirmed
    ## via mxTryHard(), which reaches the identical optimum regardless) --
    ## not a translation bug, just a harder optimization landscape than
    ## this test needs to exercise.
    set.seed(14)
    nClusters <- 60; nPerCluster <- 8; N <- nClusters*nPerCluster
    clusterID <- rep(seq_len(nClusters), each=nPerCluster)
    fb1.cl <- rnorm(nClusters, 0, 1)
    fb2.cl <- 0.6*fb1.cl + rnorm(nClusters, 0, 0.8)   # fb2 regressed on fb1
    fb1 <- rep(fb1.cl, each=nPerCluster); fb2 <- rep(fb2.cl, each=nPerCluster)
    fw1 <- rnorm(N, 0, 1); fw2 <- rnorm(N, 0, 1)
    y1 <- fb1+fw1+rnorm(N,0,.4); y2 <- fb1+fw1+rnorm(N,0,.4); y3 <- fb1+fw1+rnorm(N,0,.4)
    y4 <- fb2+fw2+rnorm(N,0,.4); y5 <- fb2+fw2+rnorm(N,0,.4); y6 <- fb2+fw2+rnorm(N,0,.4)
    df <- data.frame(clusterID=clusterID, y1=y1, y2=y2, y3=y3, y4=y4, y5=y5, y6=y6)

    model <- "level: 1
  fw1 =~ y1 + y2 + y3
  fw2 =~ y4 + y5 + y6
level: 2
  fb1 =~ y1 + y2 + y3
  fb2 =~ y4 + y5 + y6
  fb2 ~ fb1"
    RAM <- lavaan2RAM(model)
    fit <- sem(RAM=RAM, data=df, cluster="clusterID")
    expect_identical(fit$mx.fit@output$status[[1]], 0L)
    expect_equal(unname(coef(fit)["fb2ONfb1_between"]), 0.6, tolerance=0.3)

    ## Partial level-2 declaration: a level-1 manifest variable never
    ## mentioned in the "level: 2:" block gets no between-level shadow at
    ## all (matching lavaan's own semantics -- confirmed directly against
    ## lavaanify()'s output while investigating this), i.e. it is treated
    ## as purely within-level with no cluster random effect.
    model.partial <- "level: 1
  fw =~ y1 + y2 + y3
level: 2
  fb =~ y1 + y2"
    RAM.partial <- lavaan2RAM(model.partial)
    expect_false("y3" %in% colnames(RAM.partial$between$A))
    expect_true(all(c("y1","y2") %in% colnames(RAM.partial$between$A)))
})

test_that("plot.mxsem() supports multiple-group and two-level fits", {
    skip_if_not_installed("semPlot")

    ## Draw to a null device so this test has no filesystem side effect
    ## (no active device otherwise makes R auto-create ./Rplots.pdf).
    grDevices::pdf(nullfile())
    on.exit(grDevices::dev.off(), add=TRUE)

    ## Multiple-group: one diagram per group, picked via "group"
    set.seed(1)
    n1 <- 60; n2 <- 60
    x1 <- rnorm(n1); y1 <- 2 + 0.5*x1 + rnorm(n1)
    x2 <- rnorm(n2); y2 <- -1 + 0.5*x2 + rnorm(n2)
    dfmg <- rbind(data.frame(y=y1, x=x1, g="g1"), data.frame(y=y2, x=x2, g="g2"))
    RAMmg <- lavaan2RAM("y ~ c(b,b)*x + c(a1,a2)*1", ngroups=2)
    fitmg <- sem(RAM=RAMmg, data=dfmg, group="g")

    expect_silent(plot(fitmg, group="g1"))
    expect_silent(plot(fitmg, group="g2"))
    ## group set -> combine defaults to FALSE, single diagram
    expect_silent(plot(fitmg, group="g1", combine=FALSE))
    expect_error(plot(fitmg, group="nosuchgroup"))

    ## Two-level: one diagram per level, picked via "level". The CFA case
    ## exercises both a "shadow" variable (y1/y2/y3, also within-level
    ## manifests) and a genuine level-2-only latent factor (fb) at the
    ## between level -- semPlot::ramModel() errors on an all-latent model
    ## (no manifest variables), which the between level always is in this
    ## package's scope, so shadow variables are drawn as manifest for that
    ## plot specifically (see plot.mxsem()'s comments).
    set.seed(123)
    nClusters <- 40; nPerCluster <- 5
    clusterID <- rep(seq_len(nClusters), each=nPerCluster)
    fb <- rep(rnorm(nClusters, 0, 1), each=nPerCluster)
    fw <- rnorm(nClusters*nPerCluster, 0, 1)
    y1 <- fb + fw + rnorm(nClusters*nPerCluster, 0, 0.5)
    y2 <- 0.8*fb + 0.8*fw + rnorm(nClusters*nPerCluster, 0, 0.5)
    y3 <- 1.2*fb + 1.2*fw + rnorm(nClusters*nPerCluster, 0, 0.5)
    dftl <- data.frame(clusterID=clusterID, y1=y1, y2=y2, y3=y3)
    RAMtl <- lavaan2RAM("level: 1
  fw =~ y1 + y2 + y3
level: 2
  fb =~ y1 + y2 + y3")
    fittl <- sem(RAM=RAMtl, data=dftl, cluster="clusterID")

    expect_silent(plot(fittl, level="within"))
    expect_silent(plot(fittl, level="between"))
    ## level set -> combine defaults to FALSE, single diagram (level defaults
    ## to "within" whenever it is itself left NULL)
    expect_silent(plot(fittl, level="within", combine=FALSE))
    expect_error(plot(fittl, level="nosuchlevel"))

    ## Single-group/single-level plotting is unaffected by the new
    ## "group"/"level" arguments.
    RAMsg <- lavaan2RAM("y ~ x\ny ~ 1\nx ~ 1", obs.variables=c("y","x"))
    n <- 80; xx <- rnorm(n); yy <- 0.5*xx + rnorm(n)
    fitsg <- sem(RAM=RAMsg, data=data.frame(y=yy, x=xx))
    expect_silent(plot(fitsg))

    ## combine=TRUE: all groups (or both levels) as panels of one figure.
    ## plot.mxsem() restores par(mfrow=) via on.exit() once it returns, so
    ## only that its own two-panel layout doesn't error is checked here (not
    ## the transient mfrow value, which isn't observable after the call).
    op <- graphics::par(no.readonly=TRUE)
    expect_silent(plot(fitmg, combine=TRUE))
    expect_identical(graphics::par("mfrow"), op$mfrow)    # restored afterwards
    expect_silent(plot(fittl, combine=TRUE))
    expect_identical(graphics::par("mfrow"), op$mfrow)    # restored afterwards

    ## combine's default is missing(group) && missing(level): plot(fit) with
    ## neither set behaves exactly like the explicit combine=TRUE call above
    ## (all panels), while supplying either argument switches the default to
    ## FALSE (a single diagram) without needing combine=FALSE as well --
    ## already exercised by "group set"/"level set" cases above, so only the
    ## no-argument case needs checking here.
    expect_silent(plot(fitmg))
    expect_identical(graphics::par("mfrow"), op$mfrow)
    expect_silent(plot(fittl))
    expect_identical(graphics::par("mfrow"), op$mfrow)

    ## Explicit combine=TRUE still overrides a supplied group/level.
    expect_silent(plot(fitmg, group="g1", combine=TRUE))
    expect_identical(graphics::par("mfrow"), op$mfrow)
})

test_that("plot.mxsem() works when replace.constraints=TRUE has turned A/S/M into an mxAlgebra, in all three RAM dispatch paths", {
    ## Regression test for a high-severity bug (found on a careful re-
    ## review of the whole replace.constraints=TRUE feature): .plot.mxsem.
    ## one() decided which matrix-naming convention a submodel uses
    ## (single-/multiple-group's "Amatrix" etc. vs two-level's "A" etc.)
    ## by checking ONLY mx.sub@matrices for "Amatrix" -- but
    ## replace.constraints=TRUE can turn any of A/S/M into an mxAlgebra,
    ## which lives in @algebras, not @matrices. Two distinct failure
    ## modes resulted: (1) if "Amatrix" ITSELF was replaced, the
    ## convention check went to the wrong (two-level-style) branch
    ## entirely, looking for objects literally named "A" that don't exist
    ## in a single-/multiple-group model at all; (2) even when "Amatrix"
    ## was untouched, whichever OTHER matrix (S or M) was replaced still
    ## wasn't found (only @matrices was ever checked), returning NULL and
    ## erroring downstream at "seq_len(nrow(S))" ("argument must be
    ## coercible to non-negative integer"). A further wrinkle: an
    ## mxAlgebra's $result never carries dimnames the way a plain
    ## mxMatrix's $values does, so even after finding the right object,
    ## variable names had to be recovered from RAM's own (always-correct,
    ## substitution-unaffected) dimnames as a fallback.
    skip_if_not_installed("semPlot")
    grDevices::pdf(nullfile())
    on.exit(grDevices::dev.off(), add=TRUE)

    ## -- Single-group: S replaced, A untouched (the specific combination
    ## that used to slip past the (correct, by luck) convention check but
    ## still fail to locate "Smatrix"). --
    set.seed(1)
    df.sg <- data.frame(y=rnorm(400, 2, 1))
    model.sg <- "y ~ mu*1
  y ~~ sigma*y
  sigma == exp(a0)"
    fit.sg <- sem(RAM=lavaan2RAM(model.sg, obs.variables="y"), data=df.sg,
                 replace.constraints=TRUE)
    expect_error(plot(fit.sg), NA)

    ## -- Multiple-group: one group's S replaced, the other's untouched;
    ## check both the single-diagram and combine=TRUE (panelled) paths. --
    set.seed(1)
    n1 <- 200; n2 <- 220
    x1 <- rnorm(n1); y1 <- 2 + 0.5*x1 + rnorm(n1, sd=1)
    x2 <- rnorm(n2); y2 <- -1 + 0.5*x2 + rnorm(n2, sd=1.3)
    df.mg <- rbind(data.frame(y=y1, x=x1, g="g1"), data.frame(y=y2, x=x2, g="g2"))
    model.mg <- "y ~ c(b,b)*x + c(a1,a2)*1
  y ~~ c(sigma1,sigma2)*y
  sigma1 == exp(a0)"
    fit.mg <- sem(RAM=lavaan2RAM(model.mg, ngroups=2), data=df.mg, group="g",
                 replace.constraints=TRUE)
    expect_error(plot(fit.mg, group="g1"), NA)
    expect_error(plot(fit.mg, group="g2"), NA)
    expect_error(plot(fit.mg, combine=TRUE), NA)

    ## -- Two-level: within's S replaced, between untouched; both the
    ## single-diagram and combine=TRUE (panelled) paths. --
    set.seed(123)
    nClusters <- 60; nPerCluster <- 6
    clusterID <- rep(seq_len(nClusters), each=nPerCluster)
    u <- rep(rnorm(nClusters, 0, 0.7), each=nPerCluster)
    df.tl <- data.frame(clusterID=clusterID, y=1.5+u+rnorm(nClusters*nPerCluster))
    model.tl <- "level: 1
  y ~~ residW*y
level: 2
  y ~ 1
  y ~~ residB*y
residW == exp(a0)"
    fit.tl <- sem(RAM=lavaan2RAM(model.tl), data=df.tl, cluster="clusterID",
                 replace.constraints=TRUE)
    expect_error(plot(fit.tl, level="within"), NA)
    expect_error(plot(fit.tl, level="between"), NA)
    expect_error(plot(fit.tl, combine=TRUE), NA)

    ## -- The extracted matrix value must be numerically correct (not just
    ## "doesn't error") and carry the right dimnames.
    mx.sub <- fit.sg$mx.fit
    S <- if ("Smatrix" %in% names(mx.sub@matrices)) mx.sub@matrices$Smatrix$values else {
      v <- mx.sub@algebras$Smatrix$result
      dimnames(v) <- dimnames(fit.sg$RAM$S)
      v
    }
    expect_equal(unname(S[1,1]), unname(exp(coef(fit.sg)["a0"])), tolerance=1e-6)
    expect_identical(dimnames(S), dimnames(fit.sg$RAM$S))
})

context("Checking metaFIML functions")

test_that("metaFIML() works correctly", {

    ## Univariate meta-analysis without AV
    fit1a <- metaFIML(y=r, v=r_v, x=JP_alpha, data=Jaramillo05)

    ## r/JP_alpha are single-indicator markers for the latent fy/fx (loading
    ## fixed at 1), with the actual free mean modelled entirely through
    ## fy/fx's own "~1" below -- so r's and JP_alpha's OWN observed
    ## intercepts must be explicitly fixed at 0 (matching metaFIML()'s own
    ## internal Mmatrix, which hard-fixes them: R/metaFIML.R's "Mx" is a
    ## literal matrix(0, ...), never a free parameter). Previously this
    ## fixing happened implicitly (lavaan2RAM() silently let lavaanify's
    ## own auto-generated, still-fixed "~1" row stand in for it), but that
    ## was itself a bug (see lavaan2RAM()'s R/lavaan2RAM.R comments) fixed
    ## to make observed variables free-by-default like lavaan::sem()
    ## itself -- so what used to happen by accident must now be requested
    ## explicitly here, exactly like metaFIML() itself requests it.
    m1 <- "fy =~ 1*r
           r ~ 0*1
           r ~~ data.r_v*r
           fx =~ 1*JP_alpha
           JP_alpha ~ 0*1
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

    ## Same reasoning as m1 above: r/JP_alpha/IDV are single-indicator
    ## markers, so their own observed intercepts must be explicitly fixed
    ## at 0 to match metaFIML()'s own internal Mmatrix.
    m2 <- "fy =~ 1*r
           r ~ 0*1
           r ~~ data.r_v*r
           fx =~ 1*JP_alpha
           JP_alpha ~ 0*1
           JP_alpha ~~ 0*JP_alpha
           fy ~ Slope1_1*fx
           fy ~~ Tau2_1_1*fy
           fx ~~ CovX1_X1*fx
           fx ~ MeanX1*1
           fy ~ Intercept1*1

           fz =~ 1*IDV
           IDV ~ 0*1
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

    ## Same reasoning as m1/m2 above: lifesat/lifecon/gnp are
    ## single-indicator markers, so their own observed intercepts must be
    ## explicitly fixed at 0 to match metaFIML()'s own internal Mmatrix.
    m3 <- "fy1 =~ 1*lifesat
           lifesat ~ 0*1
           lifesat ~~ data.lifesat_var*lifesat
           fy2 =~ 1*lifecon
           lifecon ~ 0*1
           lifecon ~~ data.lifecon_var*lifecon
           lifesat ~~ data.inter_cov*lifecon

           fx =~ 1*gnp
           gnp ~ 0*1
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

context("Checking sem() two-level results against meta3L()")

## meta3L()'s three-level random-effects model
##   y_ij = b0 + b1*x_ij + u_(2)ij + u_(3)j + e_ij,  u_(2)~N(0,Tau2_2), u_(3)~N(0,Tau2_3)
## is the same model as a two-level sem() fit with x/y as level-1 (within)
## manifests and the known sampling variance "v" fixed via a definition
## variable, i.e. Tau2_2 = the level-1 latent's own variance ("tauW" below)
## and Tau2_3 = the level-2 (between) residual variance of y ("tauB" below).
## meta3L() itself works in wide format (one row per cluster, columns
## y_1..y_k) via a single mxExpectationNormal(); sem()'s two-level path
## works in long format (one row per effect size) via OpenMx's relational
## SEM (primaryKey/joinKey) -- two different OpenMx mechanisms computing,
## in principle, the same joint likelihood, so this is a genuine
## cross-implementation check, not just a restatement of one in the other.
test_that("sem() two-level matches meta3L() for intercept-only models", {
    m1 <- meta3L(y=logOR, v=v, cluster=Cluster, data=Bornmann07)
    model <- "level: 1
                fw =~ 1*logOR
                fw ~~ tauW*fw
                logOR ~~ data.v*logOR
              level: 2
                logOR ~ b0*1
                logOR ~~ tauB*logOR"
    s1 <- sem(RAM=lavaan2RAM(model, std.lv=FALSE), data=Bornmann07, cluster="Cluster")

    co3 <- summary(m1)$coefficients
    cosem <- summary(s1)$coefficients
    expect_equal(unname(co3["Intercept", "Estimate"]), unname(cosem["b0", "Estimate"]), tolerance=1e-6)
    expect_equal(unname(co3["Intercept", "Std.Error"]), unname(cosem["b0", "Std.Error"]), tolerance=1e-6)
    expect_equal(unname(co3["Tau2_2", "Estimate"]), unname(cosem["tauW", "Estimate"]), tolerance=1e-6)
    expect_equal(unname(co3["Tau2_3", "Estimate"]), unname(cosem["tauB", "Estimate"]), tolerance=1e-6)
    expect_equal(m1$mx.fit@output$Minus2LogLikelihood,
                 s1$mx.fit@output$Minus2LogLikelihood, tolerance=1e-6)

    m4 <- meta3L(y=y, v=v, cluster=District, data=Cooper03)
    model <- "level: 1
                fw =~ 1*y
                fw ~~ tauW*fw
                y ~~ data.v*y
              level: 2
                y ~ b0*1
                y ~~ tauB*y"
    s4 <- sem(RAM=lavaan2RAM(model, std.lv=FALSE), data=Cooper03, cluster="District")

    co3 <- summary(m4)$coefficients
    cosem <- summary(s4)$coefficients
    expect_equal(unname(co3["Intercept", "Estimate"]), unname(cosem["b0", "Estimate"]), tolerance=1e-6)
    expect_equal(unname(co3["Tau2_2", "Estimate"]), unname(cosem["tauW", "Estimate"]), tolerance=1e-6)
    expect_equal(unname(co3["Tau2_3", "Estimate"]), unname(cosem["tauB", "Estimate"]), tolerance=1e-6)
    expect_equal(m4$mx.fit@output$Minus2LogLikelihood,
                 s4$mx.fit@output$Minus2LogLikelihood, tolerance=1e-6)
})

test_that("sem() two-level matches meta3L() for a predictor, up to an analytic -2LL offset", {
    ## sem()'s within-level fixes an exogenous covariate's own mean at 0 by
    ## default (see .lavaan2RAM.twolevel()), so it jointly models the
    ## covariate's own distribution -- specifically its uncentered second
    ## moment mean(x^2), since the mean is fixed at 0 rather than free --
    ## as part of the likelihood. meta3L() never does this: it conditions
    ## on x rather than modelling its distribution at all. The regression
    ## and heterogeneity parameters themselves are unaffected (their
    ## information block is independent of x's own free-variance
    ## parameter), but -2LL differs by exactly N*(log(2*pi*mean(x^2))+1),
    ## the -2LL of a scalar normal MLE with its mean fixed at 0. This is
    ## NOT specific to an uncentered predictor -- mean-centering x first
    ## (as here, via scale(..., scale=FALSE)) does not make the offset
    ## vanish, it only means mean(x^2) and the ordinary sample variance
    ## coincide (since mean(x)=0 either way); the offset itself is still
    ## large and must still be subtracted, exactly as for the deliberately
    ## uncentered predictor in the next test.
    Cooper03$YearC <- scale(Cooper03$Year, scale=FALSE)[, 1]
    m5 <- meta3L(y=y, v=v, cluster=District, x=YearC, data=Cooper03)
    model <- "level: 1
                fw =~ 1*y
                fw ~~ tauW*fw
                y ~~ data.v*y
                y ~ b1*YearC
              level: 2
                y ~ b0*1
                y ~~ tauB*y"
    s5 <- sem(RAM=lavaan2RAM(model, std.lv=FALSE), data=Cooper03, cluster="District")

    co3 <- summary(m5)$coefficients
    cosem <- summary(s5)$coefficients
    expect_equal(unname(co3["Intercept", "Estimate"]), unname(cosem["b0", "Estimate"]), tolerance=1e-6)
    expect_equal(unname(co3["Slope_1", "Estimate"]), unname(cosem["b1", "Estimate"]), tolerance=1e-6)
    expect_equal(unname(co3["Slope_1", "Std.Error"]), unname(cosem["b1", "Std.Error"]), tolerance=1e-6)
    expect_equal(unname(co3["Tau2_2", "Estimate"]), unname(cosem["tauW", "Estimate"]), tolerance=1e-6)
    expect_equal(unname(co3["Tau2_3", "Estimate"]), unname(cosem["tauB", "Estimate"]), tolerance=1e-6)

    offset <- nrow(Cooper03) * (log(2*pi*mean(Cooper03$YearC^2)) + 1)
    expect_equal(m5$mx.fit@output$Minus2LogLikelihood,
                 s5$mx.fit@output$Minus2LogLikelihood - offset, tolerance=1e-6)
})

test_that("sem() two-level matches meta3L() for an uncentered predictor, up to an analytic -2LL offset", {
    ## Same offset as above (N*(log(2*pi*mean(x^2))+1)), but now with a
    ## genuinely uncentered 0/1 predictor, where mean(x^2) and the ordinary
    ## sample variance p*(1-p) actually differ -- confirming the offset
    ## formula must use the UNCENTERED second moment (matching the
    ## within-level mean fixed at 0), not the sample variance about x's own
    ## mean, regardless of whether x happens to be centered.
    Bornmann07$TypeNum <- as.numeric(Bornmann07$Type) - 1
    m2 <- meta3L(y=logOR, v=v, x=TypeNum, cluster=Cluster, data=Bornmann07)
    model <- "level: 1
                fw =~ 1*logOR
                fw ~~ tauW*fw
                logOR ~~ data.v*logOR
                logOR ~ b1*TypeNum
              level: 2
                logOR ~ b0*1
                logOR ~~ tauB*logOR"
    s2 <- sem(RAM=lavaan2RAM(model, std.lv=FALSE), data=Bornmann07, cluster="Cluster")

    co3 <- summary(m2)$coefficients
    cosem <- summary(s2)$coefficients
    expect_equal(unname(co3["Intercept", "Estimate"]), unname(cosem["b0", "Estimate"]), tolerance=1e-6)
    expect_equal(unname(co3["Slope_1", "Estimate"]), unname(cosem["b1", "Estimate"]), tolerance=1e-6)
    expect_equal(unname(co3["Slope_1", "Std.Error"]), unname(cosem["b1", "Std.Error"]), tolerance=1e-6)
    expect_equal(unname(co3["Tau2_2", "Estimate"]), unname(cosem["tauW", "Estimate"]), tolerance=1e-6)
    expect_equal(unname(co3["Tau2_3", "Estimate"]), unname(cosem["tauB", "Estimate"]), tolerance=1e-6)

    offset <- nrow(Bornmann07) * (log(2*pi*mean(Bornmann07$TypeNum^2)) + 1)
    expect_equal(m2$mx.fit@output$Minus2LogLikelihood,
                 s2$mx.fit@output$Minus2LogLikelihood - offset, tolerance=1e-6)
})

test_that("sem() two-level matches meta3L() when a label is shared across levels", {
    ## meta3L()'s RE2.constraints/RE3.constraints="0.2*EqTau2" uses
    ## metaSEM's own as.mxMatrix() "start*label" convention -- the "0.2*"
    ## start-value prefix is NOT valid lavaan syntax and must be dropped
    ## when translating to lavaan2RAM() (an earlier attempt at
    ## "0.2*EqTau2*fw" silently parsed as a FIXED value of 0.2 with no
    ## label at all, rather than erroring, so this is worth a regression
    ## test in its own right). The same label "EqTau2" used in both the
    ## level-1 and level-2 blocks equates the two variance components
    ## across the within/between submodels, exactly like meta3L()'s shared
    ## label equates Tau2_2 and Tau2_3.
    m6 <- meta3L(y=y, v=v, cluster=District, data=Cooper03,
                RE2.constraints="0.2*EqTau2", RE3.constraints="0.2*EqTau2")
    model <- "level: 1
                fw =~ 1*y
                fw ~~ EqTau2*fw
                y ~~ data.v*y
              level: 2
                y ~ b0*1
                y ~~ EqTau2*y"
    s6 <- sem(RAM=lavaan2RAM(model, std.lv=FALSE), data=Cooper03, cluster="District")

    expect_equal(length(coef(s6)), 2L)   # b0 and EqTau2, properly equated across levels

    co3 <- summary(m6)$coefficients
    cosem <- summary(s6)$coefficients
    expect_equal(unname(co3["Intercept", "Estimate"]), unname(cosem["b0", "Estimate"]), tolerance=1e-6)
    expect_equal(unname(co3["EqTau2", "Estimate"]), unname(cosem["EqTau2", "Estimate"]), tolerance=1e-6)
    expect_equal(unname(co3["EqTau2", "Std.Error"]), unname(cosem["EqTau2", "Std.Error"]), tolerance=1e-6)
    expect_equal(m6$mx.fit@output$Minus2LogLikelihood,
                 s6$mx.fit@output$Minus2LogLikelihood, tolerance=1e-6)
})
