test_that("FEM baseline likelihoods and RMSEA match normal-ML reference models", {
  nm <- c("x", "y")
  make <- function(a,b,c) matrix(c(a,b,b,c), 2, dimnames=list(nm,nm))
  C <- list(make(1,.1,1), make(16,.2,2), make(2,.4,3))
  n <- c(30,40,60)
  fit <- tssem1(C,n,method="FEM",cor.analysis=FALSE)
  expect_equal(fit$mx.fit@output$status$code,0)
  expect_gt(min(eigen(vcov(fit),symmetric=TRUE,only.values=TRUE)$values),0)
  baseline <- function(type) {
    models <- lapply(seq_along(C), function(i) {
      m <- OpenMx::mxMatrix(type,nrow=2,ncol=2,free=TRUE,
                            values=if(type == "Symm") OpenMx::vech(C[[i]]) else diag(C[[i]]),name="S")
      OpenMx::mxModel(paste0("g",i),m,OpenMx::mxData(C[[i]],type="cov",numObs=n[i]),
                      OpenMx::mxExpectationNormal("S",dimnames=nm),OpenMx::mxFitFunctionML())
    })
    model <- do.call(OpenMx::mxModel,c(list("baseline"),models,
                     list(OpenMx::mxFitFunctionMultigroup(paste0("g",seq_along(C))))))
    OpenMx::mxRun(model,silent=TRUE)
  }
  sat <- baseline("Symm"); ind <- baseline("Diag")
  expect_equal(sat@output$status$code,0)
  expect_equal(ind@output$status$code,0)
  expect_equal(unname(fit$baseMinus2LL[1:2]),c(sat@output$fit,ind@output$fit),tolerance=1e-6)
  ref <- lavaan::lavaan('x ~~ c(vx,vx,vx)*x\ny ~~ c(vy,vy,vy)*y\nx ~~ c(cv,cv,cv)*y',
                        sample.cov=C,sample.nobs=n,meanstructure=FALSE,auto.var=FALSE)
  expect_true(lavaan::lavInspect(ref,"converged"))
  expect_true(lavaan::lavInspect(ref,"post.check"))
  expect_equal(as.numeric(coef(fit)),as.numeric(lavaan::fitted(ref)[[1]]$cov),tolerance=1e-5)
  stat <- summary(fit)$stat
  measures <- lavaan::fitMeasures(ref,c("chisq","df","rmsea","rmsea.ci.lower","rmsea.ci.upper"))
  expect_equal(unname(stat[c("Chi-square of target model","DF of target model","RMSEA",
                              "RMSEA lower 95% CI","RMSEA upper 95% CI"),1]),
               as.numeric(measures),tolerance=1e-5)
  single <- tssem1(list(C[[1]]),n[1],method="FEM",cor.analysis=FALSE)
  expect_equal(summary(single)$stat["Chi-square of target model",1],0,tolerance=1e-6)
})

test_that("FEM covariance analysis retains variance-only observations", {
  nm <- c("x","y")
  complete <- matrix(c(2,.4,.4,3),2,dimnames=list(nm,nm))
  partial <- matrix(c(4,NA,NA,NA),2,dimnames=list(nm,nm))
  fit <- tssem1(list(complete,partial),c(60,40),method="FEM",cor.analysis=FALSE)
  expect_equal(fit$mx.fit@output$status$code,0)
  expect_equal(fit$data[[2]][1,1],4)
  expect_true(all(is.finite(fit$baseMinus2LL)))
  expect_true(is.finite(summary(fit)$stat["SRMR",1]))
  ## The normal likelihood estimates the variance from both studies.
  expect_equal(coef(fit)[1,1],(59*2+39*4)/100,tolerance=1e-5)
  one <- tssem1(list(matrix(2,1,1)),40,method="FEM",cor.analysis=FALSE)
  expect_equal(one$mx.fit@output$status$code,0)
  expect_equal(summary(one)$stat["Chi-square of target model",1],0,tolerance=1e-6)
  expect_error(tssem1(list(partial),40,method="FEM",cor.analysis=TRUE),"no observed correlations")
  expect_error(tssem1(list(matrix(NA_real_,2,2)),40,method="FEM",cor.analysis=FALSE),"no observed variances")
})

test_that("positive definiteness handles one-sided names and empty lists", {
  x <- diag(2); rownames(x) <- c("x","y")
  expect_true(is.pd(x))
  expect_identical(is.pd(list()),logical(0))
  colnames(x) <- c("y","x")
  expect_false(is.pd(x))
})

test_that("RMSEA uses the same normal-ML scaling for estimates and limits", {
  for (groups in 1:3) {
    result <- .fitIndices(tT=20,dfT=5,tB=40,dfB=10,N=150,G=groups)
    expect_equal(result$RMSEA,sqrt(groups*(20-5)/(150*5)))
    ncp <- vapply(c(.95,.05),function(p)
      uniroot(function(x) pchisq(20,df=5,ncp=x)-p,c(0,200))$root,numeric(1))
    expect_equal(unname(result$RMSEA.CI),sqrt(groups*ncp/(150*5)),tolerance=1e-6)
  }
})

test_that("meta3L I2 results preserve their requested order and values", {
  set.seed(823)
  id <- rep(1:20,each=4)
  y <- rnorm(20,sd=.4)[id]+rnorm(80,sd=.5)
  v <- seq(.01,.1,length.out=80)
  fit <- meta3L(y=y,v=v,cluster=id,I2=c("I2am","I2q"),intervals.type="z")
  expect_true(fit$mx.fit@output$status$code %in% c(0,1))
  values <- summary(fit)$I2.values
  labels <- c("I2_2 (Typical v: arithmetic mean)","I2_2 (Typical v: Q statistic)",
              "I2_3 (Typical v: arithmetic mean)","I2_3 (Typical v: Q statistic)")
  expect_identical(rownames(values),labels)
  expected <- vapply(c("I2am_2","I2q_2","I2am_3","I2q_3"),
                     function(x) as.numeric(OpenMx::mxEvalByName(x,fit$mx.fit)),numeric(1))
  expect_equal(as.numeric(values[,"Estimate"]),unname(expected),tolerance=1e-6)
})
