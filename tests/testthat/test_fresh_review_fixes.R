test_that("covariance moderators preserve variances throughout mediated paths", {
  nm <- c("x1", "x2", "m", "y", "w")
  A <- matrix("0", 5, 5, dimnames=list(nm,nm))
  A["m","x1"] <- "0.45*a"; A["m","x2"] <- "0.45*b"
  A["y","m"] <- "0.8*c"; A["w","y"] <- "0.7*d"
  S <- matrix("0", 5, 5, dimnames=list(nm,nm))
  diag(S) <- c("1", "1", "0.5*em", "0.5*ey", "0.5*ew")
  S["x1","x2"] <- S["x2","x1"] <- "0.2*r"
  Sx <- matrix(0, 5, 5, dimnames=list(nm,nm)); Sx["x1","x2"] <- Sx["x2","x1"] <- 2
  model <- OpenMx::mxModel("check", create.vechsR(A,S,Sx=Sx),
      OpenMx::mxAlgebraFromString("sum(vechsR)",name="myObjective"),
      OpenMx::mxFitFunctionAlgebra("myObjective"))
  fit <- OpenMx::mxRun(model, useOptimizer=FALSE, silent=TRUE)
  actual <- OpenMx::mxEvalByName("impliedR", fit, compute=TRUE)
  expect_equal(unname(diag(actual)), rep(1,5), tolerance=1e-10)
  expect_equal(unname(actual["y","w"]), .7, tolerance=1e-10)
})

test_that("OSMASEM fixed parameters do not depend on their labels", {
  mk <- function(r) matrix(c(1,r,r,1),2,dimnames=list(c("x","y"),c("x","y")))
  data <- Cor2DataFrame(lapply(c(.2,.3,.4,.35),mk), n=rep(100,4))
  for (label in c("beta","Cor_beta","Tau1_beta")) {
    R <- lavaan2RAM(paste0("y ~ ",label,"*x\nx ~~ 1*x"),obs.variables=c("x","y"))
    fit <- osmasem(RAM=R,data=data,RE.type="Zero")
    expect_equal(unname(coef(fit)[label]), .3125, tolerance=1e-5)
    expect_equal(rownames(vcov(fit)), label)
    expect_length(suppressWarnings(coef(fit,select="random")),0)
  }
})

test_that("legacy OSMASEM2 interval requests are not silently ignored", {
  expect_error(osmasem2(RAM=NULL,data=list(),intervals.type="LB",run=FALSE), "no longer supports")
  expect_error(osmasem2("m",NULL,list(),TRUE,"Diag","Symm","Zero",FALSE,"LB",run=FALSE), "positional")
})

test_that("correlation WLS handles no directed paths", {
  C <- matrix(c(1,.1,.3,.1,1,.6,.3,.6,1),3,dimnames=list(c("x","y","z"),c("x","y","z")))
  R <- lavaan2RAM("x ~~ 1*x\ny ~~ 1*y\nz ~~ 1*z\nx ~~ r*y\nx ~~ r*z\ny ~~ r*z",obs.variables=c("x","y","z"))
  for (dc in c(FALSE,TRUE)) {
    fit <- wls(C,diag(3)*.01,n=500,RAM=R,diag.constraints=dc)
    expect_true(fit$mx.fit@output$status$code %in% c(0,1))
    expect_equal(unname(coef(fit)["r"]), mean(c(.1,.3,.6)), tolerance=1e-5)
  }
})

test_that("WLS aligns separately supplied F and sampling covariance names", {
  vars <- c("x","y","z")
  C <- matrix(c(1,.1,.3,.1,1,.6,.3,.6,1),3,dimnames=list(vars,vars))
  R <- lavaan2RAM("y ~ b*x\nz ~ b*x\nx ~~ 1*x",obs.variables=vars)
  V <- diag(c(.01,.04,.09));dimnames(V) <- list(c("y_x","z_x","z_y"),c("y_x","z_x","z_y"))
  f <- function(F,V) wls(C,V,n=500,Amatrix=R$A,Smatrix=R$S,Fmatrix=F)
  a <- f(R$F,V); b <- f(R$F[3:1,,drop=FALSE],V[3:1,3:1])
  expect_equal(coef(b),coef(a),tolerance=1e-6)
  expect_equal(b$mx.fit@output$fit,a$mx.fit@output$fit,tolerance=1e-6)
})

test_that("asyCov aligns studies and validates sample sizes", {
  C <- matrix(c(1,.1,.3,.1,1,.6,.3,.6,1),3,dimnames=list(c("x","y","z"),c("x","y","z")))
  D <- C[c("z","x","y"),c("z","x","y")]
  for (mode in c("individual","weighted","unweighted")) {
    a <- asyCov(list(C,C),c(100,100),acov=mode)
    b <- asyCov(list(C,D),c(100,100),acov=mode)
    expect_equal(unname(b),unname(a),tolerance=1e-10)
  }
  expect_error(asyCov(list(C,C),c(50,100,200,400),acov="individual"), "Lengths")
  expect_error(asyCov(C,-100), "positive finite")
})

test_that("missing sampling covariances are diagnosed", {
  y <- cbind(c(.1,.4,.2),c(.3,.6,.8));v <- cbind(rep(.1,3),rep(0,3),rep(.2,3));v[1,2] <- NA
  expect_error(meta(y=y,v=v,run=FALSE), "Missing sampling covariance")
  y[1,2] <- NA
  expect_no_error(meta(y=y,v=v,run=FALSE))
})

test_that("generated correlation labels remain unique with underscores", {
  nm <- c("z","y_z","x_y","x");S <- diag(4);dimnames(S)<-list(nm,nm)
  S1<-S;S2<-S;S1[3,1]<-S1[1,3]<-NA;S2[3,1]<-S2[1,3]<-NA
  d <- Cor2DataFrame(list(S1,S2,S),n=c(100,120,150))
  expect_false(anyDuplicated(d$ylabels)>0)
  expect_false(anyDuplicated(colnames(d$data))>0)
  expect_equal(unname(diag(pattern.na(d,type="osmasem",show.na=FALSE))),
               unname(colSums(!is.na(d$data[,seq_along(d$ylabels),drop=FALSE]))))
  expect_no_error(.subset.Cor2DataFrame(d,c("z","x_y","x"),TRUE))
})

test_that("generated labels do not impose unintended equalities", {
  expect_warning(R <- lavaan2RAM("y ~ x\nz ~ yONx*x"), "conflicts")
  expect_false(identical(R$A["y","x"],R$A["z","x"]))
  expect_warning(R <- lavaan2RAM("aONb ~ c\na ~ bONc"), "conflicts")
  expect_false(identical(R$A["aONb","c"],R$A["a","bONc"]))
})

test_that("SMD bias adjustment follows the actual standardizer", {
  d <- smdMTS(c(0,1),c(1,4),c(5,100),homogeneity="none",bias.adjust=FALSE)$y
  g <- smdMTS(c(0,1),c(1,4),c(5,100),homogeneity="none",bias.adjust=TRUE)$y
  J <- exp(lgamma(2)-.5*log(2)-lgamma(1.5))
  expect_equal(unname(g/d),J,tolerance=1e-6)
  d <- smdMES(0,1,matrix(1),matrix(4),5,100,homogeneity="none",bias.adjust=FALSE)$y
  g <- smdMES(0,1,matrix(1),matrix(4),5,100,homogeneity="none",bias.adjust=TRUE)$y
  expect_equal(unname(g/d),J,tolerance=1e-6)
  set.seed(397)
  n <- c(5,40,100); B <- 200000
  q <- (n[1]*rchisq(B,n[1]-1)/(n[1]-1) + n[2]*rchisq(B,n[2]-1)/(n[2]-1) +
        n[3]*rchisq(B,n[3]-1)/(n[3]-1))/sum(n)
  expect_equal(.smd.bias.factor(n,pooled=TRUE), 1/mean(1/sqrt(q)), tolerance=.001)
  expect_error(.smd.bias.factor(c(2,100)), "more than one degree")
})

test_that("substitution preserves bounds and invalid bounds are rejected", {
  set.seed(146);x<-rnorm(200);d<-data.frame(x=x,y=-.5*x+rnorm(200))
  R <- lavaan2RAM("y ~ a*x\na == -0.5\na > 0.1")
  expect_error(sem(RAM=R,data=d,replace.constraints=TRUE,run=FALSE), "violates its lower bound")
  R <- lavaan2RAM("y ~ a*x")
  expect_error(sem(RAM=R,data=d,lbound=list(aa=.1),run=FALSE), "not found")
  expect_error(sem(RAM=R,data=d,lbound=list(.1),run=FALSE), "parameter names")
})

test_that("substituted parameter bounds agree with an independent constrained fit", {
  set.seed(145); n <- 600; x <- rnorm(n)
  d <- data.frame(x=x,y=-.6*x+rnorm(n),z=.6*x+rnorm(n))
  model <- 'y ~ a*x\nz ~ b*x\na == -b\na > 0.1'
  ram <- lavaan2RAM(model,auto.cov.y=TRUE)
  fit <- sem(RAM=ram,data=d,replace.constraints=TRUE)
  reference <- lavaan::sem(model,data=d,fixed.x=FALSE,meanstructure=TRUE,
                           auto.cov.y=TRUE,information='observed')
  expect_true(fit$mx.fit@output$status$code %in% c(0,1))
  expect_true(lavaan::lavInspect(reference,'converged'))
  expect_equal(unname(-coef(fit)['b']),unname(lavaan::coef(reference)['a']),tolerance=5e-5)
  expect_equal(fit$mx.fit@output$fit,-2*as.numeric(logLik(reference)),tolerance=1e-4)
})

test_that("explicit mean labels do not collide with unused default labels", {
  expect_no_warning(ram <- lavaan2RAM('y ~ ymean*1\ny ~ x'))
  expect_match(ram$M[1,'y'], 'ymean')
})

test_that("ambiguous old OSMASEM datasets fail before model construction", {
  d <- list(obslabels=c('x','y','z'), ylabels=c('y_x','y_x'),
            data=data.frame(y_x=1, check.names=FALSE))
  expect_error(osmasem2(RAM=list(),data=d,run=FALSE),'duplicate statistic')
})

test_that("legacy sampling covariance implementation validates inputs", {
  C <- matrix(c(1,.2,.3,.2,1,.4,.3,.4,1),3,dimnames=list(letters[1:3],letters[1:3]))
  expect_error(asyCovOld(list(C,C),c(10,20,30),run=FALSE),'Lengths')
  expect_error(asyCovOld(C,NA_real_,run=FALSE),'positive finite')
  a <- asyCovOld(list(C,C),c(100,100),as.matrix=FALSE)
  b <- asyCovOld(list(C,C[3:1,3:1]),c(100,100),as.matrix=FALSE)
  ## asyCovOld fits an OpenMx model for each matrix.
  expect_equal(a,b,tolerance=5e-6)
})

test_that("direct WLS node names cannot silently disagree", {
  nm <- c('x','y','z'); C <- diag(3); dimnames(C) <- list(nm,nm)
  R <- lavaan2RAM('y ~ x\nz ~ x',obs.variables=nm)
  S <- R$S[3:1,3:1,drop=FALSE]
  expect_error(wls(C,diag(3),100,Amatrix=R$A,Smatrix=S,Fmatrix=R$F,run=FALSE),
               'node dimensions')
})

test_that("FIML diagnoses missing covariances between observed outcomes", {
  y <- cbind(c(.1,.4,.2),c(.3,.6,.8))
  v <- cbind(rep(.1,3),rep(0,3),rep(.2,3)); v[1,2] <- NA
  x <- c(1,2,3)
  expect_error(metaFIML(y=y,v=v,x=x,run=FALSE),'Missing sampling covariance')
  y[1,2] <- NA
  expect_no_error(metaFIML(y=y,v=v,x=x,run=FALSE))
})

test_that("SEC warnings refer to weighted study contributions", {
  X <- matrix(1,2,2); weights <- diag(2)
  fit <- SEC(X,method='weighted',q=c(1,1),weights=weights)
  expect_equal(unname(fit$R),diag(2),tolerance=1e-12)
  X <- matrix(c(1,1,1,0,1,0,1,0,1),3,byrow=TRUE)
  W <- matrix(0,3,3); W[1,] <- 1
  expect_warning(SEC(X,method='weighted',q=rep(1,3),weights=W),'study contributions')
})

test_that("a zero-loading indicator does not anchor a latent mean", {
  m <- 'level: 1\nfw =~ y1 + y2 + y3\nlevel: 2\nfb =~ 0*y1 + y2 + y3\nfb ~ mb*1\ny1 ~ 0*1'
  expect_warning(lavaan2RAM(m), 'between-level mean.*may be unidentified')
})


test_that("meta-analysis accepts dotted R object names with valid columns", {
  effect.data <- matrix(c(.1,.2,.3,.1,.1,.1),3,2,
                        dimnames=list(NULL,c("effect","variance")))
  expect_no_error(meta(y=effect.data[,1],v=effect.data[,2],run=FALSE))
})

test_that("covariance sampling weights preserve unknown cross-covariances", {
  nm <- c("x", "y", "z")
  S <- matrix(c(1,NA,.2, NA,1,.3, .2,.3,1), 3,
              dimnames=list(nm,nm))
  V <- asyCov(S, 100, cor.analysis=FALSE)
  expect_true(is.na(V["x_x","y_y"]))
  expect_true(is.na(V["z_x","z_y"]))
  expect_equal(unname(V["z_x","z_x"]), (1 + .2^2)/100)
  expect_equal(unname(V["x_x","z_z"]), 2*.2^2/100)
  expect_true(all(is.na(V["y_x",])))
  expect_equal(V, t(V))

  C <- asyCov(S, 100, cor.analysis=TRUE)
  expect_true(is.na(C["z_x","z_y"]))
  expect_equal(unname(C["z_x","z_x"]), (1-.2^2)^2/100)

  ## A known zero determines a product even when its other factor is unknown.
  S[1,3] <- S[3,1] <- 0
  V <- asyCov(S, 100, cor.analysis=FALSE)
  expect_equal(unname(V["x_x","z_y"]), 0)
  expect_true(is.na(V["z_x","z_y"]))
  expect_true(all(is.na(V["y_x",])))
})

test_that("covariance sampling weights agree with normal-theory moments", {
  nm <- c("x", "y", "z")
  S <- matrix(c(1,.4,.2,.4,2,.3,.2,.3,3), 3, dimnames=list(nm,nm))
  V <- asyCov(S, 250, cor.analysis=FALSE)
  expect_equal(unname(V["x_x","y_y"]), 2*.4^2/250)
  expect_equal(unname(V["z_x","z_y"]), (.4*3+.2*.3)/250)
  expect_equal(unname(V["z_y","z_y"]), (2*3+.3^2)/250)
  if (requireNamespace("matrixcalc", quietly=TRUE)) {
    D <- matrixcalc::D.matrix(3)
    G <- MASS::ginv(D)
    reference <- 2 * G %*% (S %x% S) %*% t(G) / 250
    expect_equal(unname(V), unname(reference), tolerance=1e-12)
  }
  expect_equal(unname(asyCov(matrix(2,1,1),100,cor.analysis=FALSE)),
               matrix(.08,1,1))
})

test_that("pooled covariance weights use observed study averages", {
  nm <- c("x", "y", "z")
  S <- matrix(c(1,NA,.2, NA,1,.3, .2,.3,1),3,dimnames=list(nm,nm))
  T <- S; T[1,2] <- T[2,1] <- .8
  for (mode in c("weighted", "unweighted")) {
    V <- asyCov(list(S,T), c(100,200), cor.analysis=FALSE,
                acov=mode, as.matrix=FALSE)
    expect_equal(unname(V[[1]]["x_x","y_y"]),2*.8^2/100)
    expect_equal(unname(V[[2]]["x_x","y_y"]),2*.8^2/200)
    missing <- asyCov(list(S,S),c(100,200),cor.analysis=FALSE,
                      acov=mode,as.matrix=FALSE)
    expect_true(is.na(missing[[1]]["x_x","y_y"]))
  }
  V <- asyCov(list(S,T),c(100,200),cor.analysis=FALSE,
              acov="individual",as.matrix=FALSE)
  expect_true(is.na(V[[1]]["x_x","y_y"]))
  expect_equal(unname(V[[2]]["x_x","y_y"]),2*.8^2/200)
  packed <- asyCov(list(S,T),c(100,200),cor.analysis=FALSE,
                   acov="individual",as.matrix=TRUE)
  expect_true(is.na(packed[1,"C(y_y x_x)"]))
})

test_that("TSSEM requires an explicit working assumption for unknown sampling covariances", {
  nm <- c("x", "y", "z")
  S <- matrix(c(1,NA,.2, NA,1,.3, .2,.3,1),3,dimnames=list(nm,nm))
  T <- S; T[1,2] <- T[2,1] <- .5
  studies <- list(S,S,T,T)
  expect_error(tssem1(studies, rep(100,4), cor.analysis=FALSE,
                      acov="individual", RE.type="Zero", run=FALSE),
               "Missing sampling covariance")
  expect_no_error(tssem1(studies, rep(100,4), cor.analysis=FALSE,
                         acov="weighted", RE.type="Zero", run=FALSE))
})

test_that("asyCov rejects demonstrably non-positive-definite working matrices", {
  nm <- c("x", "y", "z")
  S <- diag(3); dimnames(S) <- list(nm,nm)
  S[1,2] <- S[2,1] <- 1.1
  for (cor in c(TRUE,FALSE)) {
    expect_error(asyCov(S,100,cor.analysis=cor), "not positive definite")
    expect_error(asyCov(list(S,S),c(100,100),cor.analysis=cor,acov="individual"),
                 "study 1.*not positive definite")
  }
  ## Each study reports a valid 2-by-2 submatrix, but their pairwise
  ## averages cannot form a positive-definite three-variable matrix.
  A <- diag(3); A[lower.tri(A)] <- A[upper.tri(A)] <- NA
  dimnames(A) <- list(nm,nm)
  B <- C <- A
  A[1,2] <- A[2,1] <- .9
  B[1,3] <- B[3,1] <- .9
  C[2,3] <- C[3,2] <- -.9
  for (cor in c(TRUE,FALSE)) for (mode in c("weighted","unweighted")) {
    expect_error(asyCov(list(A,B,C),rep(100,3),cor.analysis=cor,acov=mode),
                 paste0(mode," average.*not positive definite.*acov="))
    expect_error(asyCovOld(list(A,B,C),rep(100,3),cor.analysis=cor,acov=mode),
                 paste0(mode," average.*not positive definite.*acov="))
  }
  ## A missing pair remains unknown; it is not treated as zero for the check.
  A[1,2] <- A[2,1] <- NA
  expect_no_error(asyCov(A,100,cor.analysis=FALSE))
  A[1,2] <- A[2,1] <- 1.1
  expect_error(asyCov(A,100,cor.analysis=FALSE),"not positive definite")
  partial <- matrix(NA_real_,4,4); diag(partial) <- 1
  partial[1:3,1:3] <- matrix(c(1,.9,.9,.9,1,-.9,.9,-.9,1),3)
  expect_error(asyCov(partial,100,cor.analysis=FALSE),"not positive definite")
})

test_that("asyCovOld diagnoses missing observed pairs and retains missing variables", {
  nm <- c("x", "y", "z")
  S <- diag(3); dimnames(S) <- list(nm,nm)
  S[1,2] <- S[2,1] <- NA
  for (cor in c(TRUE,FALSE)) {
    expect_error(asyCovOld(S,100,cor.analysis=cor,run=FALSE),
                 "requires all covariances among retained variables")
    expect_error(asyCovOld(list(S,S),c(100,100),cor.analysis=cor,
                           acov="individual",as.matrix=FALSE),
                 "requires all covariances among retained variables")
  }
  S <- diag(3); dimnames(S) <- list(nm,nm)
  S[3,] <- S[,3] <- NA
  for (cor in c(TRUE,FALSE)) {
    expect_no_error(asyCovOld(S,100,cor.analysis=cor,dropNA=FALSE))
    expect_no_error(asyCovOld(S,100,cor.analysis=cor,dropNA=TRUE))
    expect_no_error(asyCovOld(list(S,S),c(100,120),cor.analysis=cor,
                            acov="weighted",as.matrix=TRUE))
  }
})

test_that("legacy TSSEM does not replace missing observed pairs with zero", {
  nm <- c("x", "y", "z")
  S <- matrix(c(1,NA,.2,NA,1,.3,.2,.3,1),3,dimnames=list(nm,nm))
  T <- S; T[1,2] <- T[2,1] <- .5
  for (cor in c(TRUE,FALSE)) {
    expect_error(tssem1(list(S,S,T,T),rep(100,4),cor.analysis=cor,
                        acov="individual",asyCovOld=TRUE,RE.type="Zero",run=FALSE),
                 "requires all covariances among retained variables")
    expect_no_error(fit <- tssem1(list(S,S,T,T),rep(100,4),cor.analysis=cor,
                           acov="weighted",asyCovOld=TRUE,RE.type="Zero",run=FALSE))
    if (!cor) {
      ## Cov(s_xx, s_yy) uses the reported .5, not the old zero-filled .25.
      expect_gt(fit@data@observed$v4_1[1], .004)
    }
    U <- T; U[3,] <- U[,3] <- NA
    expect_no_error(tssem1(list(U,U,T,T),rep(100,4),cor.analysis=cor,
                           acov="individual",asyCovOld=TRUE,RE.type="Zero",run=FALSE))
  }
})

test_that("TSSEM never silently imputes or drops a reported covariance", {
  nm <- c("x", "y")
  S <- matrix(c(NA,.4,.4,1),2,dimnames=list(nm,nm))
  T <- matrix(c(2,.4,.4,1),2,dimnames=list(nm,nm))
  for (legacy in c(FALSE,TRUE)) {
    for (mode in c("individual", "weighted"))
      expect_error(tssem1(list(S,S,T,T),rep(100,4),cor.analysis=TRUE,
                          acov=mode,asyCovOld=legacy,RE.type="Zero",run=FALSE),
                   "Study 1.*variable 'x'.*variance is missing.*calculate the correlation")
    message <- if (legacy) "Study 1.*variable 'x'.*variance is missing" else
                           "Study 1.*y_x.*sampling variance is unavailable"
    expect_error(tssem1(rep(list(S),4),rep(100,4),cor.analysis=FALSE,
                        acov="individual",asyCovOld=legacy,RE.type="Zero",run=FALSE),message)
    expect_error(tssem1(rep(list(S),4),rep(100,4),cor.analysis=FALSE,
                        acov="weighted",asyCovOld=legacy,RE.type="Zero",run=FALSE),
                 "Study 1.*y_x.*sampling variance is unavailable")
    ## A complete study provides the pooled working variance for x.
    model <- tssem1(list(S,S,T,T),rep(100,4),cor.analysis=FALSE,
                    acov="weighted",asyCovOld=legacy,RE.type="Zero",run=FALSE)
    expect_equal(model@data@observed$y2, rep(.4,4))
    expect_true(all(is.finite(model@data@observed$v2_2)))
  }
})

test_that("all known maximal submatrices are checked independent of order", {
  S <- matrix(NA_real_,5,5); diag(S) <- 1
  nm <- paste0("x",1:5); dimnames(S) <- list(nm,nm)
  S[3,4] <- S[4,3] <- .9; S[3,5] <- S[5,3] <- .9
  S[4,5] <- S[5,4] <- -.9
  S[1,3] <- S[3,1] <- 0; S[1,4] <- S[4,1] <- 0
  S[2,3] <- S[3,2] <- 0; S[2,5] <- S[5,2] <- 0
  expect_lt(min(eigen(S[3:5,3:5],symmetric=TRUE,only.values=TRUE)$values),0)
  for (order in list(1:5, c(3,4,5,1,2), 5:1, c(2,4,1,5,3))) {
    for (cor in c(FALSE,TRUE))
      expect_error(asyCov(S[order,order],100,cor.analysis=cor),"not positive definite")
  }
  S[4,5] <- S[5,4] <- .9
  for (order in list(1:5,5:1,c(2,4,1,5,3)))
    expect_no_error(asyCov(S[order,order],100,cor.analysis=FALSE))
  expect_warning(.check.asycov.matrix(S,work.limit=1),"validation.*incomplete.*search limit")
})

test_that("observed-submatrix search agrees with exhaustive small cases", {
  set.seed(1126)
  for (trial in seq_len(30)) {
    p <- 6L
    S <- diag(p)
    for (j in seq_len(p-1L)) for (i in (j+1L):p) {
      value <- if (runif(1) < .35) NA_real_ else runif(1,-.85,.85)
      S[i,j] <- S[j,i] <- value
    }
    invalid <- FALSE
    for (mask in 1:(2^p-1L)) {
      selected <- which(as.logical(intToBits(mask)[seq_len(p)]))
      sub <- S[selected,selected,drop=FALSE]
      if (!anyNA(sub) && inherits(tryCatch(chol(sub),error=identity),"error")) {
        invalid <- TRUE
        break
      }
    }
    result <- tryCatch(.check.asycov.matrix(S),error=identity)
    expect_identical(inherits(result,"error"),invalid,info=paste("trial",trial))
  }
})

test_that("legacy pooled TSSEM agrees with independent generalized least squares", {
  nm <- c("x","y")
  mk <- function(xx,xy,yy) matrix(c(xx,xy,xy,yy),2,dimnames=list(nm,nm))
  studies <- list(mk(NA,.3,1),mk(NA,.5,1),mk(1.8,.42,1.1),mk(2.2,.38,.9))
  n <- rep(100,4)
  set.seed(42)
  fit <- tssem1(studies,n,cor.analysis=FALSE,acov="weighted",
                asyCovOld=TRUE,RE.type="Zero")
  set.seed(42)
  V <- asyCovOld(studies,n,cor.analysis=FALSE,acov="weighted",as.matrix=FALSE)
  y <- t(vapply(studies,function(x) as.numeric(OpenMx::vech(x)),numeric(3)))
  precision <- matrix(0,3,3); score <- numeric(3)
  for (i in seq_along(n)) {
    idx <- which(!is.na(y[i,]))
    W <- solve(V[[i]][idx,idx,drop=FALSE])
    precision[idx,idx] <- precision[idx,idx] + W
    score[idx] <- score[idx] + W %*% y[i,idx]
  }
  reference <- solve(precision,score)
  objective <- 0
  for (i in seq_along(n)) {
    idx <- which(!is.na(y[i,]))
    Vi <- V[[i]][idx,idx,drop=FALSE]
    delta <- y[i,idx] - reference[idx]
    objective <- objective + length(idx)*log(2*pi) +
      as.numeric(determinant(Vi,logarithm=TRUE)$modulus) +
      drop(t(delta) %*% solve(Vi,delta))
  }
  expect_true(fit$mx.fit@output$status$code %in% c(0,1))
  expect_equal(unname(coef(fit)),as.numeric(reference),tolerance=1e-6)
  expect_equal(unname(vcov(fit)),unname(solve(precision)),tolerance=1e-6)
  expect_gt(min(eigen(vcov(fit),symmetric=TRUE,only.values=TRUE)$values),0)
  expect_equal(fit$mx.fit@output$fit,objective,tolerance=1e-6)
})
