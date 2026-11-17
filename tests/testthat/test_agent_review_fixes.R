test_that("sampling covariance inputs have symmetric missingness and finite values", {
  S <- matrix(c(1,.3,NA,1),2)
  for (cor in c(FALSE,TRUE))
    expect_error(asyCov(S,100,cor.analysis=cor),"asymmetric missing-value")
  S <- matrix(c(NA,Inf,Inf,1),2)
  expect_error(asyCov(S,100,cor.analysis=FALSE),"non-finite observed")
  S <- matrix(c(1,NA,NA,1),2)
  expect_true(is.na(asyCov(S,100,cor.analysis=FALSE)[1,3]))
})

test_that("explicit A bounds override only their specified labels", {
  nm <- c("x","y","z")
  C <- matrix(c(1,.3,.2,.3,1,.4,.2,.4,1),3,dimnames=list(nm,nm))
  data <- Cor2DataFrame(rep(list(C),4),n=rep(100,4))
  ram <- suppressWarnings(lavaan2RAM("y ~ a*x + b*z\na > 0.1\nb > 0.2\na < 0.8\nb < 0.7",obs.variables=nm))
  L <- U <- matrix(NA_real_,3,3,dimnames=dimnames(ram$A))
  L["y","x"] <- 0; U["y","x"] <- .9
  fit <- suppressWarnings(osmasem(RAM=ram,data=data,A.lbound=L,A.ubound=U,RE.type="Zero",run=FALSE))
  expect_equal(unname(omxGetParameters(fit$mx.fit,fetch="lbound")[c("a","b")]), c(0,.2))
  expect_equal(unname(omxGetParameters(fit$mx.fit,fetch="ubound")[c("a","b")]), c(.9,.7))
})

test_that("multigroup names and equality checks include all groups and levels", {
  set.seed(29); d <- data.frame(x=rnorm(60),y=rnorm(60),g=rep(c("foo","B"),each=30))
  ram <- lavaan2RAM("y ~ c(b,b)*x\nfoo := b*b",ngroups=2)
  fit <- sem(RAM=ram,data=d,group="g",run=FALSE)
  expect_false(fit$group.map[["foo"]] == "foo")
  ram <- lavaan2RAM("y ~ c(a,a)*x\na == exp(newparam)",ngroups=2)
  expect_error(sem(RAM=ram,data=d,group="g",run=FALSE),"introduced only.*newparam")
  ram <- lavaan2RAM("level: 1\ny ~ a*x\nlevel: 2\ny ~~ b*y\na == exp(newparam)")
  d$id <- rep(1:20,each=3)
  expect_error(sem(RAM=ram,data=d,cluster="id",run=FALSE),"introduced only.*newparam")
})

test_that("WLS aligns reversed pair labels and handles covariance subsets", {
  nm <- c("x","y","z"); C <- matrix(c(1,.1,.3,.1,1,.6,.3,.6,1),3,dimnames=list(nm,nm))
  ram <- lavaan2RAM("y ~ b*x\nz ~ b*x\nx ~~ 1*x",obs.variables=nm)
  V <- diag(c(.01,.04,.09));dimnames(V) <- list(c("y_x","z_x","z_y"),c("y_x","z_x","z_y"))
  alt <- V[c(3,1,2),c(3,1,2)]; dimnames(alt) <- list(c("y_z","x_y","z_x"),c("y_z","x_y","z_x"))
  a <- wls(C,V,500,RAM=ram); b <- wls(C,alt,500,RAM=ram)
  expect_true(a$mx.fit@output$status$code %in% c(0,1))
  expect_true(b$mx.fit@output$status$code %in% c(0,1))
  expect_equal(coef(a),coef(b),tolerance=1e-7)
  expect_equal(a$mx.fit@output$fit,b$mx.fit@output$fit,tolerance=1e-7)
  expect_no_error(wls(C,diag(6),100,cor.analysis=FALSE,subset.variables="x",
                      Amatrix=matrix(0,1,1),Smatrix=matrix("0.5*s",1,1),run=FALSE))
  expect_error(wls(C,diag(3),100,cor.analysis=TRUE,subset.variables="x",run=FALSE),"at least two")
  expect_error(wls(C,diag(3),100,RAM=lavaan2RAM("y~x",ngroups=2),run=FALSE),"single-group")
})

test_that("WLS likelihood intervals omit fixed S labels", {
  C <- matrix(c(1,.3,.2,.3,1,.4,.2,.4,1),3)
  S <- matrix("0",3,3);diag(S)<-1
  S[1,2]<-S[2,1]<-"0.1@fixed_cov";S[1,3]<-S[3,1]<-"0.2*b"
  m<-wls(C,diag(3)*.01,100,Amatrix=matrix(0,3,3),Smatrix=S,intervals.type="LB",run=FALSE)
  expect_false("fixed_cov" %in% names(m@intervals))
  expect_true("b" %in% names(m@intervals))
})

test_that("clustered FEM forwards run and SEC validates included data", {
  C <- matrix(c(1,.4,.4,1),2)
  fits <- tssem1FEM(rep(list(C),4),rep(100,4),cluster=c(1,1,2,2),run=FALSE)
  expect_true(all(vapply(fits,inherits,logical(1),"MxModel")))
  expect_equal(.fitIndices(1,2,2,3,100)$CFI,1)
  expect_error(SEC(diag(2),V=c(1,Inf)),"positive")
  X <- matrix(c(1,1,0,1,0,0),3,2)
  expect_no_error(SEC(X,method="weighted",q=c(1,1,-1)))
  expect_error(SEC(diag(2),method="size",n=c(10,Inf)),"finite")
})

test_that("multi-column moderators get usable names", {
  nm <- c("x","y","z"); C <- diag(3);dimnames(C)<-list(nm,nm)
  d <- list(data=rep(list(C),4),n=rep(100,4),mods=data.frame(a=1:4,b=5:8))
  obj <- Cor2DataFrame(d)
  expect_true(all(c("mods_a","mods_b") %in% colnames(obj$data)))
  expect_no_error(osmasem(RAM=lavaan2RAM("y~x\nx~~1*x\nz~~1*z",obs.variables=nm),data=obj,run=FALSE))
})

test_that("non-finite curvature and singular constrained information are diagnosed", {
  active <- list(jacobian=matrix(1,1,1),estimate=c(a=0),values=function(z) z["a"])
  expect_error(.constraint.curvature(active,matrix(1,1,1),c(a=Inf),"a"),"non-finite")
  status <- list(hessian=matrix(0,2,2),active=NULL)
  expect_error(.constrained.vcov(NULL,p=c("a","b"),bounded=c(FALSE,FALSE),status=status),"singular.*identified")
})

test_that("osmasem2 reports both failures when its retry errors", {
  C <- matrix(c(1,.2,.2,1),2,dimnames=list(c("x","y"),c("x","y")))
  data <- Cor2DataFrame(rep(list(C),4),n=rep(100,4))
  ram <- lavaan2RAM("y~x\nx~~1*x",obs.variables=c("x","y"))
  local_mocked_bindings(mxRun=function(...) stop("initial failure"),
                        mxTryHard=function(...) stop("retry failure"),.package="metaSEM")
  expect_error(osmasem2(RAM=ram,data=data,RE.type.Sigma="Zero"),
               "initial failure.*mxTryHard.*retry failure")
})

test_that("row gradients fall back for model names containing spaces", {
  set.seed(92); d <- data.frame(x=rnorm(30),y=rnorm(30))
  fit <- sem(model.name="my model",RAM=lavaan2RAM("y~x"),data=d)
  expect_true(fit$mx.fit@output$status$code %in% c(0,1))
  expect_error(.row.m2ll.differences(fit$mx.fit,coef(fit$mx.fit)))
  invisible(capture.output(actual <- .row.m2ll.gradients(fit$mx.fit)))
  invisible(capture.output(reference <- OpenMx::imxRowGradients(fit$mx.fit,robustSE=TRUE)))
  expect_equal(unname(actual),unname(as.matrix(reference)),tolerance=1e-8)
})

test_that("shared A labels receive explicit bounds consistently", {
  nm <- c("x","y","z"); C <- diag(3);dimnames(C)<-list(nm,nm)
  data <- Cor2DataFrame(rep(list(C),4),n=rep(100,4))
  ram <- suppressWarnings(lavaan2RAM("y~a*x\nx~~a*z\na>0.2",obs.variables=nm))
  L <- matrix(NA_real_,3,3,dimnames=dimnames(ram$A)); L["y","x"] <- 0
  fit <- suppressWarnings(osmasem(RAM=ram,data=data,A.lbound=L,RE.type="Zero",run=FALSE))
  expect_equal(unname(omxGetParameters(fit$mx.fit,fetch="lbound")["a"]),0)
})

test_that("as.mxMatrix preserves variable names containing spaces", {
  x <- matrix(1:4,2,dimnames=list(c("x 1","y"),NULL))
  m <- as.mxMatrix(x,name="test")
  expect_equal(rownames(m$values),c("x 1","y"))
  expect_null(colnames(m$values))
})

test_that("unknown standard errors are not silently fixed", {
  expect_warning(selected <- .rerun.invalid.tau.se(c(NA_real_,NaN,.1)),"cause is unclear")
  expect_identical(unname(selected),c(FALSE,TRUE,FALSE))
})

test_that("multigroup equality checks recognize auxiliary model parameters", {
  set.seed(92); d <- data.frame(x=rnorm(40),y=rnorm(40),g=rep(c("A","B"),each=20))
  ram <- lavaan2RAM("y ~ c(a,a)*x\na == exp(newparam)",ngroups=2)
  extra <- OpenMx::mxMatrix("Full",1,1,free=TRUE,values=0,labels="newparam",name="additional")
  expect_no_error(sem(RAM=ram,data=d,group="g",mxModel.Args=list(extra),run=FALSE))
})

test_that("multigroup defined-parameter name collisions preserve estimates and fit", {
  set.seed(919); n <- 250; d <- data.frame(x=rnorm(n),g=rep(c("foo","B"),each=n/2))
  d$y <- .4*d$x+rnorm(n)
  syntax <- "y ~ c(b,b)*x\nfoo := b*b"
  fit <- sem(RAM=lavaan2RAM(syntax,ngroups=2),data=d,group="g")
  ref <- lavaan::sem(syntax,data=d,group="g",fixed.x=FALSE,meanstructure=TRUE)
  expect_true(fit$mx.fit@output$status$code %in% c(0,1))
  expect_true(lavaan::lavInspect(ref,"converged"))
  expect_equal(unname(coef(fit)["b"]),unname(lavaan::coef(ref)["b"]),tolerance=1e-5)
  expect_equal(fit$mx.fit@output$fit,-2*as.numeric(logLik(ref)),tolerance=1e-5)
})

test_that("one-variable covariance WLS fits the observed variance", {
  C <- matrix(c(1,.2,.2,2),2,dimnames=list(c("x","y"),c("x","y")))
  fit <- wls(C,diag(3)*.01,500,cor.analysis=FALSE,subset.variables="x",
             Amatrix=matrix(0,1,1),Smatrix=matrix("0.5*s",1,1))
  expect_true(fit$mx.fit@output$status$code %in% c(0,1))
  expect_equal(unname(coef(fit)["s"]),1,tolerance=1e-6)
  expect_equal(fit$mx.fit@output$fit,0,tolerance=1e-8)
})

test_that("inherited OSMASEM bounds agree with independent weighted optimization", {
  nm <- c("x","y","z")
  C <- matrix(c(1,.25,.1,.25,1,-.3,.1,-.3,1),3,dimnames=list(nm,nm))
  data <- Cor2DataFrame(rep(list(C),6),n=rep(100,6))
  ram <- lavaan2RAM("y~a*x+b*z\nx~~1*x\nz~~1*z\nx~~r*z\nb>0.2",obs.variables=nm)
  L <- matrix(NA_real_,3,3,dimnames=dimnames(ram$A)); L["y","x"] <- -1
  fit <- osmasem(RAM=ram,data=data,A.lbound=L,RE.type="Zero")
  weight <- solve(asyCov(C,100,cor.analysis=TRUE))
  objective <- function(q) {
    a<-q[1];b<-q[2];r<-q[3]
    residual <- c(.25,.1,-.3)-c(a+b*r,r,a*r+b)
    6*drop(t(residual)%*%weight%*%residual)
  }
  reference <- optim(c(.2,.2,.1),objective,method="L-BFGS-B",
                     lower=c(-1,.2,-.99),upper=c(1,1,.99))
  expect_true(fit$mx.fit@output$status$code %in% c(0,1))
  expect_identical(reference$convergence,0L)
  expect_equal(unname(coef(fit)[c("a","b","r")]),reference$par,tolerance=1e-5)
  expect_equal(objective(coef(fit)[c("a","b","r")]),reference$value,tolerance=1e-5)
})

test_that("WLS preserves observed names containing spaces", {
  nm <- c("x 1","y"); C <- matrix(c(1,.3,.3,1),2,dimnames=list(nm,nm))
  A <- matrix("0",2,2,dimnames=list(nm,nm)); A["y","x 1"] <- "0.2*b"
  S <- F <- diag(2); dimnames(S) <- dimnames(F) <- list(nm,nm)
  fit <- suppressWarnings(wls(C,matrix(.01),100,Amatrix=A,Smatrix=S,Fmatrix=F))
  expect_true(fit$mx.fit@output$status$code %in% c(0,1))
  expect_equal(unname(coef(fit)["b"]),.3,tolerance=1e-6)
})

test_that("combined group and level models validate equality parameters", {
  syntax <- paste("group: A", "level: 1", "y ~ a*x", "level: 2", "y ~~ b*y",
                  "group: B", "level: 1", "y ~ a*x", "level: 2", "y ~~ c*y",
                  "a == exp(newparam)",sep="\n")
  ram <- lavaan2RAM(syntax,ngroups=2)
  set.seed(29); d <- data.frame(x=rnorm(60),y=rnorm(60),
                               g=rep(c("A","B"),each=30),id=rep(1:20,each=3))
  expect_error(suppressWarnings(sem(RAM=ram,data=d,group="g",cluster="id",run=FALSE)),
               "introduced only.*newparam")
})

test_that("ambiguous reverse pair aliases are rejected", {
  nm <- c("a_b","a","b_c","c"); C <- diag(4);dimnames(C)<-list(nm,nm)
  labs <- .genCorNames(nm,cor.analysis=TRUE)$ylabels
  labs[1] <- "a_b_c"
  V <- diag(6);dimnames(V)<-list(labs,labs)
  expect_error(.align.wls.acov(C,V,TRUE),"Ambiguous")
})

test_that("generated moderator name collisions are diagnosed", {
  nm <- c("x","y"); C <- diag(2);dimnames(C)<-list(nm,nm)
  data <- list(data=rep(list(C),4),n=rep(100,4),mods=data.frame(a=1:4,b=5:8),mods_a=1:4)
  expect_error(Cor2DataFrame(data),"column names collide")
})

test_that("TSSEM forwards optimizer warning and silence options", {
  C <- matrix(c(1,.4,.4,1),2)
  local({
    seen <- NULL
    local_mocked_bindings(meta=function(...){seen <<- list(...); OpenMx::mxModel("placeholder")},.package="metaSEM")
    tssem1REM(rep(list(C),4),rep(100,4),RE.type="Zero",suppressWarnings=FALSE,run=FALSE)
    expect_false(seen$suppressWarnings)
  })
  local({
    seen <- logical()
    local_mocked_bindings(mxRun=function(model,silent,...) {
      seen <<- c(seen,setNames(silent,model@name))
      OpenMx::mxRun(model,silent=TRUE,...)
    },.package="metaSEM")
    invisible(capture.output(tssem1FEM(rep(list(C),4),rep(100,4),
                                      cluster=c(1,1,2,2),silent=FALSE,model.name="myfit")))
    expect_identical(unname(seen[names(seen)=="myfit"]),c(FALSE,FALSE))
  })
})

test_that("rerun detects transformations in multiline Tau2 expressions", {
  formula <- str2lang(paste0("exp(x) + ",paste(rep("x",80),collapse=" + ")))
  expect_gt(length(deparse(formula)),1)
  expect_equal(.rerun.tau.zero(formula),-Inf)
  expect_equal(.rerun.tau.zero(quote(x^2)),0)
})

test_that("non-finite constraint values give an actionable error", {
  m <- OpenMx::mxModel("test",OpenMx::mxMatrix("Full",1,1,free=TRUE,values=1,labels="p",name="M"))
  local_mocked_bindings(.constraint.list=function(...) list(dummy=TRUE),
                        .constraint.values=function(...) NaN,.package="metaSEM")
  expect_error(.active.constraints(m,"p"),"constraint values are not finite")
})

test_that("a numeric bound on any occurrence applies to a shared A label", {
  nm <- c("x","y","z"); C <- diag(3);dimnames(C)<-list(nm,nm)
  data <- Cor2DataFrame(rep(list(C),4),n=rep(100,4))
  ram <- suppressWarnings(lavaan2RAM("y~a*x\nz~a*x\na>0.2",obs.variables=nm))
  L <- matrix(NA_real_,3,3,dimnames=dimnames(ram$A));L["z","x"] <- 0
  fit <- suppressWarnings(osmasem(RAM=ram,data=data,A.lbound=L,RE.type="Zero",run=FALSE))
  expect_equal(unname(omxGetParameters(fit$mx.fit,fetch="lbound")["a"]),0)
  L["y","x"] <- .1
  expect_error(suppressWarnings(osmasem(RAM=ram,data=data,A.lbound=L,RE.type="Zero",run=FALSE)),
               "Conflicting A bounds.*a")
})

test_that("NA scalar A bounds retain all RAM bounds", {
  nm <- c("x","y","z"); C <- diag(3);dimnames(C)<-list(nm,nm)
  data <- Cor2DataFrame(rep(list(C),4),n=rep(100,4))
  ram <- suppressWarnings(lavaan2RAM("y~a*x+b*z\na>0.1\nb>0.2",obs.variables=nm))
  fit <- suppressWarnings(osmasem(RAM=ram,data=data,A.lbound=NA,RE.type="Zero",run=FALSE))
  expect_equal(unname(omxGetParameters(fit$mx.fit,fetch="lbound")[c("a","b")]),c(.1,.2))
})

test_that("SEC ignores invalid values only for wholly excluded studies", {
  X <- rbind(diag(2),c(0,0))
  expect_no_error(SEC(X,method="weighted",q=c(1,1,NA_real_)))
  W <- rbind(matrix(1,2,2),c(NA,Inf))
  expect_no_error(SEC(X,method="weighted",q=c(1,1,1),weights=W))
  expect_error(SEC(X,method="weighted",q=c(1,Inf,1)),"finite")
})

test_that("equality checks allow qualified OpenMx group matrix references", {
  set.seed(19); d <- data.frame(x=rnorm(100),y=rnorm(100),g=rep(c("A","B"),each=50))
  ram <- lavaan2RAM("y ~ c(a,b)*x",ngroups=2)
  ram[[1]]$mxalgebras <- list(co=OpenMx::mxConstraint(a == B.Amatrix[1,2],name="co"))
  fit <- sem(RAM=ram,data=d,group="g",run=FALSE)
  expect_no_error(OpenMx::mxRun(fit$mx.fit,useOptimizer=FALSE,silent=TRUE))
})

test_that("Becker09 reported correlations retain their unit diagonals in TSSEM", {
  data(Becker09,package="metaSEM")
  expect_error(tssem1(Becker09$data,Becker09$n,method="REM",run=FALSE),
               "Study 5.*variance is missing")
  corrected <- lapply(Becker09$data,function(S) {
    observed <- vapply(seq_len(nrow(S)),function(i) any(!is.na(S[i,-i])),logical(1))
    d <- diag(S); d[is.na(d) & observed] <- 1; diag(S) <- d; S
  })
  m <- tssem1(corrected,Becker09$n,method="REM",run=FALSE)
  raw <- suppressWarnings(list2matrix(lapply(Becker09$data,cov2cor),diag=FALSE))
  fixed <- suppressWarnings(list2matrix(lapply(corrected,cov2cor),diag=FALSE))
  gained <- which(is.na(raw[5,]) & !is.na(fixed[5,]))
  expect_length(gained,2)
  expect_equal(unname(as.numeric(m@data@observed[5,gained])),unname(fixed[5,gained]))
})

test_that("Roorda11 correlation example retains reported entries with unit diagonals", {
  data(Roorda11,package="metaSEM")
  expect_error(tssem1(Roorda11$data,Roorda11$n,method="REM",run=FALSE),
               "variance is missing")
  corrected <- lapply(Roorda11$data,function(S) {
    observed <- vapply(seq_len(nrow(S)),function(i) any(!is.na(S[i,-i])),logical(1))
    d <- diag(S); d[is.na(d) & observed] <- 1; diag(S) <- d; S
  })
  m <- tssem1(corrected,Roorda11$n,method="REM",run=FALSE)
  raw <- suppressWarnings(list2matrix(lapply(Roorda11$data,cov2cor),diag=FALSE))
  fixed <- suppressWarnings(list2matrix(lapply(corrected,cov2cor),diag=FALSE))
  gained <- which(is.na(raw) & !is.na(fixed),arr.ind=TRUE)
  expect_gt(nrow(gained),0)
  expect_true(all(!is.na(as.matrix(m@data@observed[,seq_len(ncol(fixed))])[gained])))
})
