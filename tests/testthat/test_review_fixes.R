test_that("SMD functions retain zero means and single endpoints", {
  for (hom in c("variance", "none")) {
    a <- smdMTS(c(0, 1, 2), rep(1, 3), rep(30, 3),
                homogeneity=hom, bias.adjust=FALSE)
    b <- smdMTS(c(10, 11, 12), rep(1, 3), rep(30, 3),
                homogeneity=hom, bias.adjust=FALSE)
    expect_equal(unname(a$y), c(1, 2), tolerance=1e-5)
    expect_equal(a$V, b$V, tolerance=1e-5)
  }
  for (hom in c("covariance", "correlation", "none")) {
    a <- smdMES(0, 1, matrix(1), matrix(1), 30, 30,
                homogeneity=hom, bias.adjust=FALSE)
    expect_equal(unname(a$y), 1, tolerance=1e-5)
    expect_equal(dim(a$V), c(1L, 1L))
    expect_true(is.finite(a$V[1, 1]))
    b <- smdMES(c(0, NA), c(1, 2), diag(2), diag(2), 30, 30,
                homogeneity=hom, bias.adjust=FALSE)
    expect_equal(b$y[1], a$y[1])
    expect_true(is.na(b$y[2]))
    V <- diag(c(1, NA_real_))
    b <- smdMES(c(0, 1), c(1, 2), V, diag(2), 30, 30,
                homogeneity=hom, bias.adjust=FALSE)
    expect_equal(b$V[1, 1], a$V[1, 1], tolerance=1e-5)
    expect_true(is.na(b$y[2]))
  }
  expect_error(smdMES(NA, 1, matrix(1), matrix(1), 30, 30), "At least one endpoint")
})

test_that("FIML models accept intercept and random-effect constraints", {
  d <- data.frame(y=c(.1,.3,.2,.4,.1,.5), v=.05,
                  x=c(1,2,1,2,1,2), id=rep(1:3,each=2))
  m <- metaFIML(y, v, x=x, data=d, intercept.constraints=matrix("0.2*mu"), run=FALSE)
  expect_equal(m$Mmatrix$values[1,1], .2)
  expect_identical(m$Mmatrix$labels[1,1], "mu")
  m <- metaFIML(y, v, x=x, data=d, intercept.constraints=matrix(0), run=FALSE)
  expect_false(m$Mmatrix$free[1,1])
  m <- metaFIML(cbind(y,y+.1), cbind(v,0,v), x=x, data=d,
                intercept.constraints=c("0.2*mu1","0.3*mu2"), run=FALSE)
  expect_equal(unname(m$Mmatrix$values[1,1:2]), c(.2,.3))
  expect_equal(unname(m$Mmatrix$labels[1,1:2]), c("mu1","mu2"))
  for (lev in c(2, 3)) {
    args <- list(y=d$y, v=d$v, cluster=d$id, x2=d$x, run=FALSE)
    args[[paste0("RE",lev,".constraints")]] <- "0.2*tau"
    m <- do.call(meta3LFIML, args)
    expect_equal(m$Tau$values[lev-1,lev-1], .2)
    expect_identical(m$Tau$labels[lev-1,lev-1], "tau")
    args[[paste0("RE",lev,".constraints")]] <- 0
    m <- do.call(meta3LFIML, args)
    expect_false(m$Tau$free[lev-1,lev-1])
  }
})

review_osmasem_data <- function(means=FALSE) {
  C <- matrix(c(1,.2,.3,.2,1,.4,.3,.4,1),3,
              dimnames=list(c("x","y","z"), c("x","y","z")))
  if (!means) return(Cor2DataFrame(rep(list(C),6), n=rep(100,6), cor.analysis=FALSE))
  M <- matrix(rep(c(.2,.3,.4),6), ncol=3, byrow=TRUE,
              dimnames=list(NULL,colnames(C)))
  Cor2DataFrame(rep(list(C),6), n=rep(100,6), cor.analysis=FALSE, Means=M)
}

test_that("osmasem2 resolves chained replacements and remaining algebras", {
  d <- review_osmasem_data()
  ram <- lavaan2RAM("y ~ a*x\nz ~ b*y\nb == 2*a\na == exp(k)\nind := a*b\nk > -4")
  expect_equal(attr(ram,"lbound")$k,-4)
  expect_false(any(vapply(ram$mxalgebras,function(z)
    inherits(z,"MxConstraint") && identical(z$formula,str2lang("k > -4")),logical(1))))
  a <- osmasem2(RAM=ram,data=d,cor.analysis=FALSE,replace.constraints=TRUE,
                startvalues=list(k=-2),RE.type.Sigma="Zero")
  reference <- lavaan2RAM("y ~ a*x\nz ~ b*y\nb == 2*a\nind := a*b\na > 0.0183156388887342")
  b <- osmasem2(RAM=reference,data=d,cor.analysis=FALSE,replace.constraints=FALSE,
                RE.type.Sigma="Zero")
  expect_equal(OpenMx::mxEvalByName("expSigma",a$mx.fit),
               OpenMx::mxEvalByName("expSigma",b$mx.fit), tolerance=1e-4)
  expect_equal(a$mx.fit$output$fit,b$mx.fit$output$fit,tolerance=1e-4)
  k <- unname(coef(a)["k"])
  expect_equal(as.numeric(OpenMx::mxEvalByName("ind",a$mx.fit)), 2*exp(k)^2, tolerance=1e-6)
  expect_true(any(vapply(a$mx.fit@matrices, function(m) any(m$lbound == -4, na.rm=TRUE), logical(1))))
  if (requireNamespace("semPlot", quietly=TRUE)) {
    grDevices::pdf(tempfile(fileext=".pdf"))
    on.exit(grDevices::dev.off(), add=TRUE)
    expect_error(plot(a), NA)
  }
})

test_that("osmasem2 supports constant replacements and algebra means", {
  d <- review_osmasem_data(TRUE)
  ram <- lavaan2RAM("y ~ a*x\na == exp(-2)\nx ~~ vx*x\ny ~~ vy*y\nz ~~ vz*z\nvx == 1\nvy == 1\nvz == 1\nx ~ m*1\nm == exp(q)")
  a <- osmasem2(RAM=ram,data=d,cor.analysis=FALSE,mean.analysis=TRUE,
                replace.constraints=TRUE,run=FALSE)
  expect_s4_class(a$mx.fit$Amatrix,"MxMatrix")
  expect_s4_class(a$mx.fit$Smatrix,"MxMatrix")
  expect_equal(a$mx.fit$Amatrix$values["y","x"], exp(-2))
  expect_s4_class(a$mx.fit$Mmatrix,"MxAlgebra")
  expect_true("q" %in% names(OpenMx::omxGetParameters(a$mx.fit)))
  a <- osmasem2(RAM=ram,data=d,cor.analysis=FALSE,mean.analysis=TRUE,
                replace.constraints=TRUE, RE.type.Sigma="Zero", RE.type.Mu="Zero",
                startvalues=list(q=log(.2)))
  expect_equal(as.numeric(OpenMx::mxEvalByName("Mmatrix",a$mx.fit))[1],
               .3-exp(-2)*.2, tolerance=1e-4)
  if (requireNamespace("semPlot", quietly=TRUE)) {
    grDevices::pdf(tempfile(fileext=".pdf"))
    on.exit(grDevices::dev.off(), add=TRUE)
    expect_error(plot(a), NA)
  }
})

test_that("osmasem2 uses substituted diagonal expressions for correlation constraints", {
  C <- matrix(c(1,.2,.3,.2,1,.4,.3,.4,1),3,
              dimnames=list(c("x","y","z"), c("x","y","z")))
  d <- Cor2DataFrame(rep(list(C),6), n=rep(100,6))
  ram <- lavaan2RAM("y ~ a*x\nx ~~ vx*x\ny ~~ vy*y\nz ~~ vz*z\nvx == 1\nvy == 1-a*a\nvz == 1")
  a <- suppressWarnings(osmasem2(RAM=ram,data=d,replace.constraints=TRUE,run=FALSE))
  expect_length(a$mx.fit@constraints,0L)
  a <- suppressWarnings(osmasem2(RAM=ram,data=d,replace.constraints=TRUE,
                                  RE.type.Sigma="Zero"))
  expect_equal(diag(OpenMx::mxEvalByName("expSigma",a$mx.fit)), rep(1,3), tolerance=1e-5)
  expect_equal(summary(a)$degreesOfFreedom, 17)
  explicit <- suppressWarnings(osmasem2(RAM=ram,data=d,replace.constraints=FALSE,
                                         RE.type.Sigma="Zero"))
  expect_equal(summary(explicit)$degreesOfFreedom,17)
  expect_equal(coef(a)["a"],coef(explicit)["a"],tolerance=1e-4)
  expect_equal(as.numeric(vcov(a)),as.numeric(vcov(explicit)["a","a"]),tolerance=1e-4)
  if (requireNamespace("semPlot", quietly=TRUE)) {
    grDevices::pdf(tempfile(fileext=".pdf"))
    on.exit(grDevices::dev.off(), add=TRUE)
    expect_error(plot(a), NA)
  }
})

test_that("SMD sample sizes and retained covariances are validated", {
  for (n in list(NA_real_,Inf,0,1,-1,numeric(0),c(20,30),"20")) {
    expect_error(smdMES(0,1,matrix(1),matrix(1),n,30),"finite sample size greater than 1")
    expect_error(smdMES(0,1,matrix(1),matrix(1),30,n),"finite sample size greater than 1")
  }
  V <- matrix(c(1,NA,NA,1),2)
  expect_error(smdMES(c(0,1),c(1,2),V,diag(2),30,30),"Covariances among retained")
})

test_that("introduced parameter bounds are extracted end to end", {
  ram <- lavaan2RAM("y ~ a*x\na == exp(k)\nk > -4\nk < 1")
  expect_equal(attr(ram,"lbound")$k,-4)
  expect_equal(attr(ram,"ubound")$k,1)
  d <- review_osmasem_data()
  expect_error(osmasem2(RAM=ram,data=d,subset.variables=c("x","y"),cor.analysis=FALSE,
                        run=FALSE),"require 'replace.constraints=TRUE'")
  C <- matrix(c(1,.2,.2,1),2,dimnames=list(c("x","y"),c("x","y")))
  expect_error(sem(RAM=ram,Cov=C,numObs=100,run=FALSE),"require 'replace.constraints=TRUE'")
  fit <- sem(RAM=ram,Cov=C,numObs=100,replace.constraints=TRUE)
  expect_equal(coef(fit)["k"],c(k=log(.2)),tolerance=1e-4)
  expect_false(any(vapply(fit$mx.fit@constraints,function(z)
    "k" %in% all.vars(z$formula),logical(1))))
  auxiliary <- lavaan2RAM("y ~ a*x\na == auxiliary")
  extras <- list(OpenMx::mxMatrix("Full",1,1,free=TRUE,values=-2,
                                  labels="k",name="extraParameter"),
                 OpenMx::mxAlgebraFromString("exp(k)",name="auxiliary"))
  expect_s4_class(sem(RAM=auxiliary,Cov=C,numObs=100,mxModel.Args=extras,
                       run=FALSE)$mx.fit,"MxModel")
})

test_that("correlation restrictions preserve reparameterizations and real restrictions", {
  C <- matrix(c(1,.2,.2,1),2,dimnames=list(c("x","y"),c("x","y")))
  d <- Cor2DataFrame(rep(list(C),2),n=c(100,100))
  ram <- lavaan2RAM("y ~ a*x\nx ~~ 1*x\ny ~~ vy*y\nvy == exp(k)")
  fit <- osmasem2(RAM=ram,data=d,replace.constraints=TRUE,
                  RE.type.Sigma="Zero",startvalues=list(k=log(.96)))
  expect_equal(unname(coef(fit)["a"]),.2,tolerance=1e-4)
  expect_equal(summary(fit)$degreesOfFreedom,1)
  expect_equal(diag(OpenMx::mxEvalByName("expSigma",fit$mx.fit)),c(1,1),tolerance=1e-5)
  ram <- lavaan2RAM("y ~ a*x\nx ~~ 1*x\ny ~~ vy*y\nvy == 0.5*a")
  fit <- expect_warning(osmasem2(RAM=ram,data=d,replace.constraints=TRUE,RE.type.Sigma="Zero",
                                 startvalues=list(a=.7)), "no basis")
  expect_equal(unname(coef(fit)["a"]),(-.5+sqrt(4.25))/2,tolerance=1e-4)
  expect_equal(summary(fit)$degreesOfFreedom,2)
  expect_equal(as.numeric(vcov(fit)),0,tolerance=1e-8)
})

test_that("unit variances hold at each moderator pattern", {
  C <- matrix(c(1,.2,.2,1),2,dimnames=list(c("x","y"),c("x","y")))
  d <- Cor2DataFrame(rep(list(C),2),n=c(100,100));d$data$mod <- c(10,20)
  ram <- lavaan2RAM("y ~ a*x\nx ~~ 1*x\ny ~~ v*y\na == 0.01*data.mod\nv == v0+v1*data.mod")
  fit <- expect_warning(osmasem2(RAM=ram,data=d,replace.constraints=TRUE,RE.type.Sigma="Zero",
                                 startvalues=list(v0=1.02,v1=-.003)), "no basis")
  expect_length(fit$mx.fit@constraints,2L)
  for (mod in c(10,20)) {
    m <- fit$mx.fit
    m@data@observed[1,"mod"] <- mod
    expect_equal(unname(diag(OpenMx::mxEvalByName("expSigma",m,compute=TRUE))),c(1,1),tolerance=1e-5)
  }
  reversed <- d; reversed$data <- d$data[2:1,,drop=FALSE]
  other <- expect_warning(osmasem2(RAM=ram,data=reversed,replace.constraints=TRUE,RE.type.Sigma="Zero",
                                   startvalues=list(v0=1.02,v1=-.003)), "no basis")
  expect_equal(coef(fit),coef(other),tolerance=1e-5)
  expect_equal(as.numeric(vcov(fit)),rep(0,4),tolerance=1e-8)
  incompatible <- lavaan2RAM("y ~ a*x\nx ~~ 1*x\ny ~~ v*y\na == 0.01*data.mod")
  expect_error(osmasem2(RAM=incompatible,data=d,replace.constraints=TRUE,run=FALSE),
               "incompatible across studies")
})

test_that("plot matrices use moderator means for direct and substituted definitions", {
  C <- matrix(c(1,.2,.2,1),2,dimnames=list(c("x","y"),c("x","y")))
  d <- Cor2DataFrame(rep(list(C),2),n=c(100,100),cor.analysis=FALSE)
  d$data$mod <- c(10,20)
  ram <- lavaan2RAM("y ~ a*x\na == a0+a1*data.mod")
  fit <- osmasem2(RAM=ram,data=d,cor.analysis=FALSE,replace.constraints=TRUE,run=FALSE)
  m <- OpenMx::omxSetParameters(fit$mx.fit,labels=c("a0","a1"),values=c(.4,-.01),free=FALSE)
  at.mean <- metaSEM:::.osmasem.plot.at.means(m)
  A <- metaSEM:::.osmasem.plot.matrix("Amatrix",at.mean,ram$A)
  expect_equal(A["y","x"],.25)
  expect_equal(m@data@observed$mod,c(10,20))
  ram <- lavaan2RAM("y ~ data.mod*x")
  fit <- osmasem2(RAM=ram,data=d,cor.analysis=FALSE,run=FALSE)
  at.mean <- metaSEM:::.osmasem.plot.at.means(fit$mx.fit)
  A <- metaSEM:::.osmasem.plot.matrix("Amatrix",at.mean,ram$A)
  expect_equal(A["y","x"],15)
})

test_that("cyclic RAM correlation restrictions retain the matrix inverse", {
  C <- matrix(c(1,.2,.2,1),2,dimnames=list(c("x","y"),c("x","y")))
  d <- Cor2DataFrame(rep(list(C),2),n=c(100,100))
  ram <- lavaan2RAM("y ~ a*x\nx ~ b*y\nx ~~ vx*x\ny ~~ vy*y")
  fit <- expect_warning(osmasem2(RAM=ram,data=d,run=FALSE,RE.type.Sigma="Zero",
                                 startvalues=list(a=.1,b=.2,vx=.9,vy=.8)), "Non-recursive")
  m <- fit$mx.fit
  expected <- diag(OpenMx::mxEvalByName("Id_A",m,compute=TRUE) %*%
                     m$Smatrix$values %*% t(OpenMx::mxEvalByName("Id_A",m,compute=TRUE))) - 1
  actual <- metaSEM:::.constraint.values(m,metaSEM:::.constraint.list(m),
                                         OpenMx::omxGetParameters(m))
  expect_equal(as.numeric(actual),as.numeric(expected),tolerance=1e-10)
})

review_moderator_data <- function(mod, correlations=rep(.2,length(mod))) {
  matrices <- lapply(correlations,function(r)
    matrix(c(1,r,r,1),2,dimnames=list(c("x","y"),c("x","y"))))
  data <- Cor2DataFrame(matrices,n=rep(100,length(mod)))
  data$data$mod <- mod
  data
}

test_that("negative moderator powers preserve signs and nonzero slopes", {
  syntax <- paste("y ~ a*x", "x ~~ 1*x", "y ~~ v*y",
                  "a == a0+a1*data.mod",
                  "v == 1-a0^2-2*a0*a1*data.mod-a1^2*data.mod^2",sep="\n")
  ram <- lavaan2RAM(syntax)
  d <- review_moderator_data(c(-2,2),c(.1,.3))
  fit <- osmasem2(RAM=ram,data=d,replace.constraints=TRUE,RE.type.Sigma="Zero",
                  startvalues=list(a0=.2,a1=.05))
  expect_length(fit$mx.fit@constraints,0L)
  expect_equal(unname(coef(fit)[c("a0","a1")]),c(.2,.05),tolerance=1e-5)
  for (mod in c(-2,2)) {
    m <- fit$mx.fit
    m@data@observed[1,"mod"] <- mod
    expect_equal(unname(diag(OpenMx::mxEvalByName("expSigma",m,compute=TRUE))),c(1,1),tolerance=1e-10)
  }
  ## scale() returns a matrix column; preserve its numeric value and precision.
  d$data$mod <- scale(c(1,2))
  model <- osmasem2(RAM=ram,data=d,replace.constraints=TRUE,run=FALSE)
  expect_length(model$mx.fit@constraints,0L)
})

test_that("floating-point moderator patterns give independent variance constraints", {
  mod <- sqrt(c(2,3,5,7,11))
  d <- review_moderator_data(mod,.15+.02*mod)
  ram <- lavaan2RAM(paste("y ~ a*x", "x ~~ 1*x", "y ~~ v*y",
                          "a == a0+a1*data.mod", "v == v0+v1*data.mod+v2*data.mod^2",sep="\n"))
  starts <- list(a0=.15,a1=.02,v0=1-.15^2,v1=-2*.15*.02,v2=-.02^2)
  for (order in list(seq_along(mod),rev(seq_along(mod)))) {
    reordered <- d; reordered$data <- d$data[order,,drop=FALSE]
    fit <- osmasem2(RAM=ram,data=reordered,replace.constraints=TRUE,
                    RE.type.Sigma="Zero",startvalues=starts)
    expect_length(fit$mx.fit@constraints,3L)
    expect_equal(summary(fit)$degreesOfFreedom,3)
    expect_equal(unname(coef(fit)[c("a0","a1")]),c(.15,.02),tolerance=1e-5)
    X <- cbind(1,mod)
    reference <- solve(crossprod(X,X / d$data[[d$vlabels]]))
    expect_equal(unname(vcov(fit)[c("a0","a1"),c("a0","a1")]),
                 unname(reference),tolerance=1e-4)
    for (value in mod) {
      m <- fit$mx.fit; m@data@observed[1,"mod"] <- value
      expect_equal(unname(diag(OpenMx::mxEvalByName("expSigma",m,compute=TRUE))),c(1,1),tolerance=1e-8)
    }
  }
})

test_that("constant division identities do not add variance constraints", {
  d <- review_moderator_data(c(1,2))
  ram <- lavaan2RAM("y ~ a*x\nx ~~ 1*x\ny ~~ v*y\na == b/2\nv == 1-b^2/4")
  fit <- osmasem2(RAM=ram,data=d,replace.constraints=TRUE,RE.type.Sigma="Zero",
                  startvalues=list(b=.4))
  expect_length(fit$mx.fit@constraints,0L)
  expect_equal(unname(coef(fit)["b"]),.4,tolerance=1e-5)
  expect_equal(summary(fit)$degreesOfFreedom,1)
  p <- metaSEM:::.osmasem.polynomials()
  expect_error(p$parse(str2lang("b/0")),"division by zero")
  expect_error(p$parse(str2lang("b/1e-320")),"finite scalars")
})

test_that("small independent restrictions are retained and enforced", {
  d <- review_moderator_data(c(1,2))
  ram <- lavaan2RAM("y ~ a*x\nx ~~ 1*x\ny ~~ v*y\nv == 1-a*a+1e-12*t")
  fit <- osmasem2(RAM=ram,data=d,replace.constraints=TRUE,RE.type.Sigma="Zero",
                  startvalues=list(a=.2,t=1))
  expect_length(fit$mx.fit@constraints,1L)
  expect_equal(unname(coef(fit)["t"]),0,tolerance=1e-7)
  expect_equal(unname(coef(fit)["a"]),.2,tolerance=1e-5)
  expect_equal(vcov(fit)["t","t"],0,tolerance=1e-8)
  p <- metaSEM:::.osmasem.polynomials()
  expect_equal(p$parse(str2lang("a+1e-14*b-a")),c(b=1e-14))
})

test_that("older algebra fits recover plot names from F without RAM", {
  d <- review_moderator_data(c(1,2))
  ram <- lavaan2RAM("y ~ a*x\nx ~~ 1*x\ny ~~ v*y\na == exp(k)\nv == 1-a*a")
  fit <- osmasem2(RAM=ram,data=d,replace.constraints=TRUE,RE.type.Sigma="Zero",
                  startvalues=list(k=log(.2)))
  fit$RAM <- NULL
  A <- metaSEM:::.osmasem.plot.matrix("Amatrix",fit$mx.fit)
  S <- metaSEM:::.osmasem.plot.matrix("Smatrix",fit$mx.fit)
  expect_equal(dimnames(A),rep(list(colnames(fit$mx.fit$Fmatrix$values)),2))
  expect_equal(dimnames(S),dimnames(A))
  if (requireNamespace("semPlot",quietly=TRUE)) {
    grDevices::pdf(tempfile(fileext=".pdf"))
    on.exit(grDevices::dev.off(),add=TRUE)
    expect_error(plot(fit),NA)
  }
})

test_that("replacement diagnostics recognize qualified nested algebras", {
  ram <- lavaan2RAM("y ~ a*x\na == aux.value")
  child <- OpenMx::mxModel("child",
    OpenMx::mxMatrix("Full",1,1,free=TRUE,values=-2,labels="k",name="parameter"),
    OpenMx::mxAlgebraFromString("exp(k)",name="value"))
  aux <- OpenMx::mxModel("aux",child,
    OpenMx::mxAlgebraFromString("child.value",name="value"))
  C <- matrix(c(1,.2,.2,1),2,dimnames=list(c("x","y"),c("x","y")))
  fit <- sem(RAM=ram,Cov=C,numObs=100,mxModel.Args=list(aux))
  expect_equal(as.numeric(OpenMx::mxEvalByName("aux.value",fit$mx.fit)),.2,tolerance=1e-5)
  expect_equal(unname(coef(fit)["a"]),.2,tolerance=1e-5)
})

test_that("rank pruning preserves a small new monomial after elimination", {
  d <- review_moderator_data(c(0,1))
  ram <- lavaan2RAM("y ~ a*x\nx ~~ 1*x\ny ~~ v*y\nv == 1-a*a+t+1e-12*u*data.mod")
  for (order in list(1:2,2:1)) {
    reordered <- d; reordered$data <- d$data[order,,drop=FALSE]
    fit <- osmasem2(RAM=ram,data=reordered,replace.constraints=TRUE,run=FALSE)
    expect_length(fit$mx.fit@constraints,2L)
  }
})
