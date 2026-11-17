test_that("explicitly fixed latent means do not fix indicator intercepts", {
  set.seed(108); n <- 500; f <- rnorm(n)
  d <- data.frame(y1=2+f+rnorm(n), y2=3+.8*f+rnorm(n),
                  y3=4+1.2*f+rnorm(n), y4=5+.7*f+rnorm(n))
  syntax <- "f =~ y1+y2+y3+y4"
  explicit <- paste(syntax, "f ~ 0*1", sep="\n")
  base <- sem(RAM=lavaan2RAM(syntax), data=d)
  fit <- sem(RAM=lavaan2RAM(explicit), data=d)
  ref <- lavaan::sem(explicit, data=d, std.lv=TRUE, meanstructure=TRUE, information="observed")
  expect_true(base$mx.fit$output$status$code %in% c(0,1))
  expect_true(fit$mx.fit$output$status$code %in% c(0,1))
  expect_true(lavaan::lavInspect(ref, "converged"))
  expect_true(lavaan::lavInspect(ref, "post.check"))
  expect_equal(fit$mx.fit$output$fit, base$mx.fit$output$fit, tolerance=1e-6)
  expect_equal(fit$mx.fit$output$fit, -2*as.numeric(logLik(ref)), tolerance=1e-5)
  expect_equal(unname(coef(fit)["y2ONf"]), lavaan::lavInspect(ref, "est")$lambda["y2", "f"], tolerance=1e-5)
  reference <- lavaan::parameterEstimates(ref)
  se.ref <- reference$se[reference$lhs == "f" & reference$op == "=~" & reference$rhs == "y2"]
  expect_equal(unname(sqrt(vcov(fit)["y2ONf","y2ONf"])), se.ref, tolerance=1e-5)
  expect_true(all(eigen(vcov(fit),symmetric=TRUE,only.values=TRUE)$values > 0))
  grouped <- lavaan2RAM(explicit, ngroups=2)
  for (ram in grouped) expect_true(all(as.mxMatrix(ram$M)$free[1, match(paste0("y",1:4),colnames(ram$M))]))
})

test_that("free latent means retain user-specified identification", {
  set.seed(118); n <- 400; factor <- rnorm(n, .6)
  d <- data.frame(y1=factor+rnorm(n), y2=2+.8*factor+rnorm(n),
                  y3=3+1.2*factor+rnorm(n), y4=4+.7*factor+rnorm(n))
  syntax <- "f =~ y1+y2+y3+y4\nf ~ mu*1\ny1 ~ 0*1"
  ram <- lavaan2RAM(syntax)
  expect_false(as.mxMatrix(ram$M)$free[1,match("y1",colnames(ram$M))])
  expect_true(all(as.mxMatrix(ram$M)$free[1,match(paste0("y",2:4),colnames(ram$M))]))
  fit <- sem(RAM=ram,data=d)
  ref <- lavaan::sem(syntax,data=d,std.lv=TRUE,meanstructure=TRUE)
  expect_true(fit$mx.fit$output$status$code %in% c(0,1))
  expect_true(lavaan::lavInspect(ref,"converged"))
  expect_true(lavaan::lavInspect(ref,"post.check"))
  expect_equal(fit$mx.fit$output$fit,-2*as.numeric(logLik(ref)),tolerance=1e-5)
  expect_equal(unname(coef(fit)["mu"]),lavaan::parameterEstimates(ref)$est[lavaan::parameterEstimates(ref)$label=="mu"],tolerance=1e-5)
  expect_warning(lavaan2RAM("f =~ y1+y2+y3\nf ~ mu*1"),"mean structure may be unidentified")
  marker <- lavaan2RAM("f =~ 1*y\nf ~ mu*1")
  expect_false(as.mxMatrix(marker$M)$free[1,match("y",colnames(marker$M))])
  expect_warning(lavaan2RAM("f =~ 1*y\nf ~ mu*1\ny ~ m*1"), "mean structure may be unidentified")
})

test_that("cluster summary reconstruction retains small variances and likelihood", {
  S <- diag(2); dimnames(S) <- list(c("x","y"),c("x","y"))
  mats <- rep(list(S),4); names(mats) <- paste0("C",1:4)
  d <- .expand.cluster.summaries(mats,rep(list(c(x=0,y=0)),4),rep(5,4),"id")
  d$y <- d$y*1e-6
  pieces <- split(d[,c("x","y")],d$id)
  covs <- lapply(pieces,cov); means <- lapply(pieces,colMeans)
  rebuilt <- .expand.cluster.summaries(covs,means,rep(5,4),"id")
  expect_equal(var(rebuilt$y)/var(d$y),1,tolerance=1e-10)
  ram <- lavaan2RAM("level: 1\nx ~~ 1*x\ny ~~ 0.000000000001*y\nlevel: 2\nx ~~ 0.5*x\ny ~~ 0.0000000000005*y\nx ~ 0*1\ny ~ 0*1")
  raw <- sem(RAM=ram,data=d,cluster="id",run=FALSE)$mx.fit
  summaries <- sem(RAM=ram,Cov=covs,means=means,numObs=rep(5,4),cluster="id",run=FALSE)$mx.fit
  raw <- OpenMx::mxRun(raw,useOptimizer=FALSE,silent=TRUE)
  summaries <- OpenMx::mxRun(summaries,useOptimizer=FALSE,silent=TRUE)
  expect_equal(raw$output$fit,summaries$output$fit,tolerance=1e-9)
  singular <- matrix(1,2,2,dimnames=dimnames(S))
  out <- .expand.cluster.summaries(list(singular),list(c(x=0,y=0)),2,"id")
  expect_equal(cov(out[,c("x","y")]),singular,tolerance=1e-10)
  bad <- S; bad[2,2] <- -1e-12
  expect_error(.expand.cluster.summaries(list(bad),list(c(x=0,y=0)),5,"id"),"nonnegative variances")
  tiny <- diag(c(1e-12,1)); dimnames(tiny) <- dimnames(S)
  expect_error(.expand.cluster.summaries(list(tiny),list(c(x=1e12,y=0)),5,"id"),
               "cannot be reproduced accurately after adding its means")
})

test_that("covariance-only likelihood intervals use fitted parameters", {
  C <- matrix(c(1,.4,.4,1),2,dimnames=list(c("x","y"),c("x","y")))
  ram <- lavaan2RAM("y ~ b*x",obs.variables=c("x","y"))
  fit <- sem(RAM=ram,Cov=C,numObs=200,intervals.type="LB")
  expect_null(fit$error)
  expect_true(fit$mx.fit$output$status$code %in% c(0,1))
  expect_true("b" %in% rownames(fit$mx.fit$output$confidenceIntervals))
  expect_false(any(grepl("mean",names(fit$mx.fit@intervals))))
  with.means <- sem(RAM=ram,Cov=C,means=c(x=1,y=2),numObs=200,intervals.type="LB")
  expect_true(all(c("xmean","ymean") %in% names(with.means$mx.fit@intervals)))
})

test_that("SEM plots preserve F order, zero paths and parameter metadata", {
  skip_if_not_installed("semPlot")
  testthat::local_mocked_bindings(semPaths=function(object,...) object,.package="semPlot")
  C <- matrix(c(1,.4,.4,1),2,dimnames=list(c("x","y"),c("x","y")))
  ram <- lavaan2RAM("y ~ start(0)*b*x\nx ~~ 1*x\ny ~ my*1",obs.variables=c("x","y"))
  ram$F <- ram$F[2:1,,drop=FALSE]
  fit <- sem(RAM=ram,Cov=C,means=c(x=1,y=2),numObs=200,run=FALSE)
  plot <- plot(fit,manNames=c("Y name","X name"))
  path <- subset(plot@Pars,lhs=="X name" & rhs=="Y name" & edge=="->")
  expect_equal(nrow(path),1L); expect_equal(path$est,0); expect_false(path$fixed)
  variance <- subset(plot@Pars,lhs=="X name" & rhs=="X name" & edge=="<->")
  expect_true(variance$fixed)
  mean <- subset(plot@Pars,rhs=="Y name" & edge=="int")
  expect_equal(mean$est,0); expect_false(mean$fixed)
  ## Repeated labels inside one model must receive the same parameter ID.
  same <- lavaan2RAM("y ~ a*x\nz ~ a*x\nx ~~ 1*x")
  C3 <- diag(3); dimnames(C3)<-list(c("y","z","x"),c("y","z","x"))
  p <- plot(sem(RAM=same,Cov=C3,numObs=200,run=FALSE))
  paths <- p@Pars[p@Pars$label=="a" & p@Pars$edge=="->",]
  expect_equal(nrow(paths),2L); expect_equal(length(unique(paths$par)),1L)
  substituted <- lavaan2RAM("y ~ a*x\nz ~ b*x\na == b")
  p <- plot(sem(RAM=substituted,Cov=C3,numObs=200,replace.constraints=TRUE,run=FALSE))
  paths <- p@Pars[p@Pars$edge == "->", ]
  expect_equal(nrow(paths),2L); expect_true(all(paths$label == "b"))
  expect_equal(length(unique(paths$par)),1L)
  nonlinear <- lavaan2RAM("y ~ a*x\nz ~ b*x\na == 2*b")
  p <- plot(sem(RAM=nonlinear,Cov=C3,numObs=200,replace.constraints=TRUE,run=FALSE))
  derived <- subset(p@Pars,lhs=="x" & rhs=="y" & edge=="->")
  expect_equal(derived$par,0); expect_false(derived$fixed)
  expect_equal(derived$est,0.2)
  constrained <- lavaan2RAM("y ~ b*x\nb == 0.2",obs.variables=c("x","y"))
  p <- plot(sem(RAM=constrained,Cov=C,numObs=200,replace.constraints=TRUE,run=FALSE))
  fixed.path <- subset(p@Pars,lhs=="x" & rhs=="y" & edge=="->")
  expect_true(fixed.path$fixed); expect_equal(fixed.path$est,0.2)
  group.ram <- lavaan2RAM("y ~ c(b1,b2)*x",ngroups=2,M.start=0.2)
  group.fit <- sem(RAM=group.ram,Cov=list(A=C,B=C),numObs=c(100,100),run=FALSE)
  group.plot <- plot(group.fit,group="A")
  expect_false(any(group.plot@Pars$edge == "int"))
})

test_that("zero-valued free paths are drawn in estimate diagrams", {
  skip_if_not_installed("semPlot")
  C <- diag(2); dimnames(C) <- list(c("y","x"),c("y","x"))
  ram <- lavaan2RAM("y ~ start(0)*b*x\nx ~~ 1*x")
  fit <- sem(RAM=ram,Cov=C,numObs=200,run=FALSE)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(),add=TRUE)
  p <- plot(fit)
  labels <- p$graphAttributes$Nodes$labels
  path <- labels[p$Edgelist$from] == "x" & labels[p$Edgelist$to] == "y" & !p$Edgelist$bidirectional
  expect_equal(sum(path),1L)
  expect_equal(p$graphAttributes$Edges$labels[path],"0.00")
})

test_that("auxiliary OpenMx submodels are not plotted as groups", {
  skip_if_not_installed("semPlot")
  C <- matrix(c(1,.4,.4,1),2,dimnames=list(c("x","y"),c("x","y")))
  ram <- lavaan2RAM("y ~ c(b1,b2)*x",ngroups=2)
  fit <- sem(RAM=ram,Cov=list(A=C,B=C),numObs=c(100,100),
             mxModel.Args=list(OpenMx::mxModel("aux")),run=FALSE)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(),add=TRUE)
  expect_no_error(plot(fit,combine=TRUE))
  expect_no_error(plot(fit,group="B"))
})

test_that("explicit fixed latent variances are retained with std.lv", {
  ram <- lavaan2RAM("f =~ x1+x2+x3\nf ~~ 2*f",std.lv=TRUE)
  expect_equal(as.numeric(ram$S["f","f"]),2)
})
