remaining_review_correlations <- function() {
  C <- matrix(c(1,.2,.1,.2,1,.3,.1,.3,1),3,
              dimnames=list(c("x","y","z"),c("x","y","z")))
  rep(list(C),6)
}

test_that("mxsem plots forward node names and select parameter or RAM labels", {
  skip_if_not_installed("semPlot")
  testthat::local_mocked_bindings(semPaths=function(object,...) object,.package="semPlot")
  ram <- lavaan2RAM("f =~ 1*x + loading*y\nf ~~ latentVar*f\nx ~~ errorX*x\ny ~~ errorY*y",
                    obs.variables=c("y","x"),std.lv=FALSE)
  C <- matrix(c(1,.4,.4,1),2,dimnames=list(c("y","x"),c("y","x")))
  fit <- sem(RAM=ram,Cov=C,numObs=200,run=FALSE)
  ## A fitted matrix specification, without estimating an underidentified CFA.
  named <- plot(fit,manNames=c("Observed Y","Observed X"),latNames="Factor",labels="labels")
  expect_setequal(named@Vars$name,c("Observed Y","Observed X","Factor"))
  loading <- with(named@Pars,lhs == "Factor" & rhs == "Observed Y" & edge == "->")
  expect_equal(named@Pars$label[loading],"loading")
  ram.labels <- plot(fit,labels="RAM")
  expect_false("loading" %in% ram.labels@Pars$label)
  expect_true("loading" %in% plot(fit)@Pars$label)
  expect_error(plot(fit,manNames="one"),"Lengths of 'manNames'")
  expect_error(plot(fit,labels="invalid"),"arg")
  ## Shared and group-specific labels retain their own matrix positions.
  groups <- lavaan2RAM("y ~ c(shared,other)*x",ngroups=2)
  fit <- sem(RAM=groups,Cov=list(C,C),numObs=c(100,100),run=FALSE)
  expect_true("shared" %in% plot(fit,group=1,combine=FALSE)@Pars$label)
  expect_true("other" %in% plot(fit,group=2,combine=FALSE)@Pars$label)
  set.seed(143)
  d <- data.frame(id=rep(1:8,each=4),x=rnorm(32),y=rnorm(32))
  two <- sem(RAM=lavaan2RAM("level: 1\ny ~ bw*x\nlevel: 2\ny ~ bb*x"),
              data=d,cluster="id",run=FALSE)
  within <- plot(two,level="within",combine=FALSE,manNames=c("Y","X"))
  between <- plot(two,level="between",combine=FALSE,manNames=c("Y","X"))
  expect_setequal(within@Vars$name,c("Y","X"))
  expect_setequal(between@Vars$name,c("Y","X"))
  expect_true("bw" %in% within@Pars$label)
  expect_true("bb" %in% between@Pars$label)
})

test_that("cluster WLS summaries forward df adjustment and subclass TSSEM dispatches", {
  matrices <- remaining_review_correlations()
  stage1 <- tssem1(matrices,rep(100,6),method="FEM",cluster=rep(1:2,each=3))
  class(stage1) <- c("review_cluster",class(stage1))
  ram <- lavaan2RAM("y ~ a*x\nz ~ b*y\nx ~~ 1*x")
  stage2 <- tssem2(stage1,RAM=ram)
  expect_s3_class(stage2,"wls.cluster")
  adjusted <- summary(stage2,df.adjustment=2)
  for (i in seq_along(stage2)) {
    direct <- summary(stage2[[i]],df.adjustment=2)
    expect_equal(adjusted[[i]],direct)
    expect_equal(adjusted[[i]]$stat["DF manually adjusted",1],2)
    expect_equal(adjusted[[i]]$stat["DF of target model",1],
                 summary(stage2[[i]])$stat["DF of target model",1]+2)
  }
  single <- stage1[[1]]; class(single) <- c("review_stage1",class(single))
  expect_equal(coef(tssem2(single,RAM=ram)),coef(stage2[[1]]),tolerance=1e-6)
  ## Cluster reruns preserve subclass metadata and rerun each genuine WLS fit.
  class(stage2) <- c("review_stage2",class(stage2)); attr(stage2,"review") <- "keep"
  invisible(capture.output(rerun.fit <- rerun(stage2,extraTries=0)))
  expect_identical(class(rerun.fit),class(stage2))
  expect_identical(attr(rerun.fit,"review"),"keep")
  expect_equal(coef(rerun.fit[[1]]),coef(stage2[[1]]),tolerance=1e-5)
})

test_that("rerun and TSSEM random-effects dispatch recognize subclasses", {
  matrices <- remaining_review_correlations()
  stage1 <- tssem1(matrices,rep(100,6),method="REM",RE.type="Zero")
  original <- stage1
  class(stage1) <- c("review_random",class(stage1))
  ram <- lavaan2RAM("y ~ a*x\nz ~ b*y\nx ~~ 1*x")
  expect_equal(coef(tssem2(stage1,RAM=ram)),coef(tssem2(original,RAM=ram)),tolerance=1e-6)
  invisible(capture.output(rerun.fit <- rerun(stage1,extraTries=0)))
  expect_identical(class(rerun.fit),class(stage1))
  expect_equal(coef(rerun.fit),coef(stage1),tolerance=1e-5)
})

test_that("OSMASEM covariance selections preserve matrix dimensions and R2 accepts subclasses", {
  correlations <- c(.05,.2,.3,.5,.55,.6)
  matrices <- lapply(correlations,function(r)
    matrix(c(1,r,r,1),2,dimnames=list(c("x","y"),c("x","y"))))
  d <- Cor2DataFrame(matrices,n=rep(100,6))
  ram <- lavaan2RAM("y ~ a*x\nx ~~ 1*x")
  first <- osmasem(RAM=ram,data=d,RE.type="Zero")
  second <- osmasem2(RAM=lavaan2RAM("y ~ a*x\nx ~~ 1*x\ny ~~ v*y\nv == 1-a*a"),
                      data=d,RE.type.Sigma="Zero",replace.constraints=TRUE)
  for (fit in list(first,second)) {
    V <- vcov(fit)
    expect_equal(dim(V),c(1L,1L))
    expect_identical(dimnames(V),list("a","a"))
    expect_equal(V[1,1],as.numeric(vcov(fit$mx.fit)),tolerance=1e-6)
    random <- suppressWarnings(vcov(fit,select="random"))
    expect_equal(dim(random),c(0L,0L))
  }
  random <- osmasem(RAM=ram,data=d,RE.type="Diag")
  V <- suppressWarnings(vcov(random,select="random"))
  expect_equal(dim(V),c(1L,1L))
  expect_true(all(is.finite(V)))
  subclass <- random; class(subclass) <- c("review_osmasem",class(random))
  expect_equal(osmasemR2(subclass,subclass),osmasemR2(random,random))
  invisible(capture.output(rerun.fit <- rerun(subclass,autofixtau2=TRUE,extraTries=0)))
  expect_identical(class(rerun.fit),class(subclass))
})

test_that("equal REML variance components retain sample and statistic counts", {
  d <- data.frame(y=c(.1,.3,.2,.4,.6,.5,.2,.4,.3,.7,.9,.8),v=.05,
                  cluster=rep(1:4,each=3))
  equal <- reml3L(y,v,cluster=cluster,data=d,RE.equal=TRUE,run=FALSE)
  reference <- reml3L(y,v,cluster=cluster,data=d,RE.equal=FALSE,run=FALSE)
  expect_equal(equal$fitfunction$numObs,12)
  expect_equal(equal$fitfunction$numStats,11)
  expect_equal(equal$fitfunction$numObs,reference$fitfunction$numObs)
  expect_equal(equal$fitfunction$numStats,reference$fitfunction$numStats)
  reference <- OpenMx::omxSetParameters(reference,labels=c("Tau2_2","Tau2_3"),
                                        newlabels=c("Tau2","Tau2"),values=c(.1,.1))
  equal <- OpenMx::mxRun(equal,silent=TRUE)
  reference <- OpenMx::mxRun(reference,silent=TRUE)
  expect_equal(coef(equal),coef(reference),tolerance=1e-6)
  expect_equal(equal$output$fit,reference$output$fit,tolerance=1e-6)
})

test_that("create.V gives a concise invalid-dimension diagnostic", {
  expect_error(create.V(c("first","second","third"),type="Full"),
               "Length of 'x' \\(3\\) is not a square")
})
