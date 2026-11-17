## Evaluate "expr", collecting all of its warnings. Robust covariances at a
## bound must warn about the bound, and about the OpenMx status exactly when
## the fit's status code is neither 0 nor 1.
expect_robust_warnings <- function(expr, fit, pattern="at a bound") {
  messages <- character(0)
  value <- withCallingHandlers(expr, warning=function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_true(any(grepl(pattern, messages)))
  status <- fit$mx.fit@output$status$code
  expect_identical(any(grepl("OpenMx status code", messages)),
                   !(status %in% c(0, 1)))
  value
}

test_that("positive definiteness rejects singular and asymmetric matrices", {
  expect_false(is.pd(matrix(0, 2, 2)))
  expect_false(is.pd(matrix(c(1, 0, 1, 1), 2, 2)))
  expect_true(is.pd(matrix(2, 1, 1)))
  expect_identical(is.pd(matrix(NA_real_, 1, 1)), NA)
})

test_that("bare parameter labels are accepted", {
  m <- as.mxMatrix(matrix("loading", 1, 1), name="A")
  expect_true(m$free[1, 1])
  expect_identical(m$labels[1, 1], "loading")
  unlabeled <- as.mxMatrix(matrix("5*", 1, 1), name="B")
  expect_true(unlabeled$free[1, 1])
  expect_equal(unlabeled$values[1, 1], 5)
  expect_true(is.na(unlabeled$labels[1, 1]))
})

test_that("REML uses matrix starting values and bounds", {
  y <- matrix(seq_len(12)/10, 6, 2)
  v <- matrix(rep(c(1, 0, 1), 6), 6, 3, byrow=TRUE)
  starts <- matrix(c(.1, .02, .02, .2), 2)
  bounds <- matrix(c(0, NA, NA, 0), 2)
  model <- reml(y, v, RE.startvalues=starts, RE.lbound=bounds, run=FALSE)
  expect_equal(model$Tau$values[1:2, 1:2], starts)
  expect_equal(model$Tau$lbound[1:2, 1:2], bounds)
  expect_error(reml(y, v, RE.startvalues=matrix(1:4, 2), run=FALSE),
               "must be symmetric")
})

test_that("OSMASEM missingness counts require both statistics to be present", {
  x <- list(data=data.frame(a=c(1, NA, 1), b=c(NA, 1, 1)),
            ylabels=c("a", "b"))
  expect_equal(unname(pattern.na(x, show.na=FALSE, type="osmasem")),
               matrix(c(2, 1, 1, 2), 2))
  expect_equal(unname(pattern.na(x, show.na=TRUE, type="osmasem")),
               matrix(c(1, 0, 0, 1), 2))
})

test_that("fully missing covariance entries cannot be averaged", {
  one <- matrix(c(1, NA, NA, 1), 2,
                dimnames=list(c("a", "b"), c("a", "b")))
  expect_error(Cor2DataFrame(list(one, one), n=c(100, 100)),
               "missing from every covariance matrix")
})

test_that("one-variable covariance asymptotic variance follows normal theory", {
  expect_equal(unname(asyCov(matrix(2, 1, 1), n=100,
                             cor.analysis=FALSE)), matrix(.08, 1, 1))
  expect_error(asyCov(matrix(1, 1, 1), n=100), "At least two variables")
})

test_that("path levels handle isolated nodes and reject cycles", {
  none <- create.vechsR(matrix(0, 2, 2), diag(2))
  expect_true(all(c("Amatrix", "Smatrix", "vechsR") %in% names(none)))
  a <- matrix(c(0, "0.1*b", "0.1*a", 0), 2)
  expect_error(metaSEM:::.pathLevels(as.mxMatrix(a, name="A")),
               "directed cycle")
})

test_that("flat SEM rejects group and retains cluster identifiers", {
  ram <- lavaan2RAM("y ~~ v*y", obs.variables="y")
  d <- data.frame(y=rnorm(10), g=rep(1:2, each=5))
  expect_error(sem(RAM=ram, data=d, group="g", run=FALSE),
               "requires a multiple-group RAM")
  expect_identical(sem(RAM=ram, data=d, cluster="g", run=FALSE)$cluster,
                   "g")
  expect_error(sem(RAM=ram, data=d, cluster="missing", run=FALSE),
               "must name a complete column")
})

test_that("single-level clustered SEM agrees with an independent sandwich", {
  skip_if_not_installed("sandwich")
  set.seed(44)
  id <- rep(seq_len(25), each=6)
  x <- rnorm(length(id))
  d <- data.frame(id=id, x=x,
                  y=.6*x + rep(rnorm(25), each=6) + rnorm(length(id)))
  ram <- lavaan2RAM("y ~ b*x")
  fit <- sem(RAM=ram, data=d, cluster="id")
  ordinary <- sem(RAM=ram, data=d)
  direct <- clusterRobustSE(fit$mx.fit, cluster="id", details=TRUE)
  reference <- sandwich::vcovCL(lm(y ~ x, data=d), cluster=d$id,
                                type="HC0", cadjust=TRUE)
  expect_equal(unname(direct$SE["b"]),
               unname(sqrt(diag(reference))["x"]), tolerance=1e-4)
  expect_equal(vcov(fit, robust=TRUE), direct$cov, tolerance=1e-9)
  expect_equal(coef(fit), coef(ordinary), tolerance=1e-8)
  expect_equal(unname(summary(fit, robust=TRUE)$coefficients["b", "Std.Error"]),
               unname(direct$SE["b"]), tolerance=1e-8)
  corrected <- clusterRobustSE(fit$mx.fit, cluster="id", details=TRUE,
                               smallSample=TRUE)
  expect_equal(direct$smallSampleFactor, 25/24)
  expect_equal(corrected$smallSampleFactor,
               25/24 * (nrow(d)-1)/(nrow(d)-length(coef(fit))))
  expect_equal(corrected$cov,
               direct$cov * corrected$smallSampleFactor/direct$smallSampleFactor,
               tolerance=1e-8)
  expect_equal(vcov(fit, robust=TRUE, smallSample=TRUE), corrected$cov)
  expect_error(clusterRobustSE(fit$mx.fit, cluster=rep(1, nrow(d))),
               "at least two clusters")
  missing.data <- d
  missing.data$y[c(2, 8, 14)] <- NA_real_
  missing.fit <- sem(RAM=ram, data=missing.data, cluster="id")
  expect_true(is.finite(vcov(missing.fit, robust=TRUE)["b", "b"]))
})

test_that("flat multiple-group cluster scores follow contributing groups", {
  set.seed(45)
  d <- data.frame(g=rep(c("B", "A"), each=60),
                  id=rep(rep(seq_len(12), each=5), 2),
                  x=rnorm(120), y=rnorm(120))
  ram <- lavaan2RAM("y ~ c(b1,b2)*x", ngroups=2)
  fit <- sem(RAM=ram, data=d, group="g", cluster="id")
  ids <- unlist(lapply(fit$mx.fit@submodels[
    fit$mx.fit@fitfunction$groups],
    function(x) x@data@observed$id), use.names=FALSE)
  expect_equal(vcov(fit, robust=TRUE),
               clusterRobustSE(fit$mx.fit, cluster=ids, details=TRUE)$cov)
  expect_true(all(is.finite(vcov(fit, robust=TRUE))))
  ## Independent reference: lavaan MLR with clusters nested in groups. lavaan
  ## multiplies each group's meat by G_g/(G_g-1) (12/11), whereas the
  ## clustered sandwich uses the total number of clusters (24/23).
  d$id <- paste(d$g, d$id, sep=".")
  model <- "y ~ c(b1,b2)*x + c(i1,i2)*1
  y ~~ c(e1,e2)*y
  x ~~ c(v1,v2)*x
  x ~ c(m1,m2)*1"
  fit <- sem(RAM=lavaan2RAM(model, ngroups=2), data=d, group="g",
             cluster="id")
  reference <- lavaan::sem(model, data=d, group="g", cluster="id",
                           estimator="MLR", fixed.x=FALSE)
  expected <- lavaan::parameterEstimates(reference)
  expected <- expected[expected$label != "", ]
  expected <- setNames(expected$se, expected$label)
  actual <- sqrt(diag(vcov(fit, robust=TRUE)))
  expect_setequal(names(actual), names(expected))
  expect_equal(actual*sqrt((12/11)/(24/23)), expected[names(actual)],
               tolerance=1e-4)
})

test_that("cluster sandwich uses G/(G-1), as lavaan and imxRobustSE() do", {
  set.seed(46)
  n <- 40
  d <- data.frame(id=seq_len(n), x=rnorm(n))
  d$y <- .3*d$x + rnorm(n)*exp(.5*d$x)
  ram <- lavaan2RAM("y ~ b*x", obs.variables=c("y", "x"))
  ## One row per cluster: the same as the observation-level sandwich.
  expect_equal(vcov(sem(RAM=ram, data=d, cluster="id"), robust=TRUE),
               vcov(sem(RAM=ram, data=d), robust=TRUE), tolerance=1e-8)
  ## Clusters of size 4 match lavaan's single-group cluster-robust SEs.
  d$id <- rep(seq_len(n/4), each=4)
  model <- "y ~ b*x + i*1
  y ~~ e*y
  x ~~ v*x
  x ~ m*1"
  fit <- sem(RAM=lavaan2RAM(model), data=d, cluster="id")
  reference <- lavaan::sem(model, data=d, cluster="id", estimator="MLR",
                           fixed.x=FALSE)
  expected <- lavaan::parameterEstimates(reference)
  expected <- expected[expected$label != "", ]
  expected <- setNames(expected$se, expected$label)
  actual <- sqrt(diag(vcov(fit, robust=TRUE)))
  expect_equal(actual, expected[names(actual)], tolerance=1e-4)
})

test_that("WLS subsets the asymptotic covariance in requested variable order", {
  cov <- diag(4)
  dimnames(cov) <- list(c("a", "b", "c", "d"), c("a", "b", "c", "d"))
  acov <- diag(1:6)
  model <- wls(cov, acov, n=100, Smatrix=diag(3),
               subset.variables=c("c", "a", "b"), run=FALSE)
  expect_equal(diag(model$invAcov$values), c(1/2, 1/4, 1))
})

test_that("fixed-effect TSSEM aligns named studies before pooling", {
  variables <- c("a", "b", "c")
  r <- matrix(c(1, .1, .2, .1, 1, .3, .2, .3, 1), 3,
              dimnames=list(variables, variables))
  permuted <- r[c("b", "a", "c"), c("b", "a", "c")]
  fit <- tssem1FEM(list(r, permuted), n=c(200, 200))
  expect_equal(coef(fit), r, tolerance=1e-5)
  expect_equal(fit$data[[1]], fit$data[[2]])

  differently.ordered.axes <- r[c("b", "a", "c"), variables]
  fit.axes <- tssem1FEM(list(r, differently.ordered.axes),
                        n=c(200, 200), run=FALSE)
  expect_equal(fit.axes@submodels[[1]]@data@observed,
               fit.axes@submodels[[2]]@data@observed)

  wrong.names <- r
  dimnames(wrong.names) <- list(c("a", "b", "d"), c("a", "b", "d"))
  expect_error(tssem1FEM(list(r, wrong.names), c(200, 200), run=FALSE),
               "variable names differ")
  expect_error(tssem1FEM(list(r, unname(r)), c(200, 200), run=FALSE),
               "every study matrix needs")
})

test_that("fixed-effect TSSEM rejects unequal study dimensions clearly", {
  r <- diag(3)
  expect_error(tssem1FEM(list(r, diag(2)), c(100, 100), run=FALSE),
               "not the same square size")
  expect_error(tssem1FEM(list(r), c(100, 100), run=FALSE),
               "Lengths of 'Cov' and 'n'")
})

test_that("matrix stacking and random-effect TSSEM align named studies", {
  variables <- c("a", "b", "c")
  make.cor <- function(ab, ac, bc) {
    matrix(c(1, ab, ac, ab, 1, bc, ac, bc, 1), 3,
           dimnames=list(variables, variables))
  }
  studies <- list(s1=make.cor(.15, .20, .30),
                  s2=make.cor(.10, .25, .35),
                  s3=make.cor(.20, .18, .28),
                  s4=make.cor(.12, .22, .32))
  reordered <- studies
  reordered$s2 <- studies$s2[c("b", "a", "c"), c("b", "a", "c")]
  reordered$s4 <- studies$s4[c("c", "a", "b"), c("c", "a", "b")]

  expect_identical(list2matrix(reordered), list2matrix(studies))
  expect_identical(rownames(list2matrix(reordered)), names(studies))
  expect_identical(list2matrix(reordered, diag=TRUE),
                   list2matrix(studies, diag=TRUE))
  expect_equal(coef(tssem1(reordered, n=c(160, 180, 170, 190))),
               coef(tssem1(studies, n=c(160, 180, 170, 190))),
               tolerance=1e-5)

  wrong <- studies$s2
  dimnames(wrong) <- list(c("a", "b", "d"), c("a", "b", "d"))
  expect_error(list2matrix(list(studies$s1, wrong)), "variable names differ")
  expect_error(tssem1(list(studies$s1, unname(studies$s2)), c(160, 180)),
               "every study matrix needs")
})

test_that("Cor2DataFrame aligns effect sizes and sampling covariances", {
  variables <- c("a", "b", "c")
  one <- matrix(c(1, .2, .3, .2, 1, .4, .3, .4, 1), 3,
                dimnames=list(variables, variables))
  two <- matrix(c(1, .1, .25, .1, 1, .35, .25, .35, 1), 3,
                dimnames=list(variables, variables))
  permuted <- two[c("b", "a", "c"), c("b", "a", "c")]
  expected <- Cor2DataFrame(list(one, two), n=c(120, 140))
  actual <- Cor2DataFrame(list(one, permuted), n=c(120, 140))
  expect_equal(actual, expected)

  expected.cov <- Cor2DataFrame(list(one, two), n=c(120, 140),
                                cor.analysis=FALSE)
  actual.cov <- Cor2DataFrame(list(one, permuted), n=c(120, 140),
                              cor.analysis=FALSE)
  expect_equal(actual.cov, expected.cov)
})

test_that("fixed-effect TSSEM explains partially missing correlations", {
  complete <- matrix(c(1, .2, .3, .2, 1, .4, .3, .4, 1), 3)
  partial <- complete
  partial[1, 2] <- partial[2, 1] <- NA_real_
  expect_error(tssem1FEM(list(complete, partial), n=c(120, 140)),
               "requires complete correlations among variables present")
})

test_that("OSMASEM surfaces OpenMx run failures", {
  m <- create.vechsR(matrix(0, 2, 2), diag(2))
  tau <- create.Tau2(no.var=1, RE.type="Zero")
  cor <- matrix(c(1, .2, .2, 1), 2,
                dimnames=list(c("a", "b"), c("a", "b")))
  d <- Cor2DataFrame(list(cor), n=100)
  expect_s3_class(osmasem(Mmatrix=m, Tmatrix=tau, data=d, run=FALSE),
                  "osmasem")
  expect_error(osmasem(Mmatrix=m, Tmatrix=tau, data=d,
                       mxRun.Args=list(noSuchArgument=TRUE)),
               "does not accept")
})

test_that("OSMASEM saturated helper accepts Cor2DataFrame input directly", {
  cor <- matrix(c(1, .2, .2, 1), 2,
                dimnames=list(c("a", "b"), c("a", "b")))
  d <- Cor2DataFrame(list(cor, cor, cor), n=c(100, 120, 110))
  tau <- create.Tau2(no.var=1, RE.type="Zero")
  sat <- metaSEM:::.osmasemSatIndMod(model="Saturated", Tmatrix=tau,
                                     data=d, run=FALSE)
  ind <- metaSEM:::.osmasemSatIndMod(model="Independence", Tmatrix=tau,
                                     data=d, run=FALSE)
  expect_s4_class(sat, "MxModel")
  expect_s4_class(ind, "MxModel")
  expect_equal(nrow(sat@data@observed), 3)
  expect_equal(sat@expectation@dims, d$ylabels)
})

test_that("two-level regression on a cluster covariate agrees with lavaan", {
  set.seed(1)
  id <- rep(seq_len(20), each=10)
  z <- rep(rnorm(20), each=10)
  y <- .5*z + rep(rnorm(20), each=10) + rnorm(200)
  d <- data.frame(id=id, z=z, y=y)
  model <- "level: 1
  y ~~ y
level: 2
  y ~ b*z
  y ~~ y
  z ~~ z"
  reference <- lavaan::sem(model, data=d, cluster="id")
  fit <- sem(RAM=lavaan2RAM(model), data=d, cluster="id")
  expected <- lavaan::parameterEstimates(reference)
  expect_equal(unname(coef(fit)["b"]),
               expected$est[expected$lhs=="y" & expected$op=="~" & expected$rhs=="z"],
               tolerance=1e-4)
  expect_equal(-fit$mx.fit@output$Minus2LogLikelihood/2,
               as.numeric(lavaan::fitMeasures(reference, "logl")), tolerance=1e-3)
  d$z[2] <- d$z[2]+1
  expect_error(sem(RAM=lavaan2RAM(model), data=d, cluster="id", run=FALSE),
               "constant within each cluster")
})

test_that("combined group and two-level estimates agree with lavaan", {
  set.seed(12)
  g <- rep(c("A", "B"), each=100)
  id <- rep(seq_len(40), each=5)
  y <- c(rep(rnorm(20), each=5)+rnorm(100),
         rep(rnorm(20, 1), each=5)+rnorm(100))
  d <- data.frame(g=g, id=id, y=y)
  model <- "group: 1
  level: 1
    y ~~ y
  level: 2
    y ~~ y
group: 2
  level: 1
    y ~~ y
  level: 2
    y ~~ y"
  reference <- lavaan::sem(model, data=d, group="g", cluster="id")
  fit <- sem(RAM=lavaan2RAM(model, ngroups=2), data=d,
             group="g", cluster="id")
  expected <- lavaan::parameterEstimates(reference)
  expect_equal(unname(coef(fit)["yWITHy_within_g1"]),
               expected$est[expected$group==1 & expected$level==1 &
                            expected$op=="~~"], tolerance=1e-3)
  expect_equal(unname(coef(fit)["yWITHy_between_g2"]),
               expected$est[expected$group==2 & expected$level==2 &
                            expected$op=="~~"], tolerance=1e-3)
  expect_equal(-fit$mx.fit@output$Minus2LogLikelihood/2,
               as.numeric(lavaan::fitMeasures(reference, "logl")), tolerance=1e-3)
})

test_that("combined group and level cross-group constraints can be substituted", {
  set.seed(21)
  g <- rep(c("A", "B"), each=100)
  id <- rep(seq_len(40), each=5)
  y <- c(rep(rnorm(20), each=5)+rnorm(100),
         rep(rnorm(20, 1), each=5)+rnorm(100))
  d <- data.frame(g=g, id=id, y=y)
  model <- "group: 1\n level: 1\n y ~~ v1*y\n level: 2\n y ~~ y\n group: 2\n level: 1\n y ~~ v2*y\n level: 2\n y ~~ y\nv2 == 2*v1"
  ram <- lavaan2RAM(model, ngroups=2)
  fit <- sem(RAM=ram, data=d, group="g", cluster="id",
             replace.constraints=TRUE)
  expect_true(is.finite(fit$mx.fit@output$Minus2LogLikelihood))
  expect_true("v1" %in% names(coef(fit)))
  expect_false("v2" %in% names(coef(fit)))
})

test_that("multiple-group covariance summaries reproduce raw-data estimates", {
  set.seed(42)
  d <- data.frame(y=rnorm(200), x=rnorm(200),
                  g=rep(c("a", "b"), each=100))
  ram <- lavaan2RAM("y ~ c(b,b)*x + c(a1,a2)*1", ngroups=2)
  raw <- sem(RAM=ram, data=d, group="g")
  groups <- split(d, d$g)
  covs <- lapply(groups, function(x) cov(x[, c("y", "x")]))
  means <- lapply(groups, function(x) colMeans(x[, c("y", "x")]))
  summary.fit <- sem(RAM=ram, Cov=covs, means=means, numObs=c(100, 100))
  expect_equal(coef(summary.fit), coef(raw), tolerance=1e-4)
  no.means <- sem(RAM=ram, Cov=covs, numObs=c(100, 100))
  expect_false(any(grepl("mean|^a[12]$", names(coef(no.means)))))
})

test_that("two-level cluster-robust standard errors agree with lavaan", {
  skip_on_cran()
  ## Unbalanced clusters, heteroscedastic errors, a predictor at both levels,
  ## and a cluster-level covariate; parameters matched by label.
  set.seed(123)
  G <- 30
  id <- rep(seq_len(G), sample(3:10, G, replace=TRUE))
  z <- rnorm(G)
  xb <- rnorm(G)
  x <- xb[id] + rnorm(length(id))
  y <- 0.5*x + 0.3*z[id] + 0.4*xb[id] + rnorm(G)[id] +
    rnorm(length(id))*exp(0.3*z[id])
  d <- data.frame(id=id, x=x, y=y, z=z[id])
  model <- "level: 1
    y ~ bw*x
    y ~~ ew*y
    x ~~ vxw*x
  level: 2
    y ~ bb*x + bz*z
    y ~~ eb*y
    x ~~ vxb*x
    z ~~ vz*z
    x ~~ cxz*z
    y ~ my*1
    x ~ mx*1
    z ~ mz*1"
  reference <- lavaan::sem(model, data=d, cluster="id", fixed.x=FALSE,
                           se="robust.huber.white")
  expected <- lavaan::parameterEstimates(reference)
  expected <- expected[expected$label != "" & !duplicated(expected$label), ]
  expected <- setNames(expected$se, expected$label)
  fit <- sem(RAM=lavaan2RAM(model), data=d, cluster="id")
  actual <- sqrt(diag(vcov(fit, robust=TRUE)))
  expect_setequal(names(actual), names(expected))
  expect_equal(actual, expected[names(actual)], tolerance=1e-4)
  expect_equal(summary(fit, robust=TRUE)$coefficients[names(actual), "Std.Error"],
               unname(actual), tolerance=1e-8)
  expect_error(clusterRobustSE(fit$mx.fit, cluster=id),
               "Nested two-level models are unsupported")
})

test_that("two-level robust covariance treats parameters at a bound as fixed", {
  skip_on_cran()
  set.seed(3)
  id <- rep(seq_len(25), each=6)
  d <- data.frame(id=id, x=rnorm(150), y=rnorm(150))
  d$y <- d$y - ave(d$y, d$id)
  model <- "level: 1
    y ~ bw*x
    y ~~ ew*y
  level: 2
    y ~~ tb*y
    y ~ my*1"
  fit <- sem(RAM=lavaan2RAM(model), data=d, cluster="id", lbound=list(tb=0))
  expect_equal(unname(coef(fit)["tb"]), 0, tolerance=1e-6)
  covariance <- expect_robust_warnings(vcov(fit, robust=TRUE), fit, "at a bound \\(tb\\)")
  expect_true(all(is.na(covariance["tb", ])))

  ## Reference: the model refitted with tb fixed at 0, scored per cluster by
  ## R-level central differences of each cluster's -2LL.
  fixed <- sem(RAM=lavaan2RAM(sub("tb\\*y", "0*y", model)), data=d,
               cluster="id")
  mx <- fixed$mx.fit
  theta <- coef(mx)
  between <- mx@submodels[[1]]
  within.data <- mx@data@observed
  between.data <- between@data@observed
  scores <- t(vapply(unique(as.character(within.data$id)), function(i) {
    one <- mxModel(mx,
      mxModel(between, mxData(between.data[as.character(between.data$id)==i, ,
                                           drop=FALSE],
                              type="raw", primaryKey="id")),
      mxData(within.data[as.character(within.data$id)==i, , drop=FALSE],
             type="raw"))
    f <- function(j, h) mxRun(omxSetParameters(one, labels=names(theta)[j],
                                               values=theta[j]+h),
                              useOptimizer=FALSE, silent=TRUE)@output$fit
    vapply(seq_along(theta), function(j) {
      h <- 1e-5*max(1, abs(theta[j]))
      (f(j, h)-f(j, -h))/(2*h)
    }, numeric(1))
  }, numeric(length(theta))))
  bread <- vcov(mx)[names(theta), names(theta)]/2
  expected <- bread %*% crossprod(scores) %*% bread
  expect_equal(covariance[names(theta), names(theta)], expected,
               tolerance=1e-4, ignore_attr=TRUE)
})

test_that("fast row gradients reproduce OpenMx row gradients", {
  set.seed(8)
  id <- rep(seq_len(15), each=4)
  d <- data.frame(id=id, x=rnorm(60), z=rnorm(15)[id])
  d$y <- 0.3*d$x + rnorm(15)[id] + rnorm(60)
  d$y[c(3, 17, 40)] <- NA_real_
  fit <- sem(RAM=lavaan2RAM("y ~ b*x + c*z"), data=d, cluster="id")
  fast <- metaSEM:::.row.m2ll.differences(fit$mx.fit, coef(fit))
  reference <- OpenMx::imxRowGradients(fit$mx.fit, robustSE=TRUE)
  expect_equal(unname(fast), unname(as.matrix(reference)),
               tolerance=1e-6)

  ## Observation-level covariance reproduces imxRobustSE().
  plain <- sem(RAM=lavaan2RAM("y ~ b*x + c*z"), data=d)
  expect_equal(vcov(plain, robust=TRUE),
               suppressMessages(imxRobustSE(plain$mx.fit, details=TRUE)$cov),
               tolerance=1e-6, ignore_attr=TRUE)

  ## The imxRowGradients() fallback gives the same covariance.
  fast.cov <- clusterRobustSE(fit$mx.fit, cluster="id", details=TRUE)$cov
  testthat::local_mocked_bindings(
    .row.m2ll.differences=function(...) stop("unavailable"),
    .package="metaSEM")
  expect_equal(clusterRobustSE(fit$mx.fit, cluster="id", details=TRUE)$cov,
               fast.cov, tolerance=1e-6)
})

test_that("single-level robust covariance treats parameters at a bound as fixed", {
  skip_on_cran()
  set.seed(7)
  id <- rep(seq_len(60), sample(1:5, 60, replace=TRUE))
  v <- runif(length(id), 0.02, 0.1)
  d <- data.frame(id=id, v=v, yi=0.3 + 0.6*rnorm(length(id), sd=sqrt(v)))
  model <- "f =~ 1*yi
    f ~~ tau2*f
    yi ~~ data.v*yi
    f ~ mu*1"
  fit <- sem(RAM=lavaan2RAM(model, obs.variables="yi", std.lv=FALSE),
             data=d, cluster="id", lbound=list(tau2=0))
  expect_equal(unname(coef(fit)["tau2"]), 0)
  covariance <- expect_robust_warnings(vcov(fit, robust=TRUE), fit, "at a bound \\(tau2\\)")
  expect_true(is.na(covariance["tau2", "tau2"]))
  fixed <- sem(RAM=lavaan2RAM(sub("tau2\\*f", "0*f", model),
                              obs.variables="yi", std.lv=FALSE),
               data=d, cluster="id")
  expect_equal(covariance["mu", "mu"], vcov(fixed, robust=TRUE)["mu", "mu"],
               tolerance=1e-6)
  se <- expect_robust_warnings(clusterRobustSE(fit$mx.fit, cluster="id"), fit)
  expect_true(is.na(se["tau2"]))
})

test_that("combined group and level robust covariance matches separate group fits", {
  skip_on_cran()
  set.seed(31)
  d <- data.frame(g=rep(c("A", "B"), each=80), id=rep(seq_len(32), each=5),
                  x=rnorm(160))
  d$y <- 0.4*d$x + rep(rnorm(32), each=5) + rnorm(160)
  one <- function(k) gsub("#", k, "level: 1
    y ~ bw#*x
    y ~~ ew#*y
    x ~~ vxw#*x
  level: 2
    y ~~ eb#*y
    x ~~ vxb#*x
    y ~~ cb#*x
    y ~ my#*1
    x ~ mx#*1")
  ram <- lavaan2RAM(paste0("group: A\n", one(1), "\ngroup: B\n", one(2)),
                    ngroups=2)
  fit <- sem(RAM=ram, data=d, group="g", cluster="id")
  covariance <- vcov(fit, robust=TRUE)
  expect_equal(covariance, t(covariance), tolerance=1e-8)
  ## Groups share no parameters, so each block is that group's own fit.
  for (k in 1:2) {
    separate <- sem(RAM=lavaan2RAM(one(k)), data=d[d$g == c("A", "B")[k], ],
                    cluster="id")
    names.k <- names(coef(separate))
    expect_equal(covariance[names.k, names.k], vcov(separate, robust=TRUE),
                 tolerance=1e-3)
  }
  ## Clusters are nested within groups, so a reused ID is ambiguous.
  d$id <- rep(rep(seq_len(16), each=5), 2)
  expect_warning(shared <- sem(RAM=ram, data=d, group="g", cluster="id"),
                 "appear in more than one group")
  expect_error(vcov(shared, robust=TRUE), "nested within groups")
  expect_error(summary(shared, robust=TRUE), "nested within groups")
})

test_that("per-cluster summaries reproduce two-level Gaussian fits", {
  skip_on_cran()
  set.seed(18)
  id <- rep(seq_len(15), each=5)
  d <- data.frame(y=rep(rnorm(15), each=5)+rnorm(75), id=id)
  groups <- split(d, d$id)
  covs <- lapply(groups, function(x) cov(x["y"]))
  means <- lapply(groups, function(x) colMeans(x["y"]))
  n <- rep(5, length(groups))
  model <- "level: 1
    y ~~ y
  level: 2
    y ~~ y"
  ram <- lavaan2RAM(model)
  raw <- sem(RAM=ram, data=d, cluster="id")
  summary.fit <- sem(RAM=ram, Cov=covs, means=means, numObs=n,
                     cluster="id")
  expect_equal(coef(summary.fit), coef(raw), tolerance=1e-5)
  expect_equal(summary.fit$mx.fit@output$Minus2LogLikelihood,
               raw$mx.fit@output$Minus2LogLikelihood, tolerance=1e-5)
  expect_error(sem(RAM=ram, Cov=cov(d["y"]), means=colMeans(d["y"]),
                   numObs=nrow(d), cluster="id"),
               "one covariance matrix")
})

test_that("nested summaries reproduce combined group and level fits", {
  skip_on_cran()
  set.seed(19)
  g <- rep(c("A", "B"), each=75)
  id <- rep(rep(seq_len(15), each=5), 2)
  y <- c(rep(rnorm(15), each=5)+rnorm(75),
         rep(rnorm(15, 1), each=5)+rnorm(75))
  d <- data.frame(g=g, id=id, y=y)
  model <- "group: 1
  level: 1
    y ~~ y
  level: 2
    y ~~ y
group: 2
  level: 1
    y ~~ y
  level: 2
    y ~~ y"
  ram <- lavaan2RAM(model, ngroups=2)
  groups <- split(d, d$g)
  covs <- lapply(groups, function(a)
    lapply(split(a, a$id), function(x) cov(x["y"])))
  means <- lapply(groups, function(a)
    lapply(split(a, a$id), function(x) colMeans(x["y"])))
  n <- lapply(groups, function(a) rep(5, 15))
  expect_warning(raw <- sem(RAM=ram, data=d, group="g", cluster="id"),
                 "appear in more than one group")
  summary.fit <- sem(RAM=ram, Cov=covs, means=means, numObs=n,
                     group="g", cluster="id")
  expect_equal(coef(summary.fit), coef(raw), tolerance=1e-5)
  expect_equal(summary.fit$mx.fit@output$Minus2LogLikelihood,
               raw$mx.fit@output$Minus2LogLikelihood, tolerance=1e-5)
  ## Generated cluster IDs are prefixed by group position because they repeat.
  position <- match(summary.fit$data$g, names(covs))
  expect_true(all(startsWith(as.character(summary.fit$data$id),
                             paste0(position, "."))))
  ## For complete data each cluster score depends only on its summaries.
  d$id <- paste(d$g, d$id, sep=".")
  raw <- sem(RAM=ram, data=d, group="g", cluster="id")
  expect_equal(vcov(summary.fit, robust=TRUE), vcov(raw, robust=TRUE),
               tolerance=1e-4)
})

test_that("generated cluster IDs stay distinct when names contain periods", {
  ## Prefixing by group name would map group "a.b"/cluster "c" and group
  ## "a"/cluster "b.c" to the same ID "a.b.c".
  one <- function(v) matrix(v, 1, 1, dimnames=list("y", "y"))
  covs <- list("a.b"=list(c=one(1.1), d=one(0.9), e=one(1)),
               "a"=list("b.c"=one(1.2), d=one(0.8), f=one(1)))
  means <- list("a.b"=list(c(y=0), c(y=1), c(y=-1)),
                "a"=list(c(y=0.5), c(y=0), c(y=2)))
  n <- list(rep(5, 3), rep(5, 3))
  model <- "group: 1\n level: 1\n y ~~ y\n level: 2\n y ~~ y\ngroup: 2\n level: 1\n y ~~ y\n level: 2\n y ~~ y"
  expect_no_warning(fit <- sem(RAM=lavaan2RAM(model, ngroups=2), Cov=covs,
                               means=means, numObs=n, group="g",
                               cluster="id"))
  ids <- tapply(fit$data$g, fit$data$id, function(x) length(unique(x)))
  expect_length(ids, 6)
  expect_true(all(ids == 1))
  expect_true(all(is.finite(vcov(fit, robust=TRUE))))
})

test_that("robust covariance respects explicit equality constraints", {
  skip_on_cran()
  set.seed(11)
  id <- rep(seq_len(50), each=6)
  d <- data.frame(id=id, x1=rnorm(300), x2=rnorm(300))
  d$y <- 0.2*d$x1 + 0.4*d$x2 + rnorm(50)[id] + rnorm(300)
  model <- "y ~ a*x1 + b*x2
    b == 2*a"
  for (cluster in list("id", NULL)) {
    fit <- sem(RAM=lavaan2RAM(model), data=d, cluster=cluster)
    covariance <- vcov(fit, robust=TRUE)
    expect_equal(covariance["b", "b"], 4*covariance["a", "a"], tolerance=1e-6)
    expect_equal(covariance["a", "b"], 2*covariance["a", "a"], tolerance=1e-6)
  }
  ## Same as substituting the constraint (imxRobustSE() is not a reference
  ## here: it uses OpenMx's Hessian, which is wrong with mxConstraints)
  substituted <- sem(RAM=lavaan2RAM(model), data=d, replace.constraints=TRUE)
  names <- intersect(names(coef(fit)), names(coef(substituted)))
  expect_equal(covariance[names, names], vcov(substituted, robust=TRUE)[names, names],
               tolerance=1e-5)

  ## Two-level, with a constraint spanning the levels
  set.seed(3)
  id <- rep(seq_len(40), each=5)
  d <- data.frame(id=id, x=rnorm(200))
  d$y <- 0.4*d$x + rnorm(40)[id] + rnorm(200)
  fit <- sem(RAM=lavaan2RAM("level: 1
    y ~ bw*x
    y ~~ ew*y
  level: 2
    y ~~ eb*y
    y ~ my*1
    eb == 0.5*ew"), data=d, cluster="id")
  covariance <- vcov(fit, robust=TRUE)
  expect_equal(covariance["eb", "eb"], 0.25*covariance["ew", "ew"],
               tolerance=1e-6)
  expect_equal(covariance["eb", "ew"], 0.5*covariance["ew", "ew"],
               tolerance=1e-6)
})

test_that("a bound does not remove robust SEs of unrelated defined parameters", {
  skip_on_cran()
  set.seed(7)
  id <- rep(seq_len(60), sample(1:5, 60, replace=TRUE))
  v <- runif(length(id), 0.02, 0.1)
  x <- rnorm(length(id))
  d <- data.frame(id=id, v=v, x=x,
                  yi=0.3 + 0.1*x + 0.6*rnorm(length(id), sd=sqrt(v)))
  model <- "f =~ 1*yi
    f ~~ tau2*f
    yi ~~ data.v*yi
    f ~ mu*1
    twomu := 2*mu
    tau2mu := tau2 + mu"
  fit <- sem(RAM=lavaan2RAM(model, obs.variables="yi", std.lv=FALSE),
             data=d, cluster="id", lbound=list(tau2=0))
  s <- expect_robust_warnings(summary(fit, robust=TRUE), fit)
  se.mu <- s$coefficients["mu", "Std.Error"]
  expect_true(is.finite(se.mu))
  expect_equal(unname(s$mxalgebras["twomu", "ubound"] -
                        s$mxalgebras["twomu", "estimate"]),
               qnorm(0.975)*2*se.mu, tolerance=1e-2)
  expect_true(is.na(s$mxalgebras["tau2mu", "ubound"]))

  ## A constraint combined with a bound: the bounded parameter is fixed and
  ## the constraint is still respected.
  constrained <- "f =~ 1*yi
    f ~~ tau2*f
    yi ~~ data.v*yi
    f ~ mu*1
    yi ~ b*x
    b == mu/3"
  fit <- sem(RAM=lavaan2RAM(constrained, obs.variables=c("yi", "x"),
                            std.lv=FALSE),
             data=d, cluster="id", lbound=list(tau2=0))
  covariance <- expect_robust_warnings(vcov(fit, robust=TRUE), fit)
  expect_equal(covariance["b", "b"], covariance["mu", "mu"]/9, tolerance=1e-6)
  fixed <- sem(RAM=lavaan2RAM(sub("tau2\\*f", "0*f", constrained),
                              obs.variables=c("yi", "x"), std.lv=FALSE),
               data=d, cluster="id")
  expect_true(fixed$mx.fit@output$status$code %in% c(0, 1))
  names <- names(coef(fixed))
  expect_equal(coef(fit)[names], coef(fixed), tolerance=1e-4)
  expect_equal(covariance[names, names], vcov(fixed, robust=TRUE),
               tolerance=1e-3)
})

test_that("lavaan2RAM() keeps inequalities on single parameters as bounds", {
  set.seed(11)
  d <- data.frame(x1=rnorm(300), x2=rnorm(300), x3=rnorm(300),
                  g=rep(c("A", "B"), each=150))
  d$y <- 0.2*d$x1 + 0.4*d$x2 + 0.1*d$x3 + rnorm(300)
  agrees <- function(model, ..., lavaan.args=list()) {
    fit <- sem(RAM=lavaan2RAM(model, ...), data=d,
               group=lavaan.args$group)
    reference <- do.call(lavaan::sem,
                         c(list(model, data=d, fixed.x=FALSE,
                                meanstructure=TRUE), lavaan.args))
    expected <- lavaan::parameterEstimates(reference)
    expected <- expected[expected$label != "" & !duplicated(expected$label), ]
    expected <- setNames(expected$est, expected$label)
    names <- intersect(names(coef(fit)), names(expected))
    expect_equal(coef(fit)[names], expected[names], tolerance=1e-4)
    fit
  }
  ## lavaan's own "lower"/"upper" columns, with one or both present
  ram <- lavaan2RAM("y ~ a*x1 + b*x2 + c*x3\n b > 3\n c < 0.05")
  expect_identical(attr(ram, "lbound"), list(b=3))
  expect_identical(attr(ram, "ubound"), list(c=0.05))
  fit <- agrees("y ~ a*x1 + b*x2 + c*x3\n b > 3\n c < 0.05")
  expect_equal(unname(coef(fit)[c("b", "c")]), c(3, 0.05), tolerance=1e-6)
  ## A constant on the left or a constant expression
  fit <- agrees("y ~ a*x1 + b*x2 + c*x3\n 2 < c\n a > -1 + 1.3")
  expect_equal(unname(coef(fit)[c("a", "c")]), c(0.3, 2), tolerance=1e-6)
  ## An unlabeled parameter bounded with lower()
  ram <- lavaan2RAM("y ~ a*x1 + lower(0.6)*x2")
  expect_identical(attr(ram, "lbound"), list(yONx2=0.6))
  ## Multiple groups, where lavaan keeps the inequality as a constraint row,
  ## including a label shared across groups
  fit <- agrees("y ~ c(a1,a2)*x1 + c(b1,b2)*x2\n b1 > 3", ngroups=2,
                lavaan.args=list(group="g"))
  expect_equal(unname(coef(fit)["b1"]), 3, tolerance=1e-6)
  fit <- agrees("y ~ c(a1,a2)*x1 + c(b,b)*x2\n b > 3", ngroups=2,
                lavaan.args=list(group="g"))
  expect_equal(unname(coef(fit)["b"]), 3, tolerance=1e-6)
  ## Bounds given to sem() take precedence
  fit <- sem(RAM=lavaan2RAM("y ~ a*x1 + b*x2\n b > 3"), data=d,
             lbound=list(b=1))
  expect_equal(unname(coef(fit)["b"]), 1, tolerance=1e-6)
  ## An inequality between parameters remains a constraint
  ram <- lavaan2RAM("y ~ a*x1 + b*x2\n b > 3*a")
  expect_null(attr(ram, "lbound"))
  expect_length(ram$mxalgebras, 1L)
})

test_that("bounds from lavaan syntax reach two-level, wls(), and osmasem() fits", {
  set.seed(5)
  id <- rep(seq_len(60), sample(3:8, 60, replace=TRUE))
  d <- data.frame(id=id, x=rnorm(length(id)), g=ifelse(id <= 30, "A", "B"))
  d$y <- 0.4*d$x + rnorm(60, sd=0.6)[id] + rnorm(length(id))
  ## Two-level: equal to the same bounds given explicitly. (lavaan itself
  ## does not enforce bounds on level-2 parameters, so it is no reference.)
  model <- "level: 1\n y ~ bw*x\n y ~~ ew*y\nlevel: 2\n y ~~ tb*y\n y ~ my*1"
  bounded <- sem(RAM=lavaan2RAM(paste0(model, "\n tb > 0.8\n bw < 0.2")),
                 data=d, cluster="id")
  explicit <- sem(RAM=lavaan2RAM(model), data=d, cluster="id",
                  lbound=list(tb=0.8), ubound=list(bw=0.2))
  expect_equal(coef(bounded), coef(explicit), tolerance=1e-6)
  expect_equal(unname(coef(bounded)[c("tb", "bw")]), c(0.8, 0.2),
               tolerance=1e-6)
  ## Group x level, with "group: A" labels
  block <- function(k) gsub("#", k, "level: 1\n y ~ bw#*x\n y ~~ ew#*y\nlevel: 2\n y ~~ tb#*y\n y ~ my#*1")
  model <- paste0("group: A\n", block(1), "\ngroup: B\n", block(2))
  ram <- lavaan2RAM(paste0(model, "\n tb2 > 0.9"), ngroups=2)
  expect_identical(attr(ram, "lbound"), list(tb2=0.9))
  bounded <- sem(RAM=ram, data=d, group="g", cluster="id")
  explicit <- sem(RAM=lavaan2RAM(model, ngroups=2), data=d, group="g",
                  cluster="id", lbound=list(tb2=0.9))
  expect_equal(coef(bounded), coef(explicit), tolerance=1e-6)

  ## wls(): an identified correlation structure
  R <- matrix(c(1, .3, .4, .3, 1, .5, .4, .5, 1), 3,
              dimnames=list(c("x", "m", "y"), c("x", "m", "y")))
  ram <- lavaan2RAM("m ~ a*x\n y ~ b*m + c*x\n x ~~ 1*x\n c < 0.05",
                    obs.variables=c("x", "m", "y"))
  fit <- wls(Cov=R, aCov=asyCov(R, n=200), n=200, RAM=ram)
  expect_equal(unname(coef(fit)["c"]), 0.05, tolerance=1e-6)

  ## osmasem(): a bound on the A matrix
  ram <- lavaan2RAM("W2 ~ w2w*W1 + s2w*S1
    S2 ~ w2s*W1 + s2s*S1
    W1 ~~ w1WITHs1*S1
    W2 ~~ w2WITHs2*S2
    W1 ~~ 1*W1
    S1 ~~ 1*S1
    W2 ~~ Errw2*W2
    S2 ~~ Errs2*S2
    s2w < 0", obs.variables=c("W1", "S1", "W2", "S2"))
  fit <- osmasem(model.name="bounded", RAM=ram,
                 data=Cor2DataFrame(Nohe15A1, acov="weighted"))
  expect_lte(unname(coef(fit)["s2w"]), 1e-8)
})

test_that("summary.wls() handles df=0 as lavaan and df<0 without an error", {
  R <- matrix(c(1, .3, .4, .3, 1, .5, .4, .5, 1), 3,
              dimnames=list(c("x", "m", "y"), c("x", "m", "y")))
  fit <- function(model)
    wls(Cov=R, aCov=asyCov(R, n=200), n=200,
        RAM=lavaan2RAM(model, obs.variables=c("x", "m", "y")))
  stat <- function(object) summary(object)$stat[, 1]
  m0 <- "m ~ a*x\n y ~ b*m + c*x\n x ~~ 1*x"

  ## Saturated: lavaan gives p=NA, RMSEA=0 with CI (0,0), TLI=1, CFI=1
  s <- stat(fit(m0))
  expect_equal(unname(s["DF of target model"]), 0)
  expect_true(is.na(s["p value of target model"]))
  expect_equal(unname(s[c("RMSEA", "RMSEA lower 95% CI", "RMSEA upper 95% CI",
                          "TLI", "CFI")]), c(0, 0, 0, 1, 1))

  ## df=0 with an active bound: chi-square > 0, CFI from its formula
  s <- stat(fit(paste(m0, "\n c < 0.05")))
  expect_gt(s["Chi-square of target model"], 1)
  expect_true(is.na(s["p value of target model"]))
  expect_equal(unname(s[c("RMSEA", "TLI")]), c(0, 1))
  tT <- s["Chi-square of target model"]
  tB <- s["Chi-square of independence model"]
  dfB <- s["DF of independence model"]
  expect_equal(unname(s["CFI"]), unname(1 - tT/max(tT, tB-dfB)))

  ## df<0: a warning and NA fit indices instead of an error
  under <- suppressWarnings(fit("m ~ a*x\n y ~ b*m + c*x\n x ~~ vx*x"))
  expect_warning(s <- stat(under), "more free parameters")
  expect_equal(unname(s["DF of target model"]), -1)
  expect_true(all(is.na(s[c("p value of target model", "RMSEA", "RMSEA lower 95% CI",
                            "RMSEA upper 95% CI", "TLI", "CFI")])))

  ## df>0: use total N in the normal-ML RMSEA convention
  s <- stat(fit("m ~ a*x\n y ~ b*m\n x ~~ 1*x"))
  tT <- s["Chi-square of target model"]
  dfT <- s["DF of target model"]
  tB <- s["Chi-square of independence model"]
  dfB <- s["DF of independence model"]
  expect_equal(unname(s["p value of target model"]), pchisq(unname(tT), dfT, lower.tail=FALSE))
  expect_equal(unname(s["RMSEA"]), unname(sqrt(max((tT-dfT)/200, 0)/dfT)))
  expect_equal(unname(s["TLI"]), unname((tB/dfB - tT/dfT)/(tB/dfB - 1)))
  expect_equal(unname(s["CFI"]), unname(1 - max(tT-dfT, 0)/max(tT-dfT, tB-dfB, 0)))
})

test_that("summary.tssem1FEM() handles df=0 as lavaan", {
  R <- matrix(c(1, .3, .4, .3, 1, .5, .4, .5, 1), 3,
              dimnames=list(c("x", "m", "y"), c("x", "m", "y")))
  s <- summary(tssem1(list(R), 200, method="FEM"))$stat[, 1]
  expect_equal(unname(s["DF of target model"]), 0)
  expect_true(is.na(s["p value of target model"]))
  expect_equal(unname(s[c("RMSEA", "RMSEA lower 95% CI", "RMSEA upper 95% CI", "TLI")]),
               c(0, 0, 0, 1))
})

test_that("a defined parameter of several bounded parameters has an NA robust SE", {
  skip_on_cran()
  set.seed(1)
  d <- data.frame(x1=rnorm(300), x2=rnorm(300))
  d$y <- 1.5*d$x1 + 1.5*d$x2 + rnorm(300)
  fit <- sem(RAM=lavaan2RAM("y ~ a*x1 + b*x2 + c*1
    difference := a - b
    total := a + b
    twoc := 2*c", obs.variables=c("y", "x1", "x2")),
    data=d, ubound=list(a=1, b=1))
  expect_equal(unname(coef(fit)[c("a", "b")]), c(1, 1), tolerance=1e-6)
  s <- expect_robust_warnings(summary(fit, robust=TRUE), fit)
  ## The changes in a and b cancel in a - b if they are shifted together
  expect_true(is.na(s$mxalgebras["difference", "ubound"]))
  expect_true(is.na(s$mxalgebras["total", "ubound"]))
  expect_equal(unname(s$mxalgebras["twoc", "ubound"] - s$mxalgebras["twoc", "estimate"]),
               qnorm(0.975)*2*s$coefficients["c", "Std.Error"], tolerance=1e-2)
})

test_that("constrained covariances agree with lavaan and with substituted constraints", {
  skip_on_cran()
  set.seed(1)
  d <- data.frame(x=rnorm(200))
  d$m <- 0.5*d$x + rnorm(200)
  d$y <- 0.4*d$m + 0.1*d$x + rnorm(200)
  base <- "m ~ a*x\n y ~ b*m + c*x\n ind := a*b\n "
  ## OpenMx's own SEs with these constraints are wrong (e.g., 0.0418 for a
  ## with b == 2*a, against 0.0345), or NA from a negative variance.
  for (constraint in c("b == 2*a", "c < -0.5*a", "c > 2*a", "c > -5*a")) {
    model <- paste0(base, constraint)
    fit <- sem(RAM=lavaan2RAM(model, obs.variables=c("x", "m", "y")), data=d)
    s <- summary(fit)
    reference <- suppressWarnings(lavaan::parameterEstimates(
      lavaan::sem(model, data=d, fixed.x=FALSE, meanstructure=TRUE,
                  information="observed")))
    reference <- setNames(reference$se, ifelse(reference$label == "",
                                               reference$lhs, reference$label))
    expect_equal(s$coefficients[c("a", "b", "c"), "Std.Error"],
                 unname(reference[c("a", "b", "c")]), tolerance=1e-4)
    expect_equal(unname(s$mxalgebras["ind", "ubound"] - s$mxalgebras["ind", "estimate"]),
                 qnorm(0.975)*unname(reference["ind"]), tolerance=1e-4)
    expect_equal(unname(sqrt(diag(vcov(fit)))[c("a", "b", "c")]),
                 s$coefficients[c("a", "b", "c"), "Std.Error"])
  }

  ## A nonlinear equality: lavaan omits the curvature of the constraint, so
  ## the reference is the fit with the constraint substituted.
  model <- paste0(base, "b == a^2 + 0.42")
  fit <- sem(RAM=lavaan2RAM(model, obs.variables=c("x", "m", "y")), data=d)
  substituted <- sem(RAM=lavaan2RAM(model, obs.variables=c("x", "m", "y")),
                     data=d, replace.constraints=TRUE)
  s <- summary(fit)
  s.sub <- summary(substituted)
  expect_equal(s$coefficients[c("a", "c"), "Std.Error"],
               s.sub$coefficients[c("a", "c"), "Std.Error"], tolerance=1e-4)
  expect_equal(unname(s$mxalgebras["ind", "ubound"] - s$mxalgebras["ind", "estimate"]),
               unname(s.sub$mxalgebras["ind", "ubound"] -
                        s.sub$mxalgebras["ind", "estimate"]), tolerance=1e-4)

  ## Binding inequality constraints are found from the estimates; OpenMx
  ## reports a zero Jacobian row for them
  fit <- sem(RAM=lavaan2RAM(paste0(base, "c > 2*a\n c > -5*a"),
                            obs.variables=c("x", "m", "y")), data=d)
  p <- names(coef(fit$mx.fit))
  jacobian <- .active.constraint.jacobian(fit$mx.fit, p)
  expect_equal(nrow(jacobian), 1)
  expect_equal(unname(jacobian[1, c("a", "b", "c")]), c(-2, 0, 1), tolerance=1e-6)

  ## A bound together with an equality and a binding inequality: the bounded
  ## parameter is treated as fixed, as in a refit with it fixed
  model <- "m ~ a*x\n y ~ b*m + c*x\n b == 2*a\n c > 2*a"
  ## The bounded fit ends with a status code other than 0 or 1, which sem()
  ## warns about; its estimates are still close to those of the refit.
  fit <- suppressWarnings(sem(RAM=lavaan2RAM(model, obs.variables=c("x", "m", "y")),
                              data=d, ubound=list(mWITHm=0.5)))
  mx <- fit$mx.fit
  status <- .binding.status(mx, names(coef(mx)))
  expect_true(status$bounded[["mWITHm"]])
  bread <- .robust.bread(mx, names(coef(mx)), status)
  fixed <- sem(RAM=lavaan2RAM(paste0(model, "\n m ~~ 0.5*m"),
                              obs.variables=c("x", "m", "y")), data=d)
  names <- names(coef(fixed))
  ## About 1e-4 from the refit
  expect_equal(coef(mx)[names], coef(fixed), tolerance=1e-3)
  expect_equal(bread[names, names], vcov(fixed), tolerance=1e-3)
  expect_true(all(is.na(bread["mWITHm", ])))
})

test_that("robust covariances warn when OpenMx reports a problem", {
  set.seed(1)
  d <- data.frame(x=rnorm(100))
  d$y <- 0.3*d$x + rnorm(100)
  fit <- sem(RAM=lavaan2RAM("y ~ b*x", obs.variables=c("y", "x")), data=d)
  expect_silent(vcov(fit, robust=TRUE))
  fit$mx.fit@output$status$code <- 6L
  expect_warning(vcov(fit, robust=TRUE), "OpenMx status code is 6")
})

test_that("lavaan2RAM() keeps a negative right-hand side of '<' constraints", {
  formula <- function(model, ...) {
    ram <- lavaan2RAM(model, obs.variables=c("x", "m", "y"), ...)
    ## Multiple-group constraints are kept with the first group
    if (inherits(ram, "RAM_multigroup")) ram <- ram[[1]]
    deparse(ram$mxalgebras$constraint1@formula)
  }
  base <- "m ~ a*x\n y ~ b*m + c*x\n "
  ## "c<-0.5*a" would be parsed as an assignment
  expect_identical(formula(paste0(base, "c < -0.5*a")), "c < -0.5 * a")
  expect_false(exists("c", envir=environment(), inherits=FALSE))
  expect_identical(formula(paste0(base, "c > -2*a")), "c > -2 * a")
  expect_identical(formula(paste0(base, "b == -a")), "b == -a")
  groups <- "group: 1\n m ~ a1*x\n y ~ b1*m + c1*x\n group: 2\n m ~ a2*x\n y ~ b2*m + c2*x\n "
  expect_identical(formula(paste0(groups, "c2 < -0.5*a1"), ngroups=2), "c2 < -0.5 * a1")

})

test_that("active inequality constraints are satisfied and agree with lavaan", {
  skip_on_cran()
  set.seed(5)
  d <- data.frame(x=rnorm(400))
  d$m <- 0.5*d$x + rnorm(400)
  d$y <- 0.4*d$m + 0.1*d$x + rnorm(400)
  base <- "m ~ a*x\n y ~ b*m + c*x\n "
  ## Under OpenMx's default feasibility tolerance (0.05), these constraints
  ## were violated by 0.047 and 0.021, with different estimates.
  for (constraint in c("c < -0.5*a", "c > 2*a", "b == 2*a")) {
    model <- paste0(base, constraint)
    fit <- sem(RAM=lavaan2RAM(model, obs.variables=c("x", "m", "y")), data=d)
    expect_equal(fit$mx.fit@options[["Feasibility tolerance"]], 1e-6)
    expect_lte(max(abs(fit$mx.fit@output$constraintFunctionValues)), 1e-5)
    reference <- lavaan::parameterEstimates(lavaan::sem(model, data=d, fixed.x=FALSE,
                                                        meanstructure=TRUE))
    reference <- setNames(reference$est, reference$label)[c("a", "b", "c")]
    expect_equal(coef(fit)[c("a", "b", "c")], reference, tolerance=1e-4)
  }

  ## Models without constraints keep OpenMx's settings
  fit <- sem(RAM=lavaan2RAM("m ~ a*x\n y ~ b*m", obs.variables=c("x", "m", "y")), data=d)
  expect_null(fit$mx.fit@options[["Feasibility tolerance"]])

  ## A global tolerance set by the user is respected
  old <- mxOption(NULL, "Feasibility tolerance")
  on.exit(mxOption(NULL, "Feasibility tolerance", old), add=TRUE)
  mxOption(NULL, "Feasibility tolerance", 0.01)
  fit <- sem(RAM=lavaan2RAM(paste0(base, "c < -0.5*a"), obs.variables=c("x", "m", "y")),
             data=d)
  expect_null(fit$mx.fit@options[["Feasibility tolerance"]])
  mxOption(NULL, "Feasibility tolerance", old)

  ## A remaining violation gives a warning
  fit <- sem(RAM=lavaan2RAM(paste0(base, "c < -0.5*a"), obs.variables=c("x", "m", "y")),
             data=d)$mx.fit
  fit@output$constraintFunctionValues[] <- 1e-3
  expect_warning(.check.constraints(fit, rerun=FALSE), "constraint is violated by 0.001")
})

test_that("wls() and osmasem() satisfy active inequality constraints", {
  skip_on_cran()
  R <- matrix(c(1, .3, .4, .3, 1, .5, .4, .5, 1), 3,
              dimnames=list(c("x", "m", "y"), c("x", "m", "y")))
  ram <- lavaan2RAM("m ~ a*x\n y ~ b*m + c*x\n x ~~ 1*x\n c < -0.5*a",
                    obs.variables=c("x", "m", "y"))
  fit <- wls(Cov=R, aCov=asyCov(R, n=200), n=200, RAM=ram)
  expect_equal(fit$mx.fit@options[["Feasibility tolerance"]], 1e-6)
  expect_lte(unname(coef(fit)["c"] + 0.5*coef(fit)["a"]), 1e-5)

  ram <- lavaan2RAM("W2 ~ w2w*W1 + s2w*S1
    S2 ~ w2s*W1 + s2s*S1
    W1 ~~ w1WITHs1*S1
    W2 ~~ w2WITHs2*S2
    W1 ~~ 1*W1
    S1 ~~ 1*S1
    W2 ~~ Errw2*W2
    S2 ~~ Errs2*S2
    s2w < -0.5*w2s", obs.variables=c("W1", "S1", "W2", "S2"))
  data <- Cor2DataFrame(Nohe15A1, acov="weighted")
  for (fit in list(osmasem(model.name="ineq", RAM=ram, data=data),
                   osmasem2(model.name="ineq2", RAM=ram, data=data))) {
    expect_equal(fit$mx.fit@options[["Feasibility tolerance"]], 1e-6)
    coefs <- coef(fit$mx.fit)
    expect_lte(unname(coefs["s2w"] + 0.5*coefs["w2s"]), 1e-5)
  }
})

test_that("osmasem() and osmasem2() SEs with an equality constraint match equal labels", {
  skip_on_cran()
  base <- "W2 ~ w2w*W1 + s2w*S1
    S2 ~ w2s*W1 + s2s*S1
    W1 ~~ w1WITHs1*S1
    W2 ~~ w2WITHs2*S2
    W1 ~~ 1*W1
    S1 ~~ 1*S1
    W2 ~~ Errw2*W2
    S2 ~~ Errs2*S2"
  variables <- c("W1", "S1", "W2", "S2")
  constrained <- lavaan2RAM(paste(base, "\n s2w == w2s"), obs.variables=variables)
  labeled <- lavaan2RAM(sub("w2s\\*W1", "s2w*W1", base), obs.variables=variables)
  data <- Cor2DataFrame(Nohe15A1, acov="weighted")
  for (f in list(osmasem, osmasem2)) {
    a <- f(model.name="constrained", RAM=constrained, data=data)
    b <- f(model.name="labeled", RAM=labeled, data=data)
    ## OpenMx's own SEs of the constrained osmasem() fit differ by about 3%
    for (robust in c(FALSE, TRUE)) {
      sa <- summary(a, robust=robust)$parameters
      sb <- summary(b, robust=robust)$parameters
      names <- intersect(sa$name, sb$name)
      expect_equal(sa$Std.Error[match(names, sa$name)],
                   sb$Std.Error[match(names, sb$name)], tolerance=1e-4)
    }
    names <- intersect(rownames(vcov(a)), rownames(vcov(b)))
    expect_equal(vcov(a)[names, names], vcov(b)[names, names], tolerance=1e-4)
  }
})

test_that("labeled group blocks are numbered in syntax order, as in lavaan", {
  set.seed(1)
  d <- data.frame(g=rep(c("A", "B"), each=200), x=rnorm(400))
  d$y <- ifelse(d$g == "A", 0.8, -0.3)*d$x + rnorm(400)
  d$gn <- ifelse(d$g == "A", 1, 2)
  agrees <- function(model, group) {
    ram <- lavaan2RAM(model, obs.variables=c("y", "x"), ngroups=2)
    fit <- sem(RAM=ram, data=d, group=group)
    reference <- lavaan::parameterEstimates(lavaan::sem(model, data=d, group=group,
                                                        fixed.x=FALSE, meanstructure=TRUE))
    reference <- reference[reference$label != "" & !duplicated(reference$label), ]
    reference <- setNames(reference$est, reference$label)
    names <- intersect(names(coef(fit)), names(reference))
    expect_equal(coef(fit)[names], reference[names], tolerance=1e-5)
    list(ram=ram, fit=fit)
  }
  ## Labels in data order: no warning; a cross-group constraint is kept
  expect_no_warning(out <- agrees("group: A\n y ~ bA*x\n group: B\n y ~ bB*x\n bA == -2*bB", "g"))
  expect_identical(attr(out$ram, "group.labels"), c("A", "B"))
  expect_false(is.null(out$ram[[1]]$mxalgebras$constraint1))
  expect_equal(unname(coef(out$fit)["bA"]), -2*unname(coef(out$fit)["bB"]), tolerance=1e-6)

  ## Labels in the other order are still paired by position, with a warning
  expect_warning(out <- agrees("group: B\n y ~ bB*x\n group: A\n y ~ bA*x\n bB > -0.1", "g"),
                 "model group 'B' is fitted to data group 'A'")
  expect_identical(attr(out$ram, "lbound"), list(bB=-0.1))
  expect_true(grepl("bB", out$ram[[1]]$A["y", "x"]))
  expect_warning(agrees("group: 2\n y ~ b2*x\n group: 1\n y ~ b1*x", "gn"),
                 "model group '2' is fitted to data group '1'")

  ## lavaan's c() syntax has no block labels
  expect_no_warning(out <- agrees("y ~ c(b1, b2)*x", "gn"))
  expect_null(attr(out$ram, "group.labels"))

  ## Combined group and level blocks, with a cross-group constraint
  set.seed(2)
  id <- rep(seq_len(60), each=6)
  d <- data.frame(id=id, g=ifelse(id <= 30, "A", "B"), x=rnorm(360))
  d$y <- ifelse(d$g == "A", 0.6, 0.2)*d$x +
    rnorm(60, sd=rep(c(0.8, 0.4), each=30))[id] + rnorm(360)
  block <- function(L) sprintf("group: %s\n level: 1\n y ~ bw%s*x\n y ~~ ew%s*y\n level: 2\n y ~~ eb%s*y\n y ~ m%s*1\n",
                               L, L, L, L, L)
  model <- paste0(block("B"), block("A"), " ebA == 4*ebB\n")
  ram <- lavaan2RAM(model, ngroups=2)
  expect_identical(attr(ram, "group.labels"), c("B", "A"))
  expect_warning(fit <- sem(RAM=ram, data=d, group="g", cluster="id"),
                 "model group 'B' is fitted to data group 'A'")
  expect_equal(unname(coef(fit)["ebA"]), 4*unname(coef(fit)["ebB"]), tolerance=1e-5)
  reference <- lavaan::parameterEstimates(lavaan::sem(model, data=d, group="g", cluster="id"))
  reference <- reference[reference$label != "" & !duplicated(reference$label), ]
  reference <- setNames(reference$est, reference$label)
  names <- intersect(names(coef(fit)), names(reference))
  expect_equal(coef(fit)[names], reference[names], tolerance=1e-3)
})

test_that("osmasem() applies lavaan bounds outside the A matrix", {
  base <- "W2 ~ w2w*W1 + s2w*S1
    S2 ~ w2s*W1 + s2s*S1
    W1 ~~ w1WITHs1*S1
    W2 ~~ w2WITHs2*S2
    W1 ~~ 1*W1
    S1 ~~ 1*S1
    W2 ~~ Errw2*W2
    S2 ~~ Errs2*S2"
  variables <- c("W1", "S1", "W2", "S2")
  data <- Cor2DataFrame(Nohe15A1, acov="weighted")
  ## Unbounded estimates: w1WITHs1 = 0.38, w2WITHs2 = 0.17. Errw2 is computed
  ## by osmasem(), so its bound is ignored with a warning.
  expect_warning(fit <- osmasem(RAM=lavaan2RAM(paste(base, "\n w2WITHs2 > 0.3
    w1WITHs1 < 0.2\n Errw2 > 0.5"), obs.variables=variables), data=data),
    "Errw2 are ignored")
  expect_equal(unname(coef(fit$mx.fit)[c("w1WITHs1", "w2WITHs2")]), c(0.2, 0.3),
               tolerance=1e-6)
  fixed <- osmasem(RAM=lavaan2RAM(sub("w1WITHs1\\*S1", "0.2*S1",
                                      sub("w2WITHs2\\*S2", "0.3*S2", base)),
                                  obs.variables=variables), data=data)
  expect_equal(fit$mx.fit@output$fit, fixed$mx.fit@output$fit, tolerance=1e-6)
  ## A.lbound takes precedence over a bound on the A matrix in RAM
  fit <- osmasem(RAM=lavaan2RAM(paste(base, "\n s2w > 0.5"), obs.variables=variables),
                 data=data, A.lbound=matrix(-1, 4, 4))
  expect_lt(unname(coef(fit$mx.fit)["s2w"]), 0.5)
})

test_that("internal mxTryHard() fallbacks run quietly", {
  ram <- lavaan2RAM("W2 ~ w2w*W1 + s2w*S1
    S2 ~ w2s*W1 + s2s*S1
    W1 ~~ w1WITHs1*S1
    W2 ~~ w2WITHs2*S2
    W1 ~~ 1*W1
    S1 ~~ 1*S1
    W2 ~~ Errw2*W2
    S2 ~~ Errs2*S2
    s2w < -0.5*w2s", obs.variables=c("W1", "S1", "W2", "S2"))
  ## This fit ends with a status code other than 0 or 1, so osmasem2() falls
  ## back to mxTryHard(), which printed its progress
  expect_silent(fit <- osmasem2(model.name="quiet", RAM=ram,
                                data=Cor2DataFrame(Nohe15A1, acov="weighted")))
  expect_lte(unname(coef(fit$mx.fit)["s2w"] + 0.5*coef(fit$mx.fit)["w2s"]), 1e-5)
})

test_that("sem(run=FALSE) returns models with the tighter feasibility tolerance", {
  set.seed(5)
  d <- data.frame(x=rnorm(400), g=rep(1:2, each=200), id=rep(1:80, each=5))
  d$m <- 0.5*d$x + rnorm(400)
  d$y <- 0.4*d$m + 0.1*d$x + rnorm(80)[d$id] + rnorm(400)
  ## The single-level model violated c < -0.5*a by 0.047 when the returned
  ## model was run with OpenMx's default tolerance
  models <- list(
    single=list(RAM=lavaan2RAM("m ~ a*x\n y ~ b*m + c*x\n c < -0.5*a",
                               obs.variables=c("x", "m", "y")), data=d),
    multigroup=list(RAM=lavaan2RAM("group: 1\n m ~ a*x\n y ~ b1*m + c*x
      group: 2\n m ~ a*x\n y ~ b2*m + c*x\n c < -0.5*a",
                                   obs.variables=c("x", "m", "y"), ngroups=2),
                    data=d, group="g"),
    twolevel=list(RAM=lavaan2RAM("level: 1\n y ~ c*x\n level: 2\n y ~~ eb*y\n c < -0.5*eb"),
                  data=d, cluster="id"),
    combined=list(RAM=lavaan2RAM("group: 1\n level: 1\n y ~ c*x\n level: 2\n y ~~ eb1*y
      group: 2\n level: 1\n y ~ c*x\n level: 2\n y ~~ eb2*y\n c < -0.5*eb1",
                                 ngroups=2),
                  data=d, group="g", cluster="id"))
  for (args in models) {
    model <- do.call(sem, c(args, run=FALSE))$mx.fit
    expect_equal(model@options[["Feasibility tolerance"]], 1e-6)
    fit <- mxRun(model, silent=TRUE, suppressWarnings=TRUE)
    expect_lte(max(abs(fit@output$constraintFunctionValues)), 1e-5)
  }
})

test_that("nonlinear constraints use the Hessian of the Lagrangian", {
  skip_on_cran()
  ## Explicit constraints and their substitution into the model are the same
  ## model, so their model-based and robust covariances must agree. Without
  ## the curvature of a nonlinear constraint they differed by about 10%.
  same.se <- function(explicit, substituted, par, tolerance) {
    expect_equal(coef(explicit)[par], coef(substituted)[par], tolerance=tolerance)
    expect_equal(sqrt(diag(vcov(explicit)))[par],
                 sqrt(diag(vcov(substituted)))[par], tolerance=tolerance)
    expect_equal(sqrt(diag(vcov(explicit, robust=TRUE)))[par],
                 sqrt(diag(vcov(substituted, robust=TRUE)))[par],
                 tolerance=tolerance)
  }
  set.seed(5)
  n <- 300
  d <- data.frame(id=rep(seq_len(60), each=5), x1=rnorm(n), x2=rnorm(n))
  d$y <- .5*d$x1 + .4*d$x2 + rnorm(n) + rep(rnorm(60), each=5)
  obs <- c("y", "x1", "x2")
  fit <- function(model, ...)
    sem(RAM=lavaan2RAM(model, obs.variables=obs), data=d, ...)

  ## Nonlinear equality (false in the population), clustered.
  model <- "y ~ a*x1 + b*x2\n b == a^2 + 0.1"
  same.se(fit(model, cluster="id"),
          fit(model, cluster="id", replace.constraints=TRUE), "a", 1e-4)

  ## Linear equality: no curvature.
  model <- "y ~ a*x1 + b*x2\n b == 2*a"
  same.se(fit(model, cluster="id"),
          fit(model, cluster="id", replace.constraints=TRUE), "a", 1e-4)

  ## Binding nonlinear inequality, conditional on the active constraint.
  explicit <- fit("y ~ a*x1 + b*x2\n b < a^2 - 0.1")
  expect_equal(unname(coef(explicit)["b"]),
               unname(coef(explicit)["a"]^2 - 0.1), tolerance=1e-5)
  same.se(explicit, fit("y ~ a*x1 + b*x2\n b == a^2 - 0.1",
                        replace.constraints=TRUE), "a", 1e-3)

  ## Two-level model with a nonlinear cross-level constraint.
  d$x1 <- d$x1 + rep(rnorm(60), each=5)
  model <- "level: 1
  y ~ bw*x1
  y ~~ ew*y
  x1 ~~ vw*x1
level: 2
  y ~ bb*x1
  y ~~ eb*y
  x1 ~~ vb*x1
  y ~ my*1
  x1 ~ mx*1
  bb == bw^2"
  same.se(sem(RAM=lavaan2RAM(model), data=d, cluster="id"),
          sem(RAM=lavaan2RAM(model), data=d, cluster="id",
              replace.constraints=TRUE), "bw", 1e-4)
})

test_that("tssem2() SEs do not depend on diag.constraints", {
  first <- tssem1(Digman97$data, Digman97$n, method="REM", RE.type="Diag")
  ram <- lavaan2RAM("Alpha =~ a1*A + a2*C + a3*ES
  Beta =~ b1*E + b2*I
  Alpha ~~ r*Beta", obs.variables=c("A", "C", "ES", "E", "I"))
  p <- c("a1", "a2", "a3", "b1", "b2", "r")
  with.constraints <- tssem2(first, RAM=ram, diag.constraints=TRUE,
                             intervals.type="z")
  without <- tssem2(first, RAM=ram, diag.constraints=FALSE,
                    intervals.type="z")
  expect_equal(sqrt(diag(vcov(with.constraints)))[p],
               sqrt(diag(vcov(without)))[p], tolerance=1e-5)
})

test_that("summary.mxsem(t.df=) uses a t reference distribution", {
  set.seed(47)
  id <- rep(seq_len(12), each=4)
  d <- data.frame(id=id, x=rnorm(48))
  d$m <- 0.4*d$x + rnorm(48)
  d$y <- 0.3*d$m + rep(rnorm(12), each=4) + rnorm(48)
  ram <- lavaan2RAM("m ~ a*x\n y ~ b*m\n ind := a*b",
                    obs.variables=c("y", "m", "x"))
  fit <- sem(RAM=ram, data=d, cluster="id")
  z <- summary(fit, robust=TRUE)
  t <- summary(fit, robust=TRUE, t.df=10)
  expect_true(all(c("z value", "Pr(>|z|)") %in% names(z$coefficients)))
  expect_true(all(c("t value", "Pr(>|t|)") %in% names(t$coefficients)))
  expect_null(z$t.df)
  expect_equal(t$t.df, 10)
  ## Standard errors are unchanged; intervals and p values use t(10).
  expect_equal(t$coefficients$Std.Error, z$coefficients$Std.Error)
  expect_equal(t$coefficients$ubound - t$coefficients$Estimate,
               qt(.975, 10)*t$coefficients$Std.Error)
  expect_equal(t$coefficients$"Pr(>|t|)",
               2*pt(-abs(t$coefficients$Estimate/t$coefficients$Std.Error), 10))
  ## Defined parameters use the same reference distribution.
  expect_equal(unname((t$mxalgebras[, "ubound"] - t$mxalgebras[, "estimate"]) /
                        (z$mxalgebras[, "ubound"] - z$mxalgebras[, "estimate"])),
               qt(.975, 10)/qnorm(.975), tolerance=1e-6)
  ## Model-based standard errors can also be used.
  expect_equal(summary(fit, t.df=10)$coefficients$Std.Error,
               summary(fit)$coefficients$Std.Error)
  expect_output(print(t), "t distribution with 10 df")
  expect_error(summary(fit, t.df=0), "single positive, finite number")
  expect_error(summary(fit, t.df=c(5, 6)), "single positive, finite number")
  expect_error(summary(fit, t.df=Inf), "single positive, finite number")
  lb <- sem(RAM=ram, data=d, intervals.type="LB")
  expect_error(summary(lb, t.df=10), "intervals.type='z'")
})

test_that("summary.mxsem() reports the number of independent clusters", {
  set.seed(48)
  d <- data.frame(id=rep(seq_len(15), each=4), x=rnorm(60))
  d$y <- 0.3*d$x + rep(rnorm(15), each=4) + rnorm(60)
  ram <- lavaan2RAM("y ~ b*x", obs.variables=c("y", "x"))
  s <- summary(sem(RAM=ram, data=d, cluster="id"))
  expect_equal(s$no.clusters, 15)
  expect_equal(s$no.studies, 60)
  expect_output(print(s), "Number of clusters: 15")
  ## No cluster argument: nothing is reported.
  s <- summary(sem(RAM=ram, data=d))
  expect_null(s$no.clusters)
  expect_false(any(grepl("Number of clusters", capture.output(print(s)))))

  ## Flat multiple groups: an ID shared by two groups is one unit, as in
  ## the cluster-robust covariance.
  d$g <- rep(c("A", "B"), 30)
  fit <- sem(RAM=lavaan2RAM("y ~ c(b1,b2)*x", ngroups=2), data=d, group="g",
             cluster="id")
  expect_equal(summary(fit)$no.clusters, 15)

  ## Two levels, and two levels within two groups (clusters nested in
  ## groups, summed over groups).
  model <- "level: 1
  y ~ bw*x
  level: 2
  y ~~ eb*y"
  fit <- sem(RAM=lavaan2RAM(model), data=d, cluster="id")
  expect_equal(summary(fit)$no.clusters, 15)
  d$g <- ifelse(d$id <= 7, "A", "B")
  grouped <- "group: A
  level: 1
  y ~ bw*x
  level: 2
  y ~~ eb*y
  group: B
  level: 1
  y ~ bw*x
  level: 2
  y ~~ eb*y"
  fit <- sem(RAM=lavaan2RAM(grouped, ngroups=2), data=d, group="g",
             cluster="id")
  expect_equal(summary(fit)$no.clusters, 15)
})

test_that("replace.constraints=TRUE keeps further equalities on a substituted label", {
  set.seed(4)
  n <- 300
  d <- data.frame(x=rnorm(n))
  d$m <- 0.5*d$x + rnorm(n)
  d$y <- 0.6*d$m + 0.3*d$x + rnorm(n)
  ov <- c("y", "m", "x")
  compare <- function(model, par) {
    explicit <- sem(RAM=lavaan2RAM(model, obs.variables=ov), data=d)
    substituted <- sem(RAM=lavaan2RAM(model, obs.variables=ov), data=d,
                       replace.constraints=TRUE)
    expect_equal(coef(substituted)[par], coef(explicit)[par], tolerance=1e-6)
    expect_equal(substituted$mx.fit@output$Minus2LogLikelihood,
                 explicit$mx.fit@output$Minus2LogLikelihood, tolerance=1e-8)
    substituted
  }
  ## b == c + 0.3 is kept as 2*a == c + 0.3
  fit <- compare("m ~ a*x\n y ~ b*m + c*x\n b == 2*a\n b == c + 0.3", c("a", "c"))
  expect_equal(2*coef(fit)[["a"]], coef(fit)[["c"]] + 0.3, tolerance=1e-6)
  ## Together with a chain
  compare("m ~ a*x\n y ~ b*m + c*x\n b == 2*a\n a == c - 0.1\n b == c + 0.3", "c")
  ## An exact repetition always holds and is dropped
  fit <- compare("m ~ a*x\n y ~ b*m + c*x\n b == 2*a\n b == 2*a", c("a", "c"))
  expect_length(fit$mx.fit@constraints, 0)

  ## Multiple groups
  d2 <- rbind(cbind(d, g="A"), cbind(d[sample(n), ], g="B"))
  model <- "m ~ c(a1,a2)*x\n y ~ c(b1,b2)*m + c(c1,c2)*x
  b2 == 2*a1\n b2 == c1 + 0.5"
  explicit <- sem(RAM=lavaan2RAM(model, ngroups=2), data=d2, group="g")
  substituted <- sem(RAM=lavaan2RAM(model, ngroups=2), data=d2, group="g",
                     replace.constraints=TRUE)
  expect_equal(substituted$mx.fit@output$Minus2LogLikelihood,
               explicit$mx.fit@output$Minus2LogLikelihood, tolerance=1e-6)
  expect_equal(2*coef(substituted)[["a1"]], coef(substituted)[["c1"]] + 0.5,
               tolerance=1e-6)

  ## Two levels, and two levels within two groups
  set.seed(5)
  G <- 80
  dt <- data.frame(id=rep(seq_len(G), each=5), x=rnorm(5*G))
  dt$y <- 0.4*dt$x + rnorm(G, sd=0.6)[dt$id] + rnorm(5*G)
  model <- "level: 1\n y ~ bw*x\n y ~~ ew*y\n level: 2\n y ~~ eb*y\n x ~~ vxb*x
  eb == 0.5*ew\n eb == vxb"
  explicit <- sem(RAM=lavaan2RAM(model), data=dt, cluster="id")
  substituted <- sem(RAM=lavaan2RAM(model), data=dt, cluster="id",
                     replace.constraints=TRUE)
  expect_equal(coef(substituted)[c("bw", "ew")], coef(explicit)[c("bw", "ew")],
               tolerance=1e-6)
  dt$g <- ifelse(dt$id <= 40, "A", "B")
  model <- "group: A\n level: 1\n y ~ bwA*x\n y ~~ ewA*y\n level: 2\n y ~~ ebA*y
  group: B\n level: 1\n y ~ bwB*x\n y ~~ ewB*y\n level: 2\n y ~~ ebB*y
  ebB == 0.5*ewA\n ebB == ebA"
  explicit <- sem(RAM=lavaan2RAM(model, ngroups=2), data=dt, group="g",
                  cluster="id")
  substituted <- sem(RAM=lavaan2RAM(model, ngroups=2), data=dt, group="g",
                     cluster="id", replace.constraints=TRUE)
  par <- c("bwA", "ewA", "bwB", "ewB")
  expect_equal(coef(substituted)[par], coef(explicit)[par], tolerance=1e-6)
})

test_that("lavaan2RAM() accepts any valid name for a defined parameter", {
  set.seed(1)
  d <- data.frame(x=rnorm(100))
  d$m <- 0.5*d$x + rnorm(100)
  d$y <- 0.5*d$m + rnorm(100)
  ov <- c("y", "m", "x")
  reference <- summary(sem(RAM=lavaan2RAM("y ~ a*m\n m ~ b*x\n ind := a*b",
                                          obs.variables=ov), data=d))$mxalgebras
  ## Names of local variables in the former implementation, and a variable
  ## name, which lavaan also accepts
  for (name in c("y", "out", "mxalgebra", "x")) {
    model <- paste0("y ~ a*m\n m ~ b*x\n ", name, " := a*b\n a > b - 3")
    ram <- lavaan2RAM(model, obs.variables=ov)
    expect_named(ram$mxalgebras, c(name, "constraint1"))
    s <- summary(sem(RAM=ram, data=d))
    expect_equal(unname(s$mxalgebras), unname(reference), tolerance=1e-6)
  }
  ## Names that clash with parameters or model entities
  for (name in c("a", "constraint1", "Smatrix"))
    expect_error(lavaan2RAM(paste0("y ~ a*m\n m ~ b*x\n ", name, " := a*b")),
                 "clash")
})

test_that("p values of large z statistics are not rounded to zero", {
  set.seed(1)
  d <- data.frame(x=rnorm(500))
  d$y <- d$x + rnorm(500, sd=1.8)
  s <- summary(sem(RAM=lavaan2RAM("y ~ b*x", obs.variables=c("y", "x")), data=d))
  ## 2*(1 - pnorm(z)) is exactly 0 for z > 8.3
  z <- s$coefficients["b", "z value"]
  expect_gt(z, 9)
  expect_lt(z, 30)
  expect_equal(s$coefficients["b", "Pr(>|z|)"], 2*pnorm(-z))
  expect_gt(s$coefficients["b", "Pr(>|z|)"], 0)
})

test_that("estimates close to a bound or an inequality are not taken as restricted", {
  skip_on_cran()
  set.seed(7)
  G <- 60
  d <- data.frame(id=rep(seq_len(G), each=4), x=rnorm(4*G))
  d$m <- 0.5*d$x + rnorm(G, sd=0.4)[d$id] + rnorm(4*G)
  d$y <- 0.4*d$m + 0.2*d$x + rnorm(G, sd=0.4)[d$id] + rnorm(4*G)
  ov <- c("x", "m", "y")
  base <- "m ~ a*x\n y ~ b*m + c*x\n"
  free <- sem(RAM=lavaan2RAM(base, obs.variables=ov), data=d, cluster="id")
  theta <- coef(free)
  p <- names(theta)
  check <- function(fit) {
    expect_equal(coef(fit)[p], theta, tolerance=1e-5)
    expect_equal(sqrt(diag(vcov(fit)))[p], sqrt(diag(vcov(free))),
                 tolerance=1e-4)
    expect_equal(sqrt(diag(vcov(fit, robust=TRUE)))[p],
                 sqrt(diag(vcov(free, robust=TRUE))), tolerance=1e-4)
  }
  ## An inequality satisfied with a slack of about 5e-6
  limit <- theta[["c"]] - theta[["a"]] + 5e-6
  fit <- sem(RAM=lavaan2RAM(paste0(base, sprintf("c < a + %.12f", limit)),
                            obs.variables=ov), data=d, cluster="id")
  expect_null(.active.constraint.jacobian(fit$mx.fit, p))
  check(fit)
  ## An estimate about 5e-7 above its lower bound
  fit <- sem(RAM=lavaan2RAM(base, obs.variables=ov), data=d, cluster="id",
             lbound=list(b=theta[["b"]] - 5e-7))
  expect_true(.at.bound(fit$mx.fit)[["b"]])
  expect_silent(check(fit))
  ## A binding bound is still treated as fixed
  fit <- sem(RAM=lavaan2RAM(base, obs.variables=ov), data=d, cluster="id",
             lbound=list(b=theta[["b"]] + 0.1))
  expect_warning(se <- sqrt(diag(vcov(fit, robust=TRUE))), "at a bound")
  expect_true(is.na(se[["b"]]))
})

test_that("sem() retries a fit with a status code other than 0 or 1", {
  set.seed(4)
  d <- data.frame(x=rnorm(300))
  d$y <- 0.3*d$x + rnorm(300)
  ram <- lavaan2RAM("y ~ b*x", obs.variables=c("y", "x"))
  run <- OpenMx::mxRun
  status6 <- function(...) {
    fit <- run(...)
    fit@output$status$code <- 6L
    fit
  }
  ## The retry reaches status 0 and is kept
  local({
    local_mocked_bindings(mxRun=status6,
                          mxTryHard=function(model, ...) run(model, silent=TRUE),
                          .package="metaSEM")
    expect_no_warning(fit <- sem(RAM=ram, data=d))
    expect_equal(fit$mx.fit@output$status$code, 0)
  })
  ## The retry fails: the fit is kept with a warning
  local({
    local_mocked_bindings(mxRun=status6, mxTryHard=function(model, ...) model,
                          .package="metaSEM")
    expect_warning(fit <- sem(RAM=ram, data=d), "status code is 6")
    expect_equal(fit$mx.fit@output$status$code, 6)
  })
})

test_that("osmasem() shows warnings after the fit unless suppressWarnings=TRUE", {
  data(Nohe15A1, package="metaSEM")
  model <- "W2 + S2 ~ W1 + S1\n W1 ~~ S1\n W2 ~~ S2\n W1 ~~ 1*W1\n S1 ~~ 1*S1"
  RAM <- lavaan2RAM(model, obs.variables=c("W1", "S1", "W2", "S2"))
  my.df <- Cor2DataFrame(Nohe15A1)
  run <- OpenMx::mxRun
  local_mocked_bindings(mxRun=function(...) {
    warning("first warning")
    warning("second warning")
    run(...)
  }, .package="metaSEM")
  expect_silent(fit <- osmasem(RAM=RAM, data=my.df))
  expect_equal(fit$warning, c("first warning", "second warning"))
  expect_true(fit$mx.fit@output$status$code %in% c(0, 1))
  warnings <- character(0)
  fit <- withCallingHandlers(osmasem(RAM=RAM, data=my.df, suppressWarnings=FALSE),
                             warning=function(w) {
                               warnings <<- c(warnings, conditionMessage(w))
                               invokeRestart("muffleWarning")
                             })
  expect_equal(warnings, c("first warning", "second warning"))
  expect_true(fit$mx.fit@output$status$code %in% c(0, 1))
})

test_that("clusterRobustSE() accepts a cluster column name for multiple-group models", {
  set.seed(3)
  d <- data.frame(id=rep(seq_len(40), each=5), x=rnorm(200))
  d$y <- 0.3*d$x + rnorm(40)[d$id] + rnorm(200)
  d$g <- ifelse(d$id %% 2 == 0, "A", "B")
  d$id <- factor(paste0("c", d$id))
  fit <- sem(RAM=lavaan2RAM("y ~ c(b,b)*x", ngroups=2), data=d, group="g",
             cluster="id")
  mx <- fit$mx.fit
  groups <- sub("\\.fitfunction$", "", mx@fitfunction$groups)
  ids <- unlist(lapply(groups, function(k)
    as.character(mx[[k]]@data@observed$id)), use.names=FALSE)
  expect_equal(clusterRobustSE(mx, cluster="id"), clusterRobustSE(mx, cluster=ids))
  expect_equal(vcov(fit, robust=TRUE),
               clusterRobustSE(mx, cluster=ids, details=TRUE)$cov)
  expect_error(clusterRobustSE(mx, cluster="school"), "absent")
})

test_that("osmasem() and osmasem2() match the variables of RAM to the data by name", {
  data(Nohe15A1, package="metaSEM")
  my.df <- Cor2DataFrame(Nohe15A1)
  model <- "W2 ~ w2w*W1 + s2w*S1\n S2 ~ w2s*W1 + s2s*S1
    W1 ~~ 1*W1\n S1 ~~ 1*S1\n W1 ~~ wis1*S1\n W2 ~~ w2s2*S2"
  ram <- lavaan2RAM(model, obs.variables=c("W1", "S1", "W2", "S2"))
  ram.r <- lavaan2RAM(model, obs.variables=c("W2", "S2", "W1", "S1"))
  p <- c("w2w", "s2w", "w2s", "s2s", "wis1", "w2s2")
  expect_equal(coef(osmasem(RAM=ram.r, data=my.df))[p],
               coef(osmasem(RAM=ram, data=my.df))[p], tolerance=1e-6)
  expect_equal(coef(osmasem2(RAM=ram.r, data=my.df))[p],
               coef(osmasem2(RAM=ram, data=my.df))[p], tolerance=1e-6)

  ## RAM with other variables than the data or subset.variables
  expect_error(osmasem(RAM=ram, data=my.df, subset.variables=c("W1", "S1", "W2")),
               "do not match")
  ram3 <- lavaan2RAM("W2 ~ a*W1 + b*S1\n W1 ~~ 1*W1\n S1 ~~ 1*S1\n W1 ~~ c*S1",
                     obs.variables=c("W1", "S1", "W2"))
  expect_error(osmasem2(RAM=ram3, data=my.df), "do not match")
  ## RAM names share none with the data: matched by position with a warning
  renamed <- lapply(Nohe15A1$data, function(x) {
    dimnames(x) <- list(paste0("x", 1:4), paste0("x", 1:4)); x })
  my.df.x <- Cor2DataFrame(renamed, Nohe15A1$n)
  expect_warning(fit.x <- osmasem(RAM=ram, data=my.df.x), "matched by position")
  expect_equal(coef(fit.x)[p], coef(osmasem(RAM=ram, data=my.df))[p], tolerance=1e-6)
  expect_warning(fit.x2 <- osmasem2(RAM=ram, data=my.df.x), "matched by position")
  expect_equal(coef(fit.x2)[p], coef(osmasem2(RAM=ram, data=my.df))[p], tolerance=1e-6)
  ram.x <- ram
  rownames(ram.x$F)[1] <- "x1"
  expect_error(osmasem(RAM=ram.x, data=my.df.x), "do not match")
  ## subset.variables in another order than the data
  expect_equal(coef(osmasem(RAM=ram3, data=my.df, subset.variables=c("W2", "W1", "S1"))),
               coef(osmasem(RAM=ram3, data=my.df, subset.variables=c("W1", "S1", "W2"))),
               tolerance=1e-6)
})

test_that("osmasem() and osmasem2() fit models of correlations only as tssem1()", {
  skip_on_cran()
  data(Nohe15A1, package="metaSEM")
  v <- c("W1", "S1", "W2")
  dat <- lapply(Nohe15A1$data, function(x) x[v, v])
  my.df <- Cor2DataFrame(dat, Nohe15A1$n)
  ram <- lavaan2RAM("W1 ~~ 1*W1\n S1 ~~ 1*S1\n W2 ~~ 1*W2
    W1 ~~ r21*S1\n W1 ~~ r31*W2\n S1 ~~ r32*W2", obs.variables=v)
  ref <- tssem1(dat, Nohe15A1$n, method="REM", RE.type="Diag")
  r <- c("r21", "r31", "r32")
  fit1 <- osmasem(RAM=ram, data=my.df)
  fit2 <- osmasem2(RAM=ram, data=my.df)
  expect_equal(unname(coef(fit1)[r]), unname(coef(ref, select="fixed")), tolerance=1e-5)
  expect_equal(unname(coef(fit2)[r]), unname(coef(ref, select="fixed")), tolerance=1e-5)
  expect_equal(unname(diag(VarCorr(fit1))), unname(coef(ref, select="random")),
               tolerance=1e-4)
  expect_equal(unname(coef(fit2, select="random")), unname(coef(ref, select="random")),
               tolerance=1e-4)
  expect_equal(fit2$mx.fit$output$Minus2LogLikelihood,
               ref$mx.fit$output$Minus2LogLikelihood, tolerance=1e-6)
})

test_that("heterogeneity matrices that are not positive semidefinite give a warning", {
  expect_warning(.warn.Tau.not.psd(matrix(c(1, 2, 2, 1), 2)), "not positive semidefinite")
  expect_warning(.warn.Tau.not.psd(diag(c(0.2, -0.001))), "-0.001")
  expect_silent(.warn.Tau.not.psd(matrix(c(1, 1, 1, 1), 2)))
  expect_silent(.warn.Tau.not.psd(diag(c(0.2, -1e-9))))
  expect_silent(.warn.Tau.not.psd(matrix(0.1)))
})

test_that("osmasem2() methods identify random effects by the heterogeneity matrices", {
  data(Nohe15A1, package="metaSEM")
  my.df <- Cor2DataFrame(Nohe15A1)
  model <- "W2 ~ tau2w*W1 + s2w*S1\n S2 ~ w2s*W1 + s2s*S1
    W1 ~~ 1*W1\n S1 ~~ 1*S1\n W1 ~~ wis1*S1\n W2 ~~ w2s2*S2"
  ram <- lavaan2RAM(model, obs.variables=c("W1", "S1", "W2", "S2"))
  fit <- osmasem2(RAM=ram, data=my.df)
  expect_setequal(names(coef(fit, select="random")), paste0("tau2Cov", 1:6, "_", 1:6))
  expect_true("tau2w" %in% names(coef(fit, select="fixed")))
  expect_setequal(rownames(vcov(fit, select="fixed")), names(coef(fit, select="fixed")))
  ## No random effects
  fit0 <- osmasem2(RAM=ram, data=my.df, RE.type.Sigma="Zero")
  expect_equal(coef(fit0, select="fixed"), coef(fit0, select="all"))
  expect_length(coef(fit0, select="random"), 0)

  ## summary(fitIndices=TRUE) when the independence model cannot be fitted
  fit$mx.ind <- OpenMx::mxModel(fit$mx.ind, OpenMx::mxAlgebra(nonexistent + 1, name="z"))
  expect_warning(s <- summary(fit, fitIndices=TRUE), "independence model")
  expect_true(is.finite(s$Chi))
})

test_that("osmasem() retries a fit with a status code other than 0 or 1", {
  data(Nohe15A1, package="metaSEM")
  my.df <- Cor2DataFrame(Nohe15A1)
  RAM <- lavaan2RAM("W2 ~ W1 + S1\n W1 ~~ 1*W1\n S1 ~~ 1*S1\n W1 ~~ S1",
                    obs.variables=c("W1", "S1", "W2"))
  run <- OpenMx::mxRun
  status6 <- function(...) {
    fit <- run(...)
    fit@output$status$code <- 6L
    fit
  }
  local({
    local_mocked_bindings(mxRun=status6,
                          mxTryHard=function(model, ...) run(model, silent=TRUE),
                          .package="metaSEM")
    expect_no_warning(fit <- osmasem(RAM=RAM, data=my.df, subset.variables=c("W1", "S1", "W2")))
    expect_equal(fit$mx.fit@output$status$code, 0)
  })
  local({
    local_mocked_bindings(mxRun=status6, mxTryHard=function(model, ...) model,
                          .package="metaSEM")
    expect_warning(fit <- osmasem(RAM=RAM, data=my.df, subset.variables=c("W1", "S1", "W2")),
                   "status code is 6")
    expect_equal(fit$mx.fit@output$status$code, 6)
  })
})

test_that("osmasem2(subset.variables=) equals a fit to data with these variables only", {
  skip_on_cran()
  data(Nohe15A1, package="metaSEM")
  v <- c("W1", "S1", "W2")
  ram <- lavaan2RAM("W2 ~ a*W1 + b*S1\n W1 ~~ 1*W1\n S1 ~~ 1*S1\n W1 ~~ c*S1",
                    obs.variables=v)
  full <- Cor2DataFrame(Nohe15A1)
  reduced <- Cor2DataFrame(lapply(Nohe15A1$data, function(x) x[v, v]), Nohe15A1$n)
  fit <- osmasem2(RAM=ram, data=full, subset.variables=c("W2", "W1", "S1"))
  ref <- osmasem2(RAM=ram, data=reduced)
  expect_equal(coef(fit, select="all"), coef(ref, select="all"), tolerance=1e-6)
  expect_equal(fit$mx.fit$output$Minus2LogLikelihood,
               ref$mx.fit$output$Minus2LogLikelihood, tolerance=1e-8)
  expect_equal(fit$data$obslabels, v)

  ## Covariances and means
  data(issp05, package="metaSEM")
  v <- c("JA1", "JA2", "CS1", "CS2")
  ram2 <- lavaan2RAM("JA =~ JA1 + JA2\n CS =~ CS1 + CS2\n JA ~~ CS
    JA1 + JA2 + CS1 + CS2 ~ 1", obs.variables=v)
  full2 <- Cor2DataFrame(issp05$data, n=issp05$n, Means=issp05$means,
                         cor.analysis=FALSE)
  reduced2 <- Cor2DataFrame(lapply(issp05$data, function(x) x[v, v]), n=issp05$n,
                            Means=issp05$means[, v], cor.analysis=FALSE)
  fit2 <- osmasem2(RAM=ram2, data=full2, subset.variables=v, cor.analysis=FALSE,
                   mean.analysis=TRUE)
  ref2 <- osmasem2(RAM=ram2, data=reduced2, cor.analysis=FALSE, mean.analysis=TRUE)
  expect_equal(coef(fit2, select="all"), coef(ref2, select="all"), tolerance=1e-6)
  expect_equal(fit2$mx.fit$output$Minus2LogLikelihood,
               ref2$mx.fit$output$Minus2LogLikelihood, tolerance=1e-8)

  ## A study without any of the selected correlations does not contribute
  dat <- Nohe15A1$data
  dat[[1]][c("W2", "S2"), ] <- NA
  dat[[1]][, c("W2", "S2")] <- NA
  v <- c("W1", "W2", "S2")
  ram3 <- lavaan2RAM("W2 ~ a*W1\n S2 ~ b*W1\n W1 ~~ 1*W1\n W2 ~~ c*S2", obs.variables=v)
  ref3 <- osmasem2(RAM=ram3, data=Cor2DataFrame(lapply(Nohe15A1$data[-1], function(x) x[v, v]),
                                                Nohe15A1$n[-1]))
  for (f in list(osmasem, osmasem2)) {
    fit3 <- f(RAM=ram3, data=Cor2DataFrame(dat, Nohe15A1$n), subset.variables=v)
    expect_equal(coef(fit3$mx.fit)[c("a", "b", "c")], coef(ref3$mx.fit)[c("a", "b", "c")],
                 tolerance=1e-6)
  }

  ## Errors
  expect_error(osmasem2(RAM=ram, data=full, subset.variables=c("W1", "X9")),
               "X9 are not in the data")
  expect_error(osmasem2(RAM=ram, data=full, subset.variables=c("W1", "S1")),
               "do not match")
  expect_error(osmasem2(RAM=ram, data=full, subset.variables=c("W1", "S1", "W2", "W1")),
               "W1 more than once")
  expect_error(osmasem(RAM=ram, data=full, subset.variables=c("W1", "S1", "W2", "W1")),
               "W1 more than once")
  full$VyMlabels <- "x"
  expect_error(osmasem2(RAM=ram, data=full, subset.variables=c("W1", "S1", "W2")),
               "'VyMlabels' does not match")

  ## Sampling covariances between the covariances and the means (VyMlabels),
  ## added here as 0.05 of the product of the standard deviations
  add.VyM <- function(d) {
    nm <- outer(d$obslabels, d$ylabels, function(p, y) paste0("VyM(", p, " ", y, ")"))
    for (i in seq_along(d$obslabels)) for (j in seq_along(d$ylabels)) {
      p <- d$obslabels[i]
      y <- d$ylabels[j]
      d$data[[nm[i, j]]] <- 0.05*sqrt(d$data[[paste0("C(", p, " ", p, ")")]] *
                                        d$data[[paste0("C(", y, " ", y, ")")]])
    }
    d$VyMlabels <- c(nm)
    d
  }
  v <- c("JA1", "JA2", "CS1", "CS2")
  full4 <- add.VyM(full2)
  reduced4 <- add.VyM(reduced2)
  expect_identical(.subset.Cor2DataFrame(full4, rev(v), cor.analysis=FALSE)$VyMlabels,
                   reduced4$VyMlabels)
  fit4 <- osmasem2(RAM=ram2, data=full4, subset.variables=rev(v), cor.analysis=FALSE,
                   mean.analysis=TRUE)
  ref4 <- osmasem2(RAM=ram2, data=reduced4, cor.analysis=FALSE, mean.analysis=TRUE)
  expect_equal(coef(fit4, select="all"), coef(ref4, select="all"), tolerance=1e-6)
  expect_equal(fit4$mx.fit$output$Minus2LogLikelihood,
               ref4$mx.fit$output$Minus2LogLikelihood, tolerance=1e-8)
})

test_that("two-level models estimate the means of within-only variables", {
  ## x appears only at level 1 and has a mean of 3. Its mean used to be
  ## fixed at 0, so its variance absorbed the squared mean.
  set.seed(1)
  id <- rep(seq_len(50), each=6)
  d <- data.frame(id=id, x=rnorm(300, mean=3))
  d$y <- 1 + .5*d$x + rnorm(50)[id] + rnorm(300)
  model <- "level: 1
  y ~ b*x
level: 2
  y ~ 1
  y ~~ y"
  ram <- lavaan2RAM(model)
  expect_identical(ram$within$M[1, "x"], "0*xmean_within")
  expect_identical(ram$within$M[1, "y"], "0")
  fit <- sem(RAM=ram, data=d, cluster="id")
  expect_equal(unname(coef(fit)["xmean_within"]), mean(d$x), tolerance=1e-4)
  expect_equal(unname(coef(fit)["xWITHx_within"]), mean((d$x-mean(d$x))^2),
               tolerance=1e-4)
  ## ML estimates are invariant to shifting x; with x centred at its mean,
  ## they agree with lavaan.
  d$x <- d$x - mean(d$x)
  centred <- sem(RAM=ram, data=d, cluster="id")
  expect_equal(coef(centred)["b"], coef(fit)["b"], tolerance=1e-5)
  reference <- lavaan::sem(model, data=d, cluster="id", fixed.x=FALSE)
  expect_equal(unname(coef(centred)["b"]), unname(lavaan::coef(reference)["b"]),
               tolerance=1e-4)
  expect_equal(centred$mx.fit@output$Minus2LogLikelihood,
               -2*as.numeric(lavaan::logLik(reference)), tolerance=1e-6)
})

test_that("two-level fits of Demo.twolevel agree with lavaan", {
  model <- "level: 1
  fw =~ y1 + y2 + y3
  fw ~ x1 + x2 + x3
level: 2
  fb =~ y1 + y2 + y3
  fb ~ w1 + w2"
  reference <- lavaan::sem(model, data=lavaan::Demo.twolevel, cluster="cluster",
                           fixed.x=FALSE, std.lv=TRUE)
  fit <- sem(RAM=lavaan2RAM(model), data=lavaan::Demo.twolevel, cluster="cluster")
  expect_equal(length(coef(fit)), as.integer(lavaan::lavInspect(reference, "npar")))
  expect_equal(fit$mx.fit@output$Minus2LogLikelihood,
               -2*as.numeric(lavaan::logLik(reference)), tolerance=1e-6)
  ## The number of observations and BIC use the level-1 rows, as in lavaan.
  s <- summary(fit)
  expect_equal(s$no.studies, nrow(lavaan::Demo.twolevel))
  expect_equal(unname(s$informationCriteria["BIC:", "Parameters Penalty"]),
               BIC(reference), tolerance=1e-6)
  ## An explicit "~1" makes lavaanify() add fixed means for all other
  ## variables; these must not override the free between-level means.
  explicit <- sem(RAM=lavaan2RAM(paste(model, "\n  y1 ~ 1")),
                  data=lavaan::Demo.twolevel, cluster="cluster")
  expect_setequal(names(coef(explicit)), names(coef(fit)))
  expect_equal(explicit$mx.fit@output$Minus2LogLikelihood,
               fit$mx.fit@output$Minus2LogLikelihood, tolerance=1e-6)
})

test_that("likelihood-based intervals match defined parameters by their full names", {
  set.seed(3)
  x <- rnorm(200)
  m <- .5*x + rnorm(200)
  y <- .4*m + .2*x + rnorm(200)
  ram <- lavaan2RAM("m ~ a*x\n y ~ b*m + c*x\n ab := a*b\n ab2 := a*b + c",
                    obs.variables=c("y", "m", "x"))
  fit <- sem(RAM=ram, data=data.frame(x, m, y), intervals.type="LB")
  s <- summary(fit)
  expect_identical(rownames(s$mxalgebras), c("ab", "ab2"))
  expect_equal(s$mxalgebras["ab", "estimate"], unname(coef(fit)["a"]*coef(fit)["b"]),
               tolerance=1e-6)
  expect_equal(s$mxalgebras["ab2", "estimate"],
               unname(coef(fit)["a"]*coef(fit)["b"] + coef(fit)["c"]), tolerance=1e-6)
  ## The coefficients keep the order of the parameters in the model.
  expect_identical(rownames(s$coefficients), names(coef(fit)))
})

test_that("plot() draws combined group and two-level fits", {
  skip_if_not_installed("semPlot")
  set.seed(4)
  d <- do.call(rbind, lapply(c("A", "B"), function(g) {
    id <- rep(seq_len(20), each=5)
    x <- rnorm(100)
    data.frame(g=g, id=paste0(g, id), x=x, y=1 + .4*x + rnorm(20)[id] + rnorm(100))
  }))
  model <- "group: A
  level: 1
    y ~ x
  level: 2
    y ~~ y
group: B
  level: 1
    y ~ x
  level: 2
    y ~~ y"
  fit <- sem(RAM=lavaan2RAM(model, ngroups=2), data=d, group="g", cluster="id")
  grDevices::pdf(nullfile())
  on.exit(grDevices::dev.off(), add=TRUE)
  expect_error(plot(fit), NA)
  expect_s3_class(plot(fit, group="B"), "qgraph")
  expect_s3_class(plot(fit, group="A", level="between"), "qgraph")
  expect_error(plot(fit, group="C"), "'group' must be one of")
})

test_that("between-level plots draw cluster-level covariates as observed", {
  skip_if_not_installed("semPlot")
  model <- "level: 1
  fw =~ y1 + y2 + y3
level: 2
  fb =~ y1 + y2 + y3
  fb ~ w1"
  fit <- sem(RAM=lavaan2RAM(model), data=lavaan::Demo.twolevel, cluster="cluster")
  grDevices::pdf(nullfile())
  on.exit(grDevices::dev.off(), add=TRUE)
  p <- plot(fit, level="between")
  expect_setequal(p$graphAttributes$Nodes$labels[p$graphAttributes$Nodes$shape == "square"],
                  c("y1", "y2", "y3", "w1"))
})

test_that("per-cluster summaries are reconstructed whatever the units", {
  ## A variance of 1e-12 used to be dropped as 0 by an absolute tolerance.
  for (s in c(1, 1e-8, 1e-12, 1e-20)) {
    cov.list <- list(a=matrix(c(s, s/2, s/2, s), 2, dimnames=list(c("y1", "y2"), c("y1", "y2"))),
                     b=matrix(c(2*s, 0, 0, s), 2, dimnames=list(c("y1", "y2"), c("y1", "y2"))))
    ## A change of units scales the means by sqrt(s) too.
    means <- list(c(y1=0, y2=1)*sqrt(s), c(y1=2, y2=3)*sqrt(s))
    out <- .expand.cluster.summaries(cov.list, means, c(5, 6), "id")
    for (i in 1:2) {
      rows <- out$id == names(cov.list)[i]
      ## Divided by s: expect_equal() compares values below its tolerance
      ## absolutely, so 0 would pass for 1e-12.
      expect_equal(cov(out[rows, c("y1", "y2")])/s, cov.list[[i]]/s, tolerance=1e-8)
      expect_equal(colMeans(out[rows, c("y1", "y2")])/sqrt(s), means[[i]]/sqrt(s),
                   tolerance=1e-8)
    }
  }
  ## The summary-input likelihood equals the raw-data likelihood in small units.
  set.seed(5)
  id <- rep(seq_len(30), each=4)
  d <- data.frame(id=id, y1=rnorm(30)[id] + rnorm(120), y2=rnorm(120))
  d[, c("y1", "y2")] <- d[, c("y1", "y2")]*1e-6
  ram <- lavaan2RAM("level: 1\n y1 ~~ y2\n level: 2\n y1 ~~ y2")
  raw <- sem(RAM=ram, data=d, cluster="id")
  pieces <- split(d[, c("y1", "y2")], d$id)
  summ <- sem(RAM=ram, Cov=lapply(pieces, cov), means=lapply(pieces, colMeans),
              numObs=sapply(pieces, nrow), cluster="id", run=FALSE)$mx.fit
  summ <- mxRun(omxSetParameters(summ, labels=names(coef(raw)), values=coef(raw)),
                useOptimizer=FALSE, silent=TRUE)
  expect_equal(summ@output$fit, raw$mx.fit@output$fit, tolerance=1e-8)
})

test_that("Cor2DataFrame() selects the moderators by name", {
  x1 <- Nohe15A1[c("data", "n", "Lag", "Country")]
  x2 <- Nohe15A1[c("Lag", "n", "Country", "data")]
  out1 <- Cor2DataFrame(x1)
  out2 <- Cor2DataFrame(x2)
  expect_identical(out1$data$Lag, Nohe15A1$Lag)
  expect_identical(out2$data$Lag, Nohe15A1$Lag)
  expect_identical(out1$data[, names(out1$data)], out2$data[, names(out1$data)])
  ## Row names are the names of the studies, made unique.
  named <- Cor2DataFrame(Nohe15A1, row.names.unique=TRUE)
  expect_identical(rownames(named$data), make.names(names(Nohe15A1$data), unique=TRUE))
})

test_that("defined parameters named like constraints are reported", {
  ## lavaan2RAM() accepts "constraint1_effect"; sem() used to take it for a
  ## constraint because of its name and dropped it from summaries.
  set.seed(1)
  d <- data.frame(x=rnorm(100), g=rep(c("A", "B"), each=50), id=rep(1:20, each=5))
  d$y <- .5*d$x + rnorm(20)[d$id] + rnorm(100)
  fit <- sem(RAM=lavaan2RAM("y ~ b*x\n y ~~ s*y\n constraint1_effect := b*b\n s < 10*b + 100",
                            obs.variables=c("y", "x")), data=d)
  expect_length(fit$mx.fit@constraints, 1L)
  expect_identical(rownames(summary(fit)$mxalgebras), "constraint1_effect")
  expect_equal(summary(fit)$mxalgebras[1, "estimate"], unname(coef(fit)["b"]^2),
               tolerance=1e-6)
  group <- sem(RAM=lavaan2RAM("y ~ c(b, b2)*x\n constraint1_effect := b - b2", ngroups=2),
               data=d, group="g")
  expect_identical(rownames(summary(group)$mxalgebras), "constraint1_effect")
  twolevel <- sem(RAM=lavaan2RAM("level: 1\n y ~ b*x\n level: 2\n y ~~ y\n constraint1_effect := 2*b"),
                  data=d, cluster="id")
  expect_identical(rownames(summary(twolevel)$mxalgebras), "constraint1_effect")
  ## wls() (tssem2) and osmasem() take the same lavaan2RAM() algebras.
  ram <- lavaan2RAM("Agreeableness =~ a1*A + a2*C + ES\n Plasticity =~ E + I
                     Agreeableness ~~ Plasticity\n constraint1_effect := a1*a2",
                    obs.variables=colnames(Digman97$data[[1]]))
  stage1 <- tssem1(Digman97$data, Digman97$n, method="FEM")
  stage2 <- tssem2(stage1, RAM=ram)
  expect_true("constraint1_effect" %in% stage2$mx.algebras)
  one <- osmasem(RAM=ram, data=Cor2DataFrame(Digman97), run=FALSE)
  expect_true("constraint1_effect" %in% names(one$mx.fit@intervals))
})

test_that("plot() shows paths set by definition variables at their means", {
  skip_if_not_installed("semPlot")
  grDevices::pdf(nullfile())
  on.exit(grDevices::dev.off(), add=TRUE)
  ## Weight of the variance edge of y in a plot
  y.var <- function(p) {
    labels <- p$graphAttributes$Nodes$labels
    unname(p$Edgelist$weight[p$Edgelist$bidirectional & labels[p$Edgelist$from] == "y" &
                               labels[p$Edgelist$to] == "y"])
  }
  set.seed(1)
  id <- rep(seq_len(40), each=6)
  d <- data.frame(id=id, z=rnorm(40)[id], w=runif(240))
  d$y <- rnorm(40)[id] + rnorm(240)*sqrt(exp(.5*d$w))
  ## A substituted constraint with a definition variable at either level used
  ## to be shown at the first data row.
  between <- sem(RAM=lavaan2RAM("level: 1\n y ~~ sw*y\n level: 2\n y ~~ sb*y
                                 sb == exp(a0 + a1*data.z)"),
                 data=d, cluster="id", replace.constraints=TRUE)
  est <- coef(between)
  expect_equal(y.var(plot(between, level="between")),
               unname(exp(est["a0"] + est["a1"]*mean(d$z))), tolerance=1e-6)
  within <- sem(RAM=lavaan2RAM("level: 1\n y ~~ sw*y\n level: 2\n y ~~ sb*y
                                sw == exp(b0 + b1*data.w)"),
                data=d, cluster="id", replace.constraints=TRUE)
  est <- coef(within)
  expect_equal(y.var(plot(within, level="within")),
               unname(exp(est["b0"] + est["b1"]*mean(d$w))), tolerance=1e-6)
  ## In a multiple-group plot, the mean is that of the group plotted.
  d$g <- rep(c("A", "B"), each=120)
  d$w[d$g == "B"] <- d$w[d$g == "B"] + 2
  plain <- sem(RAM=lavaan2RAM("y ~~ c(data.w, data.w)*y", obs.variables="y", ngroups=2),
               data=d, group="g")
  expect_equal(y.var(plot(plain, group="B")), mean(d$w[d$g == "B"]))
  substituted <- sem(RAM=lavaan2RAM("y ~~ c(s1, s2)*y
                                     s1 == exp(c0 + c1*data.w)
                                     s2 == exp(c0 + c1*data.w)", ngroups=2),
                     data=d, group="g", replace.constraints=TRUE)
  est <- coef(substituted)
  expect_equal(y.var(plot(substituted, group="B")),
               unname(exp(est["c0"] + est["c1"]*mean(d$w[d$g == "B"]))),
               tolerance=1e-6)
})

test_that("WLS and TSSEM align named RAM observations and subset weights", {
  vars <- c("x", "y", "z")
  S <- matrix(c(1, .5, .2, .5, 1, .1, .2, .1, 1), 3,
              dimnames=list(vars, vars))
  ram <- lavaan2RAM("y ~ a*x\nz ~ 0*x\ny ~~ 0*z", obs.variables=vars)
  reversed <- ram
  reversed$F <- ram$F[rev(vars), , drop=FALSE]

  for (cor in c(FALSE, TRUE)) {
    ## Fix exogenous variances for the correlation model.
    model <- ram
    if (cor) {
      model$S["x", "x"] <- "1"
      model$S["z", "z"] <- "1"
    }
    reversed <- model
    reversed$F <- model$F[rev(vars), , drop=FALSE]
    k <- if (cor) 3 else 6
    V <- diag(seq_len(k)/1000)
    stage <- tssem1(rep(list(S), 3), n=rep(200, 3), method="FEM",
                    cor.analysis=cor)
    for (method in c("wls", "tssem2")) {
      fit <- function(R, subset=NULL) {
        if (method == "wls")
          wls(S, V, n=600, RAM=R, cor.analysis=cor, subset.variables=subset)
        else tssem2(stage, RAM=R, subset.variables=subset)
      }
      reference <- fit(model)
      reordered <- fit(reversed)
      subset <- fit(model, rev(vars))
      for (other in list(reordered, subset)) {
        expect_true(reference$mx.fit$output$status$code %in% c(0, 1))
        expect_true(other$mx.fit$output$status$code %in% c(0, 1))
        expect_equal(coef(other), coef(reference), tolerance=1e-5)
        expect_equal(other$mx.fit$output$fit, reference$mx.fit$output$fit,
                     tolerance=1e-5)
      }
    }
  }
  other.names <- ram
  rownames(other.names$F) <- c("a", "b", "c")
  expect_warning(wls(S, diag(6), 600, RAM=other.names, cor.analysis=FALSE,
                     run=FALSE), "matching by position")
})

test_that("direct OSMASEM matrices align observations without losing parameters", {
  data(Nohe15A1, package="metaSEM")
  d <- Cor2DataFrame(Nohe15A1)
  syntax <- "W2 ~ a*W1 + b*S1\nS2 ~ c*W1 + e*S1\nW1 ~~ 1*W1\nS1 ~~ 1*S1\nW1 ~~ f*S1\nW2 ~~ g*S2"
  ram <- lavaan2RAM(syntax, obs.variables=rev(d$obslabels))
  m <- create.vechsR(ram$A, ram$S)
  tau <- create.Tau2(RAM=ram, RE.type="Diag")
  reference <- osmasem(RAM=ram, data=d)
  direct <- osmasem(Mmatrix=m, Tmatrix=tau, data=d)
  expect_true(direct$mx.fit$output$status$code %in% c(0, 1))
  expect_equal(coef(direct)[names(coef(reference))], coef(reference), tolerance=1e-5)
  expect_equal(direct$mx.fit$output$fit, reference$mx.fit$output$fit, tolerance=1e-5)
  m$Fmatrix$free[1, 1] <- TRUE
  m$Fmatrix$labels[1, 1] <- "filter_parameter"
  m$Fmatrix$lbound[1, 1] <- -.5
  m$Fmatrix$ubound[1, 1] <- 2
  output <- osmasem(Mmatrix=m, Tmatrix=tau, data=d, run=FALSE)$mx.fit$Fmatrix
  i <- match(d$obslabels, rownames(m$Fmatrix$values))
  expect_equal(unname(output$free), unname(m$Fmatrix$free[i, , drop=FALSE]))
  expect_equal(unname(output$labels), unname(m$Fmatrix$labels[i, , drop=FALSE]))
  expect_equal(unname(output$lbound), unname(m$Fmatrix$lbound[i, , drop=FALSE]))
  expect_equal(unname(output$ubound), unname(m$Fmatrix$ubound[i, , drop=FALSE]))
})

test_that("WLS plot preserves endpoints when RAM variables are reordered", {
  skip_if_not_installed("semPlot")
  vars <- c("x", "y", "z")
  S <- matrix(c(1, .5, .2, .5, 1, .1, .2, .1, 1), 3,
              dimnames=list(vars, vars))
  syntax <- "y ~ a*x\nz ~ b*y\nx ~~ 0*z"
  ram <- lavaan2RAM(syntax, obs.variables=vars)
  revram <- lavaan2RAM(syntax, obs.variables=rev(vars))
  fit <- wls(S, diag(6)/1000, 600, RAM=ram, cor.analysis=FALSE)
  reversed <- wls(S, diag(6)/1000, 600, RAM=revram, cor.analysis=FALSE)
  expect_true(fit$mx.fit$output$status$code %in% c(0, 1))
  expect_true(reversed$mx.fit$output$status$code %in% c(0, 1))
  p1 <- meta2semPlot(fit)@Pars
  p2 <- meta2semPlot(reversed)@Pars
  i <- match(c("a", "b"), p1$label)
  j <- match(c("a", "b"), p2$label)
  expect_false(anyNA(i)); expect_false(anyNA(j))
  left <- p1[i, c("lhs", "edge", "rhs", "est")]
  right <- p2[j, c("lhs", "edge", "rhs", "est")]
  rownames(left) <- rownames(right) <- NULL
  expect_equal(left, right, tolerance=1e-5)
})

test_that("OSMASEM plots use data and custom names for RAM nodes", {
  skip_if_not_installed("semPlot")
  data(Nohe15A1, package="metaSEM")
  model <- "W2 ~ a*W1 + b*S1\nS2 ~ c*W1 + e*S1\nW1 ~~ 1*W1\nS1 ~~ 1*S1\nW1 ~~ f*S1\nW2 ~~ g*S2"
  ram <- lavaan2RAM(model, obs.variables=c("W2", "S2", "W1", "S1"))
  data <- Cor2DataFrame(Nohe15A1)
  renamed <- lapply(Nohe15A1$data, function(x) {
    dimnames(x) <- rep(list(paste0("x", 1:4)), 2); x
  })
  renamed <- Cor2DataFrame(renamed, Nohe15A1$n)
  testthat::local_mocked_bindings(semPaths=function(object, ...) object, .package="semPlot")
  for (method in c("osmasem", "osmasem2")) {
    fit <- do.call(method, list(RAM=ram, data=data))
    fit.x <- suppressWarnings(do.call(method, list(RAM=ram, data=renamed)))
    expect_true(fit$mx.fit$output$status$code %in% c(0, 1))
    expect_true(fit.x$mx.fit$output$status$code %in% c(0, 1))
    positional <- plot(fit.x)
    expect_setequal(positional@Vars$name, paste0("x",1:4))
    expect_true(any(positional@Pars$lhs == "x3" & positional@Pars$edge == "->" &
                    positional@Pars$rhs == "x1"))
    for (spec in list(list(fit=fit.x, names=paste0("x",1:4)),
                      list(fit=fit, names=paste0("custom",1:4)))) {
      plotted <- plot(spec$fit, manNames=spec$names)
      expect_true(all(c(plotted@Pars$lhs, plotted@Pars$rhs) %in% c(spec$names, "")))
      expect_setequal(plotted@Vars$name, spec$names)
      source <- if (identical(spec$names, paste0("x",1:4))) 3 else 1
      target <- if (source == 3) 1 else 3
      expect_true(any(plotted@Pars$lhs == spec$names[source] & plotted@Pars$edge == "->" &
                      plotted@Pars$rhs == spec$names[target]))
    }
    if (method == "osmasem2") {
      fit$RAM <- NULL
      expect_s4_class(plot(fit, manNames=paste0("custom",1:4)), "semPlotModel")
    }
  }
})

test_that("different observation counts fail early and partial names warn", {
  vars <- c("x", "y", "z")
  C <- matrix(c(1,.5,.2,.5,1,.1,.2,.1,1),3,dimnames=list(vars,vars))
  ram <- lavaan2RAM("y ~ a*x\nz ~ 0*x\ny ~~ 0*z",obs.variables=vars)
  bad <- ram; bad$F <- bad$F[-1,,drop=FALSE]
  expect_error(wls(C,diag(6),600,RAM=bad,cor.analysis=FALSE,run=FALSE), "RAM has 2 observed")
  partial <- ram; rownames(partial$F)[1] <- "typo"
  expect_warning(wls(C,diag(6),600,RAM=partial,cor.analysis=FALSE,run=FALSE), "possible typos")
  data(Nohe15A1,package="metaSEM")
  d <- Cor2DataFrame(Nohe15A1)
  mram <- lavaan2RAM("W2 ~ a*W1 + b*S1\nS2 ~ c*W1 + e*S1\nW1 ~~ 1*W1\nS1 ~~ 1*S1",obs.variables=d$obslabels)
  m <- create.vechsR(mram$A,mram$S)
  tau <- create.Tau2(RAM=mram,RE.type="Diag")
  badm <- m; badm$Fmatrix <- m$Fmatrix[-1,,drop=FALSE]
  expect_error(osmasem(Mmatrix=badm,Tmatrix=tau,data=d,run=FALSE), "Mmatrix has 3 observed")
  partialm <- m; dimnames(partialm$Fmatrix@values)[[1]][1] <- "typo"
  expect_warning(osmasem(Mmatrix=partialm,Tmatrix=tau,data=d,run=FALSE), "possible typos")
})

test_that("reordered RAM models preserve WLS and TSSEM fits", {
  skip_on_cran()
  vars <- c("x","y","z")
  C <- matrix(c(1,.3,.12,.3,1,.4,.12,.4,1),3,dimnames=list(vars,vars))
  syntax <- "y ~ a*x\nz ~ b*y\nx ~~ 1*x"
  r1 <- lavaan2RAM(syntax,obs.variables=vars)
  r2 <- lavaan2RAM(syntax,obs.variables=rev(vars))
  stage <- tssem1(rep(list(C),3),n=rep(200,3),method="FEM")
  for (method in c("wls","tssem2")) for (dc in c(FALSE,TRUE)) {
    fit <- function(r) if (method == "wls") wls(C,diag(3)/1000,600,RAM=r,diag.constraints=dc) else tssem2(stage,RAM=r,diag.constraints=dc)
    a <- fit(r1); b <- fit(r2)
    expect_true(a$mx.fit$output$status$code %in% c(0,1))
    expect_true(b$mx.fit$output$status$code %in% c(0,1))
    expect_equal(coef(b)[names(coef(a))],coef(a),tolerance=1e-5)
    expect_equal(b$mx.fit$output$fit,a$mx.fit$output$fit,tolerance=1e-5)
  }
})

test_that("OSMASEM moderator and means respect reordered RAM variables", {
  skip_on_cran()
  data(Nohe15A1,package="metaSEM")
  d <- Cor2DataFrame(Nohe15A1)
  d$data$Lag <- as.numeric(scale(Nohe15A1$Lag))
  syntax <- "W2 ~ a*W1 + b*S1\nS2 ~ c*W1 + e*S1\nW1 ~~ 1*W1\nS1 ~~ 1*S1\nW1 ~~ f*S1\nW2 ~~ g*S2"
  r1 <- lavaan2RAM(syntax,obs.variables=d$obslabels)
  r2 <- lavaan2RAM(syntax,obs.variables=rev(d$obslabels))
  a <- osmasem(RAM=r1,Ax=create.modMatrix(r1,output="A",mod="Lag"),RE.type="Zero",data=d)
  b <- osmasem(RAM=r2,Ax=create.modMatrix(r2,output="A",mod="Lag"),RE.type="Zero",data=d)
  expect_true(a$mx.fit$output$status$code %in% c(0,1)); expect_true(b$mx.fit$output$status$code %in% c(0,1))
  expect_equal(coef(b)[names(coef(a))],coef(a),tolerance=1e-5)
  expect_equal(b$mx.fit$output$fit,a$mx.fit$output$fit,tolerance=1e-5)
  vars <- c("x","y","z")
  C <- matrix(c(1,.2,.3,.2,1,.4,.3,.4,1),3,dimnames=list(vars,vars))
  M <- matrix(rep(c(.2,.5,.8),6),ncol=3,byrow=TRUE,dimnames=list(NULL,vars))
  d <- Cor2DataFrame(rep(list(C),6),n=rep(100,6),cor.analysis=FALSE,Means=M)
  syntax <- "y ~ a*x\nz ~ b*y\nx ~ mx*1\ny ~ my*1\nz ~ mz*1"
  r1 <- lavaan2RAM(syntax,obs.variables=vars); r2 <- lavaan2RAM(syntax,obs.variables=rev(vars))
  a <- osmasem2(RAM=r1,data=d,cor.analysis=FALSE,mean.analysis=TRUE,RE.type.Sigma="Zero",RE.type.Mu="Zero")
  b <- osmasem2(RAM=r2,data=d,cor.analysis=FALSE,mean.analysis=TRUE,RE.type.Sigma="Zero",RE.type.Mu="Zero")
  expect_true(a$mx.fit$output$status$code %in% c(0,1)); expect_true(b$mx.fit$output$status$code %in% c(0,1))
  expect_equal(coef(b)[names(coef(a))],coef(a),tolerance=1e-5)
  expect_equal(b$mx.fit$output$fit,a$mx.fit$output$fit,tolerance=1e-5)
})


test_that("OSMASEM plots honor custom names for latent and observed nodes", {
  skip_if_not_installed("semPlot")
  data(Nohe15A1, package="metaSEM")
  data <- Cor2DataFrame(Nohe15A1)
  ram <- lavaan2RAM("F =~ 1*W1 + l2*S1 + l3*W2 + l4*S2\nF ~~ 1*F",
                    obs.variables=rev(data$obslabels))
  testthat::local_mocked_bindings(semPaths=function(object, ...) object, .package="semPlot")
  for (method in c("osmasem", "osmasem2")) {
    fit <- suppressWarnings(if (method == "osmasem") osmasem(RAM=ram, data=data, RE.type="Zero")
                            else osmasem2(RAM=ram, data=data, RE.type.Sigma="Zero"))
    expect_true(fit$mx.fit$output$status$code %in% c(0, 1))
    p <- plot(fit, manNames=paste0("M", 1:4), latNames="Lat")
    expect_setequal(p@Vars$name, c(paste0("M", 1:4), "Lat"))
    expect_true(all(c(p@Pars$lhs, p@Pars$rhs) %in% c(paste0("M", 1:4), "Lat", "")))
  }
})
test_that("smdMTS does not replace an incomplete control standardizer", {
  m <- c(0, 0, 2)
  v <- c(100, 1, 1)
  n <- rep(100, 3)
  complete <- suppressWarnings(smdMTS(m, v, n, homogeneity="none",
                                      all.comparisons=TRUE, bias.adjust=FALSE))
  expect_equal(unname(complete$y["y3_2"]), 0.2, tolerance=1e-6)
  for (bad in list(list(m=c(NA, 0, 2), v=v, n=n),
                   list(m=m, v=c(NA, 1, 1), n=n),
                   list(m=m, v=v, n=c(NA, 100, 100)))) {
    expect_error(smdMTS(bad$m, bad$v, bad$n, homogeneity="none",
                        all.comparisons=TRUE), "first group")
  }
  pooled <- smdMTS(c(NA, 0, 2), v, n, homogeneity="variance",
                   all.comparisons=TRUE, bias.adjust=FALSE)
  expect_true(is.finite(pooled$y["y3_2"]))
})

test_that("uniR2mx aligns RAM observations with the pooled matrix", {
  vars <- c("x", "y", "z")
  C <- matrix(c(1, .5, .05, .5, 1, .3, .05, .3, 1), 3,
              dimnames=list(vars, vars))
  first <- uniR1(list(C, C), c(200, 250))
  syntax <- "y ~ b*x\nz ~ c*y"
  ordered <- lavaan2RAM(syntax, obs.variables=vars)
  reversed <- lavaan2RAM(syntax, obs.variables=rev(vars))
  unnamed <- reversed
  dimnames(unnamed$F) <- NULL
  fits <- lapply(list(ordered, reversed, unnamed), function(ram) uniR2mx(first, RAM=ram))
  reference <- uniR2lavaan(first, syntax)
  for (fit in fits) {
    expect_equal(fit$output$status$code, 0L)
    expect_equal(unname(fit$output$estimate[c("b", "c")]),
                 unname(lavaan::coef(reference)[c("b", "c")]), tolerance=1e-5)
    expect_equal(fit$output$fit, fits[[1]]$output$fit, tolerance=1e-6)
  }
  direct <- uniR2mx(first, Amatrix=reversed$A, Smatrix=reversed$S,
                    Fmatrix=reversed$F)
  expect_equal(unname(direct$output$estimate[c("b", "c")]),
               unname(fits[[1]]$output$estimate[c("b", "c")]), tolerance=1e-5)
  expect_equal(direct$output$fit, fits[[1]]$output$fit, tolerance=1e-6)
  unnamed.data <- first
  dimnames(unnamed.data$r.mean) <- NULL
  expect_equal(uniR2mx(unnamed.data, RAM=ordered)$output$fit,
               fits[[1]]$output$fit, tolerance=1e-6)
  different.names <- first
  dimnames(different.names$r.mean) <- list(letters[1:3], letters[1:3])
  expect_warning(uniR2mx(different.names, RAM=ordered, run=FALSE), "matching by position")
  expect_error(uniR2mx(first, RAM=lavaan2RAM("y ~ x", obs.variables=c("x", "y")),
                       run=FALSE), "observed variables")
})

test_that("uniR2mx aligns observed variables in models with latent factors", {
  vars <- paste0("x", 1:4)
  loadings <- c(.4, .5, .6, .7)
  C <- tcrossprod(loadings) + diag(1-loadings^2)
  dimnames(C) <- list(vars, vars)
  first <- uniR1(list(C, C), c(200, 250))
  syntax <- "f =~ x1 + x2 + x3 + x4"
  ordered <- uniR2mx(first, RAM=lavaan2RAM(syntax, obs.variables=vars))
  reversed <- uniR2mx(first, RAM=lavaan2RAM(syntax, obs.variables=rev(vars)))
  expect_equal(c(ordered$output$status$code, reversed$output$status$code), c(0L, 0L))
  expect_equal(reversed$output$fit, ordered$output$fit, tolerance=1e-6)
  expect_equal(reversed$output$estimate[names(ordered$output$estimate)],
               ordered$output$estimate, tolerance=1e-5)
})

test_that("rCorSam rejects unintended recycling", {
  A <- diag(2)
  B <- matrix(c(1, .8, .8, 1), 2)
  expect_error(rCorSam(list(A, B), rep(20, 4)), "equal lengths")
  expect_error(rCorSam(list(A, B), numeric()), "equal lengths")
  set.seed(42)
  one <- rCorSam(A, c(20, 30))
  set.seed(42)
  reference <- list(rCorSam(A, 20), rCorSam(A, 30))
  expect_equal(one, reference)
  set.seed(42)
  one.n <- rCorSam(list(A, B), 20)
  set.seed(42)
  reference.n <- list(rCorSam(A, 20), rCorSam(B, 20))
  expect_equal(one.n, reference.n)
})
