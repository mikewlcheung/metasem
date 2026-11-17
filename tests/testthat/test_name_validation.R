test_that("lavaan syntax rejects dotted variables and labels with useful advice", {
  expect_error(lavaan2RAM("y.1 ~ x"), "y\\.1.*y_1")
  expect_error(lavaan2RAM("f.1 =~ x1 + x2"), "f\\.1.*f_1")
  expect_error(lavaan2RAM("y ~ b.1*x"), "b\\.1.*b_1")
  expect_error(lavaan2RAM("y ~ bad.label*x"), "bad\\.label.*bad_label")
  expect_error(lavaan2RAM("y ~ b*x\nq.1 := b*b"), "q\\.1.*q_1")
  expect_error(lavaan2RAM("y ~~ data.v.1*y"), "v\\.1.*v_1")
  expect_error(lavaan2RAM("y ~ x", obs.variables=c("y", "x.1")), "x\\.1.*x_1")
  expect_match(lavaan2RAM("y ~~ data.v_1*y")$S[1, 1], "data.v_1", fixed=TRUE)
})

test_that("SEM and WLS reject dotted data and manually specified RAM names", {
  ram <- lavaan2RAM("y ~ b*x", obs.variables=c("x", "y"))
  cov <- matrix(c(1, .3, .3, 1), 2, dimnames=list(c("x", "y"), c("x", "y")))
  bad.cov <- cov
  dimnames(bad.cov) <- list(c("x.1", "y"), c("x.1", "y"))
  expect_error(wls(bad.cov, matrix(1, 1), n=100, RAM=ram, run=FALSE), "x\\.1.*x_1")
  expect_error(wls(cov, matrix(1, 1), n=100, RAM=ram, subset.variables="x.1", run=FALSE), "x\\.1.*x_1")
  bad.ram <- ram
  bad.ram$A["y", "x"] <- "0.1*b.1"
  expect_error(wls(cov, matrix(1, 1), n=100, RAM=bad.ram, run=FALSE), "b\\.1.*b_1")
  expect_error(sem(RAM=bad.ram, Cov=cov, numObs=100, run=FALSE), "b\\.1.*b_1")
  bad.names <- ram
  colnames(bad.names$F)[1] <- "x.1"
  expect_error(wls(cov, matrix(1, 1), n=100, RAM=bad.names, run=FALSE), "x\\.1.*x_1")
  expect_error(wls(cov, matrix(1, 1), n=100, Amatrix=bad.ram$A,
                   Smatrix=ram$S, Fmatrix=ram$F, run=FALSE), "b\\.1.*b_1")
  expect_error(as.mxMatrix(bad.ram$A), "b\\.1.*b_1")
  old.stage <- structure(list(original.names=c("x.1", "y")), class="tssem1REM")
  expect_error(tssem2(old.stage, RAM=ram, run=FALSE), "x\\.1.*x_1")
  data <- data.frame(x=1:4, y=2:5, check.names=FALSE)
  data[["unused.1"]] <- 1:4
  expect_error(sem(RAM=ram, data=data, run=FALSE), "unused\\.1.*unused_1")
  expect_error(sem(RAM=ram, Cov=bad.cov, numObs=100, run=FALSE), "x\\.1.*x_1")
  expect_error(sem(RAM=ram, Cov=cov, numObs=100, group="site.1", run=FALSE), "site\\.1.*site_1")
  expect_no_error(suppressWarnings(wls(cov, matrix(1, 1), n=100, RAM=ram, run=FALSE)))
})

test_that("study and one-stage input checks preserve identifiers and definition syntax", {
  cov <- matrix(c(1, .2, .2, 1), 2, dimnames=list(c("x", "y"), c("x", "y")))
  studies <- list("study.1"=cov, "study.2"=cov)
  bad <- studies
  dimnames(bad[[2]]) <- list(c("x.1", "y"), c("x.1", "y"))
  expect_error(tssem1(bad, c(100, 100), run=FALSE), "x\\.1.*x_1")
  expect_error(Cor2DataFrame(bad, c(100, 100)), "x\\.1.*x_1")
  expect_error(asyCov(bad, c(100, 100)), "x\\.1.*x_1")
  expect_error(rCorPop(as.data.frame(bad[[2]]), diag(1), k=1), "x\\.1.*x_1")
  expect_error(uniR1(bad, c(100, 100)), "x\\.1.*x_1")
  pooled <- uniR1(studies, c(100, 100))
  bad.pooled <- pooled
  dimnames(bad.pooled$r.mean) <- list(c("x.1", "y"), c("x.1", "y"))
  expect_error(uniR2mx(bad.pooled, RAM=lavaan2RAM("y ~ x"), run=FALSE), "x\\.1.*x_1")
  expect_error(uniR2lavaan(bad.pooled, "y ~ x"), "x\\.1.*x_1")
  expect_error(uniR2lavaan(pooled, "y ~ b.1*x"), "b\\.1.*b_1")
  expect_error(bootuniR2("y ~ x", bad, n=100), "x\\.1.*x_1")
  expect_error(calEffSizes("y ~ b.1*x", n=100, Cov=cov), "b\\.1.*b_1")
  df <- Cor2DataFrame(studies, c(100, 100))
  expect_equal(df$obslabels, c("x", "y"))
  ram <- lavaan2RAM("y ~ x", obs.variables=c("x", "y"))
  bad.df <- df
  bad.df$obslabels[1] <- "x.1"
  expect_error(osmasem(RAM=ram, data=bad.df, run=FALSE), "x\\.1.*x_1")
  expect_error(osmasem2(RAM=ram, data=bad.df, run=FALSE), "x\\.1.*x_1")
  bad.df <- df
  bad.df$data[["v.1"]] <- 1
  expect_error(osmasem2(RAM=ram, data=bad.df, run=FALSE), "v\\.1.*v_1")
  expect_error(Cor2DataFrame(list(data=studies, n=c(100, 100), "mod.1"=c(0, 1))), "mod\\.1.*mod_1")
  expect_no_error(Cor2DataFrame(list(data=studies, n=c(100, 100), "mod.1"=c(0, 1)), append.vars=FALSE))
})

test_that("dotted group values are distinct from dotted column names", {
  ram <- lavaan2RAM("y ~ c(b1,b2)*x", ngroups=2)
  data <- data.frame(x=1:6, y=c(1, 3, 2, 4, 3, 5), site=rep(c("site.1", "site.2"), each=3))
  expect_no_error(sem(RAM=ram, data=data, group="site", run=FALSE))
})

test_that("consistently renamed variables retain the regression estimate", {
  data <- data.frame(x_1=seq(-2, 2, length.out=30))
  data$y_1 <- 1 + .7 * data$x_1 + sin(seq_len(30)) / 5
  ram <- lavaan2RAM("y_1 ~ b*x_1")
  fit <- sem(RAM=ram, data=data, silent=TRUE)
  expect_equal(unname(coef(fit)["b"]), unname(coef(lm(y_1 ~ x_1, data=data))[2]),
               tolerance=1e-5)
})

test_that("meta-analysis and direct matrix helpers check variable names", {
  data <- data.frame(y=1:4, v=rep(.1, 4), check.names=FALSE)
  data[["mod.1"]] <- 1:4
  expect_error(meta(y, v, data=data, run=FALSE), "mod\\.1.*mod_1")
  expect_error(reml(y, v, data=data, run=FALSE), "mod\\.1.*mod_1")
  Y <- matrix(1:4, ncol=1, dimnames=list(NULL, "effect.1"))
  V <- matrix(rep(.1, 4), ncol=1)
  expect_error(meta(Y, V, run=FALSE), "effect\\.1.*effect_1")
  A <- matrix("0", 2, 2, dimnames=list(c("x.1", "y"), c("x.1", "y")))
  S <- diag(2)
  expect_error(create.vechsR(A, S), "x\\.1.*x_1")
  expect_error(impliedR(Amatrix=A, Smatrix=S), "x\\.1.*x_1")
  expect_error(create.modMatrix(list(A=A, S=S, F=S), mod="mod"), "x\\.1.*x_1")
  expect_error(create.Tau2(no.var=1, Tau1.labels="tau.1"), "tau\\.1.*tau_1")
  expect_no_error(create.Tau2(no.var=2))
  expect_error(smdMES(c(1, 2), c(0, 0), A, S, 20, 20), "x\\.1.*x_1")
})
