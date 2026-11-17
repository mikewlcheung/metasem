test_that("supported lavaanify options are represented across model structures", {
  cases <- list(
    list("y ~ x1 + x2", list(auto.cov.x=TRUE), c("x1", "~~", "x2")),
    list("y ~ x\nz ~ x", list(auto.cov.y=TRUE), c("y", "~~", "z")),
    list("f1 =~ x1+x2+x3\nf2 =~ x4+x5+x6", list(auto.cov.lv.x=TRUE), c("f1", "~~", "f2")),
    list("f =~ x1+x2+x3", list(auto.fix.first=TRUE, std.lv=FALSE), c("f", "=~", "x1")),
    list("f =~ x1", list(auto.fix.single=TRUE), c("x1", "~~", "x1"))
  )
  for (case in cases) for (kind in c("single", "multigroup", "twolevel", "combined")) {
    syntax <- case[[1]]
    if (kind %in% c("twolevel", "combined"))
      syntax <- paste("level: 1", syntax, "level: 2", syntax, sep="\n")
    groups <- if (kind %in% c("multigroup", "combined")) 2 else 1
    if (kind == "combined") syntax <- paste("group: A", syntax, "group: B", syntax, sep="\n")
    args <- c(list(model=syntax, ngroups=groups), case[[2]])
    reference.args <- c(args, list(auto.var=TRUE, fixed.x=FALSE))
    if (is.null(reference.args$std.lv)) reference.args$std.lv <- TRUE
    table <- suppressWarnings(do.call(lavaan::lavaanify, reference.args))
    ram <- suppressWarnings(do.call(lavaan2RAM, args))
    key <- case[[3]]
    rows <- table[table$lhs==key[1] & table$op==key[2] & table$rhs==key[3], , drop=FALSE]
    expect_gt(nrow(rows), 0)
    for (i in seq_len(nrow(rows))) {
      row <- rows[i, ]; part <- ram
      if (groups > 1) part <- part[[if (kind == "combined") match(as.character(row$group), c("A", "B")) else as.integer(row$group)]]
      if (kind %in% c("twolevel", "combined")) part <- part[[if (row$level == 1) "within" else "between"]]
      value <- if (row$op == "=~") part$A[row$rhs, row$lhs] else part$S[row$lhs, row$rhs]
      expect_identical(grepl("*", value, fixed=TRUE), row$free > 0, info=paste(kind, key[1]))
      if (row$free == 0 && !is.na(row$ustart)) expect_equal(as.numeric(value), row$ustart)
    }
  }
})

test_that("unsupported options and syntax fail before translation", {
  opts <- c("orthogonal", "warn", "debug", "effect.coding", "group.equal",
            "correlation", "int.ov.free", "marker.int.zero", "auto.th",
            "fixed.x", "fixed_x", "auto_var", "std_lv", "as.data.frame")
  for (opt in opts) {
    args <- setNames(list(TRUE), opt)
    expect_error(do.call(lavaan2RAM, c(list(model="y ~ x"), args)), "Unsupported argument", info=opt)
  }
  expect_error(lavaan2RAM("y ~ x", auto=TRUE), "'auto' is not supported")
  expect_error(lavaan2RAM("y ~ x", auto=FALSE), "'auto' is not supported")
  forward <- function(...) lavaan2RAM("y ~ x", ...)
  expect_error(forward(auto=TRUE), "'auto' is not supported")
  expect_error(lavaan2RAM("y ~ x", std=FALSE), "exact documented")
  expect_error(lavaan2RAM("y ~ x", auto.cov.x=NA), "must be TRUE or FALSE")
  expect_error(lavaan2RAM("f <~ x1 + x2"), "Unsupported lavaan operator")
  expect_error(lavaan2RAM("y | t1"), "Unsupported lavaan operator")
  expect_error(lavaan2RAM('efa("block")*f =~ x1+x2+x3'), "Exploratory factor")
  expect_error(lavaan2RAM("y ~ b*x\n.p1. == b"), "automatic .p labels|generated equality")
})

test_that("explicit constraints, factor anchors and mean defaults are preserved", {
  expect_warning(lavaan2RAM("f =~ x1+x2+x3", std.lv=FALSE), "scales may be unidentified")
  expect_no_warning(lavaan2RAM("f =~ x1+x2+x3", std.lv=FALSE, auto.fix.first=TRUE))
  expect_no_warning(lavaan2RAM("f =~ 1*x1+x2+x3", std.lv=FALSE))
  ram <- lavaan2RAM("y ~ b*x", constraints="b == 0.5", M.notation="_mu", M.start=0.3)
  expect_true(length(ram$mxalgebras) > 0)
  expect_match(ram$M[1, "y"], "0.3\\*y_mu")
  multi <- suppressWarnings(lavaan2RAM("y ~ c(b,b)*x", ngroups=2, M.notation="_mu", M.start=0.3))
  expect_match(multi[[1]]$M[1, "y"], "0.3\\*y_mu")
  expect_match(multi[[2]]$M[1, "y"], "0.3\\*y_mu")
})

test_that("supported automatic covariance and named constraints fit like lavaan", {
  set.seed(731)
  n <- 400
  x <- rnorm(n); e <- rnorm(n); y <- 0.4*x + e
  z <- -0.2*x + 0.4*e + rnorm(n)
  data <- data.frame(x=x, y=y, z=z)
  for (constrained in c(FALSE, TRUE)) {
    syntax <- if (constrained) "y ~ b*x" else "y ~ x\nz ~ x"
    options <- if (constrained) list(constraints="b == 0.5") else list(auto.cov.y=TRUE)
    ram <- do.call(lavaan2RAM, c(list(model=syntax), options))
    fit <- sem(RAM=ram, data=data)
    ref <- do.call(lavaan::sem, c(list(model=syntax, data=data, meanstructure=TRUE,
                                      fixed.x=FALSE, std.lv=TRUE), options))
    expect_true(fit$mx.fit$output$status$code %in% c(0, 1))
    expect_true(lavaan::lavInspect(ref, "converged"))
    expect_equal(as.numeric(fit$mx.fit$output$fit), -2*as.numeric(logLik(ref)), tolerance=1e-5)
    if (!constrained) {
      pars <- lavaan::parameterEstimates(ref)
      row <- pars$lhs == "y" & pars$op == "~~" & pars$rhs == "z"
      expect_equal(unname(coef(fit)["yWITHz"]), pars$est[row], tolerance=1e-5)
      expect_equal(unname(sqrt(vcov(fit)["yWITHz", "yWITHz"])), pars$se[row], tolerance=1e-5)
      expect_true(all(eigen(vcov(fit), symmetric=TRUE, only.values=TRUE)$values > 0))
      expect_true(lavaan::lavInspect(ref, "post.check"))
    }
  }
})

test_that("correlation models constrain observed total variances with fixed residuals", {
  vars <- c("x", "y")
  C <- matrix(c(1, sqrt(0.5), sqrt(0.5), 1), 2, dimnames=list(vars, vars))
  d <- Cor2DataFrame(rep(list(C), 6), n=rep(100, 6))
  ram <- lavaan2RAM("y ~ b*x\nx ~~ 1*x\ny ~~ 0.5*y", obs.variables=vars)
  model <- suppressWarnings(osmasem2(RAM=ram, data=d, RE.type.Sigma="Zero", run=FALSE))
  expect_true("unitVariance1" %in% names(model$mx.fit@constraints))
  fit <- suppressWarnings(osmasem2(RAM=ram, data=d, RE.type.Sigma="Zero"))
  expect_true(fit$mx.fit$output$status$code %in% c(0, 1))
  implied <- OpenMx::mxEvalByName("expSigma", fit$mx.fit)
  expect_equal(unname(diag(implied)), c(1, 1), tolerance=1e-5)
  expect_equal(unname(implied[1, 2]), sqrt(0.5), tolerance=1e-5)
  impossible <- lavaan2RAM("y ~ 0.4*x\nx ~~ 1*x\ny ~~ 0.5*y", obs.variables=vars)
  expect_error(suppressWarnings(osmasem2(RAM=impossible, data=d, RE.type.Sigma="Zero", run=FALSE)),
               "unit-variance restrictions are incompatible")
})

test_that("correlation-only osmasem rejects covariance data and warns about overrides", {
  vars <- c("x", "y")
  C <- matrix(c(1, .3, .3, 1), 2, dimnames=list(vars, vars))
  corr <- Cor2DataFrame(rep(list(C), 4), n=rep(100, 4))
  cov <- Cor2DataFrame(rep(list(C), 4), n=rep(100, 4), cor.analysis=FALSE)
  ram <- lavaan2RAM("y ~ b*x\nx ~~ 1*x\ny ~~ 0.5*y", obs.variables=vars)
  expect_error(osmasem(RAM=ram, data=cov, run=FALSE), "requires correlation data")
  expect_error(osmasem(RAM=ram, data=corr, cor.analysis=FALSE, run=FALSE), "correlation analysis only")
  expect_warning(suppressWarnings(checkRAM(ram$A, ram$S)), NA)
  expect_warning(suppressWarnings(create.vechsR(ram$A, ram$S, ram$F), classes="nonexistent"), "fixed values in S0")
  expect_warning(suppressWarnings(wls(C, matrix(0.01,1,1), 400, RAM=ram, run=FALSE), classes="nonexistent"), "fixed values in Smatrix")
})

test_that("WLS warns when an explicit unit variance constraint is omitted", {
  vars <- c("x", "y", "z")
  C <- matrix(c(1,.3,.2,.3,1,.1,.2,.1,1), 3, dimnames=list(vars,vars))
  ram <- lavaan2RAM("y ~ b*x\nz ~ c*x\nx ~~ 1*x\ny ~~ 0.5*y", obs.variables=vars)
  expect_warning(wls(C, diag(3)/100, 400, RAM=ram, diag.constraints=TRUE, run=FALSE),
                 "does not impose an explicit unit-variance constraint")
})

test_that("existing variance options and shared group labels remain available", {
  free <- lavaan2RAM("y ~ x", auto.var=TRUE)
  fixed <- lavaan2RAM("y ~ x", auto.var=FALSE)
  expect_true(as.mxMatrix(free$S)$free["y", "y"])
  expect_false(as.mxMatrix(fixed$S)$free["y", "y"])
  latent <- lavaan2RAM("f =~ x1+x2+x3", std.lv=TRUE)
  expect_equal(as.numeric(latent$S["f", "f"]), 1)
  shared <- suppressWarnings(lavaan2RAM("y ~ c(b,b)*x", ngroups=2))
  expect_equal(as.mxMatrix(shared[[1]]$A)$labels["y", "x"], "b")
  expect_equal(as.mxMatrix(shared[[2]]$A)$labels["y", "x"], "b")
})

test_that("forwarded argument names are checked without breaking positional calls", {
  wrapper <- function(...) lavaan2RAM("f =~ x1+x2+x3", ...)
  nested <- function(...) wrapper(...)
  expect_error(wrapper(std=FALSE), "Unsupported argument\\(s\\): std")
  expect_error(nested(auto.fix.f=TRUE), "auto.fix.f")
  expect_error(do.call(wrapper, list(std=FALSE)), "Unsupported argument\\(s\\): std")
  expect_error(wrapper(std=stop("must not be evaluated")), "Unsupported argument\\(s\\): std")
  expect_error(nested(auto=TRUE), "'auto' is not supported")
  expect_no_warning(nested(std.lv=FALSE, auto.fix.first=TRUE))
  positional <- lavaan2RAM("y ~ x", NULL, "on")
  omitted <- lavaan2RAM("y ~ x", , "on")
  expect_equal(positional$A, omitted$A)
  expect_equal(as.mxMatrix(positional$A)$labels["y", "x"], "yonx")
})

test_that("group counts are valid and character counts warn", {
  for (value in list(0, -1, 1.2, NA_real_, Inf, NaN, numeric(), c(1,2),
                     TRUE, "zero", "1.2", "0", "2e100")) {
    expect_error(lavaan2RAM("y ~ x", ngroups=value), "single positive integer")
  }
  expect_equal(length(lavaan2RAM("y ~ x", ngroups=2)), 2)
  expect_warning(result <- lavaan2RAM("y ~ x", ngroups="2"), "deprecated")
  expect_equal(length(result), 2)
})

test_that("simple nonzero named constraints count as scale anchors", {
  syntax <- "f =~ a*x1 + x2 + x3 + x4"
  expect_no_warning(ram <- lavaan2RAM(paste(syntax, "a == 1", sep="\n"), std.lv=FALSE))
  expect_no_warning(lavaan2RAM(syntax, std.lv=FALSE, constraints="1 == a"))
  expect_no_warning(lavaan2RAM("f =~ x1+x2+x3\nf ~~ v*f\nv == 2", std.lv=FALSE))
  expect_warning(lavaan2RAM(paste(syntax, "a == 0", sep="\n"), std.lv=FALSE), "scales may be unidentified")
  expect_warning(lavaan2RAM(paste(syntax, "a == b", sep="\n"), std.lv=FALSE), "scales may be unidentified")
  set.seed(84)
  n <- 400; factor <- rnorm(n)
  d <- as.data.frame(sapply(1:4, function(i) (0.7+i*.1)*factor+rnorm(n)))
  names(d) <- paste0("x", 1:4)
  fit <- sem(RAM=ram, data=d)
  ref <- lavaan::sem(syntax, data=d, std.lv=FALSE, auto.fix.first=FALSE,
                     constraints="a == 1", meanstructure=TRUE)
  expect_true(fit$mx.fit$output$status$code %in% c(0,1))
  expect_true(lavaan::lavInspect(ref, "converged"))
  expect_true(lavaan::lavInspect(ref, "post.check"))
  expect_equal(as.numeric(fit$mx.fit$output$fit), -2*as.numeric(logLik(ref)), tolerance=1e-5)
  expect_equal(unname(coef(fit)["x2ONf"]), unname(lavaan::coef(ref)["f=~x2"]), tolerance=1e-5)
})

test_that("named scale anchors are checked within each group and level", {
  multi <- "f =~ c(a,b)*x1 + x2 + x3\na == 1"
  expect_warning(lavaan2RAM(multi, ngroups=2, std.lv=FALSE), "at least one group/level")
  expect_no_warning(lavaan2RAM(paste(multi, "b == 1", sep="\n"), ngroups=2, std.lv=FALSE))
  levels <- paste("level: 1", "f =~ a*x1+x2+x3", "level: 2", "f =~ b*x1+x2+x3", sep="\n")
  expect_warning(lavaan2RAM(paste(levels, "a == 1", sep="\n"), std.lv=FALSE), "at least one group/level")
  expect_no_warning(lavaan2RAM(paste(levels, "a == 1", "b == 1", sep="\n"), std.lv=FALSE))
})
