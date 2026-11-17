# Cluster-sandwich inference for a single-level OpenMx ML model.
#
# The construction follows OpenMx::imxRobustSE(), originally written by
# Michael D. Hunter in OpenMx's R/MxRobustSE.R (Apache License 2.0): obtain
# derivatives of the -2 log likelihood for each observation and combine
# them with the model's information matrix. The extension here sums scores
# within clusters before forming the meat. Row derivatives are obtained by
# central differences of OpenMx's row likelihoods (2p+1 model evaluations
# in total) rather than one mxRun() per row; OpenMx::imxRowGradients() is
# used as a fallback. No OpenMx source code is copied.
# Source: https://github.com/OpenMx/OpenMx/blob/master/R/MxRobustSE.R

#' Cluster-robust covariance for a single-level OpenMx model
#'
#' Calculates a sandwich covariance from per-row derivatives of the log
#' likelihood, summed within independent clusters. The model is not
#' refitted. This function is for single-level raw-data maximum-likelihood
#' models, including flat multiple-group models; relational two-level
#' models use the separate cluster-score calculation in \code{vcov.mxsem()}.
#'
#' The method follows OpenMx's \code{imxRobustSE()} by Michael D. Hunter.
#' Row derivatives are computed by central differences of the row
#' likelihoods returned by OpenMx, which needs \eqn{2p+1} model evaluations
#' regardless of the number of rows; \code{OpenMx::imxRowGradients()} is
#' used if the row likelihoods cannot be evaluated this way. No OpenMx
#' source code is copied. OpenMx's original \code{R/MxRobustSE.R} is
#' distributed under the Apache License, Version 2.0.
#'
#' The meat is multiplied by \eqn{G/(G-1)}, where \eqn{G} is the number of
#' clusters. With one row per cluster, this is OpenMx's observation-level
#' factor \eqn{N/(N-1)}, so the result equals that of \code{imxRobustSE()}.
#' The bread is the inverse of the observed information. The result matches
#' \code{lavaan}'s cluster-robust standard errors of single-group models
#' (\code{se="robust.huber.white"} with \code{cluster=}) when \code{lavaan}
#' also uses the observed information (\code{information="observed"}).
#' \code{lavaan}'s default for \code{cluster=} is the expected information,
#' which gives the same result for complete data and saturated models but
#' can differ otherwise. With data missing at random, the expected
#' information can make the standard errors much too small, and the observed
#' information is the appropriate choice (Kenward & Molenberghs, 1998).
#'
#' Parameters estimated at a binding lower or upper bound have no valid
#' sandwich covariance. Their rows and columns are \code{NA}, with a
#' warning, and the covariance of the remaining parameters treats them as
#' fixed, while respecting equality constraints and binding inequality
#' constraints. A bound or an inequality constraint binds if the estimates
#' would move past it without it: the Newton step of -2LL under the other
#' binding constraints would violate it by more than 0.001 of its standard
#' error. An estimate that is merely close to a bound or to the boundary of
#' an inequality is not restricted by it and has the usual covariance. A
#' warning is also given when the OpenMx status code is neither 0 nor 1,
#' because the fit may not have converged.
#'
#' For models with \code{mxConstraint}s, the bread is not taken from OpenMx,
#' whose covariance for these models can be wrong. It is computed from the
#' Hessian of the Lagrangian, i.e., a numerical Hessian of -2LL without the
#' constraints minus the curvature of nonlinear constraints weighted by
#' their Lagrange multipliers, restricted to the null space of the binding
#' constraints (all equality constraints and the binding inequality
#' constraints). This equals the covariance of the same model with the
#' constraints substituted into it. A warning is
#' given if the gradient of -2LL is not orthogonal to these constraints,
#' because the fit may then not be at a constrained optimum.
#' @references OpenMx Project, \code{R/MxRobustSE.R}:
#' \url{https://github.com/OpenMx/OpenMx/blob/master/R/MxRobustSE.R}.
#'
#' Kenward, M. G., & Molenberghs, G. (1998). Likelihood based frequentist
#' inference when data are missing at random. \emph{Statistical Science},
#' \bold{13}(3), 236--247. \doi{10.1214/ss/1028905886}
#'
#' @param model A fitted, single-level \code{MxModel} with raw data and a
#' maximum-likelihood fit function.
#' @param cluster A vector of cluster identifiers in the same order as the
#' model's data rows, or a single character string naming a column in the
#' model's raw data. For a multiple-group model, a vector follows the order
#' of the contributing submodels, and a named column is taken from the data
#' of each of them. An ID that appears in more than one group is one
#' cluster.
#' @param details If \code{TRUE}, return a list containing the covariance,
#' standard errors, bread, meat (including the correction factor), cluster
#' count, and correction factor applied to the meat.
#' Otherwise return the standard errors.
#' @param smallSample If \code{TRUE}, also multiply the meat by
#' \code{(N-1)/(N-p)}, where \code{N} is the number of rows and \code{p}
#' the number of free parameters, so that the total factor is
#' \code{G/(G-1) * (N-1)/(N-p)}. This conventional correction is optional
#' and is not a general solution for inference with few clusters. \code{p}
#' counts all free parameters, including those determined by equality
#' constraints or estimated at a bound, so the factor is slightly larger
#' than with the number of parameters that are actually free to vary.
#' @return A named vector of standard errors, or a list when
#' \code{details=TRUE}. The \code{cov} element is the full covariance matrix.
#' @seealso \code{\link[OpenMx]{imxRobustSE}}, \code{\link[OpenMx]{imxRowGradients}},
#' \code{\link{sem}}
#' @examples
#' \donttest{
#' set.seed(1)
#' id <- rep(seq_len(20), each=5)
#' d <- data.frame(id=id, x=rnorm(100), y=rnorm(100))
#' ram <- lavaan2RAM("y ~ b*x", obs.variables=c("y", "x"))
#' fit <- sem(RAM=ram, data=d, cluster="id")
#' clusterRobustSE(fit$mx.fit, cluster="id")
#' vcov(fit, robust=TRUE)
#' }
clusterRobustSE <- function(model, cluster, details=FALSE,
                            smallSample=FALSE) {
  if (!is.logical(details) || length(details) != 1L || is.na(details) ||
      !is.logical(smallSample) || length(smallSample) != 1L ||
      is.na(smallSample))
    stop("'details' and 'smallSample' must be single logical values.")
  if (missing(cluster) || is.null(cluster))
    stop("'cluster' must be provided.")
  n <- .check.flat.ml(model)
  if (is.character(cluster) && length(cluster) == 1L && n > 1L) {
    ## A column name: its values in each contributing group, stacked in the
    ## order of the row scores (see .row.m2ll.gradients()).
    groups <- if (length(model@submodels))
      model@submodels[sub("\\.fitfunction$", "", model@fitfunction$groups)]
    else list(model)
    if (!all(vapply(groups, function(x) cluster %in% colnames(x@data@observed),
                    logical(1))))
      stop("The named cluster column is absent from the model data.")
    cluster <- unlist(lapply(groups, function(x)
      as.character(x@data@observed[[cluster]])), use.names=FALSE)
  }
  if (length(cluster) != n || anyNA(cluster))
    stop("'cluster' must contain one non-missing ID per model data row.")
  if (length(unique(cluster)) < 2L)
    stop("Cluster-robust covariance requires at least two clusters.")

  ## Log-likelihood scores, summed within clusters
  scores <- .row.m2ll.gradients(model) / -2
  cluster.scores <- rowsum(scores, group=as.factor(cluster), reorder=FALSE)
  g <- nrow(cluster.scores)
  p <- ncol(scores)
  correction <- g/(g-1)
  if (smallSample) {
    if (n <= p) stop("Small-sample correction requires more rows than parameters.")
    correction <- correction * ((n-1)/(n-p))
  }
  out <- .score.sandwich(model, cluster.scores, correction)
  if (details)
    list(SE=sqrt(diag(out$cov)), cov=out$cov, bread=out$bread,
         meat=out$meat, numClusters=g, smallSampleFactor=correction)
  else sqrt(diag(out$cov))
}

## Check that a fitted model is a flat raw-data ML model (single group or
## flat multiple groups) and return its total number of data rows.
.check.flat.ml <- function(model) {
  if (!methods::is(model, "MxModel") || is.null(model@output) ||
      !length(model@output$vcov))
    stop("'model' must be a fitted OpenMx model with a covariance matrix.")
  if (length(model@submodels)) {
    if (!methods::is(model@fitfunction, "MxFitFunctionMultigroup"))
      stop("Nested two-level models are unsupported by row gradients; use vcov() on the mxsem fit.")
    group.names <- sub("\\.fitfunction$", "", model@fitfunction$groups)
    groups <- model@submodels[group.names]
    if (!length(groups) || any(vapply(groups, is.null, logical(1))) ||
        any(vapply(groups, function(x) length(x@submodels) > 0L ||
                     is.null(x@data) || x@data$type != "raw" ||
                     !methods::is(x@fitfunction, "MxFitFunctionML"),
                   logical(1))))
      stop("All contributing groups must be flat raw-data ML models.")
    sum(vapply(groups, function(x) nrow(x@data@observed), integer(1)))
  } else {
    if (is.null(model@data) || model@data$type != "raw" ||
        !methods::is(model@fitfunction, "MxFitFunctionML"))
      stop("'model' must be a flat raw-data ML model.")
    nrow(model@data@observed)
  }
}

## Per-row derivatives of the -2 log likelihood (rows x parameters) for a
## fitted flat raw-data ML model, with rows stacked in the order of the
## contributing groups. Row likelihoods come from mxFitFunctionML(vector=TRUE)
## and are differentiated by central differences, using a one-sided
## difference when one side cannot be evaluated: 2p+1 model evaluations in
## total. Falls back to OpenMx::imxRowGradients() (one mxRun() per row) if
## the model cannot be evaluated this way.
.row.m2ll.gradients <- function(model) {
  theta <- coef(model)
  out <- tryCatch(.row.m2ll.differences(model, theta), error=function(e) NULL)
  if (is.null(out)) {
    out <- as.matrix(OpenMx::imxRowGradients(model, robustSE=TRUE))
    if (is.null(colnames(out))) colnames(out) <- names(theta)
    if (!setequal(colnames(out), names(theta)))
      stop("OpenMx gradient columns do not match model parameters.")
    out <- out[, names(theta), drop=FALSE]
  }
  if (!is.numeric(out) || any(!is.finite(out)))
    stop("Could not evaluate finite row gradients for this model.")
  out
}

.row.m2ll.differences <- function(model, theta) {
  strip <- function(m) {
    m@intervals <- list()
    m@compute <- NULL
    m@output <- list()
    m
  }
  if (length(model@submodels)) {
    groups <- sub("\\.fitfunction$", "", model@fitfunction$groups)
    evaluation <- strip(model)
    for (g in groups)
      evaluation <- mxModel(evaluation,
                            mxModel(strip(evaluation[[g]]),
                                    mxFitFunctionML(vector=TRUE)))
    fit.names <- groups
    n.rows <- vapply(groups, function(g) nrow(model[[g]]@data@observed),
                     integer(1))
  } else {
    inner <- mxModel(strip(model), mxFitFunctionML(vector=TRUE))
    evaluation <- mxModel(paste0(model@name, "_rowScores"), inner)
    fit.names <- model@name
    n.rows <- nrow(model@data@observed)
  }
  total <- paste0("-2*(", paste0("sum(log(", fit.names, ".fitfunction))",
                                  collapse=" + "), ")")
  evaluation <- mxModel(evaluation,
                        mxAlgebraFromString(total, name="rowScoresM2LL"),
                        mxFitFunctionAlgebra("rowScoresM2LL"))
  row.m2ll <- function(values) {
    m <- omxSetParameters(evaluation, labels=names(values), values=values)
    r <- mxRun(m, useOptimizer=FALSE, silent=TRUE, suppressWarnings=TRUE)
    out <- lapply(fit.names, function(x) -2*log(as.vector(r[[x]]$fitfunction$result)))
    if (!identical(lengths(out), unname(n.rows)))
      stop("Row likelihoods do not match the data rows.")
    unlist(out, use.names=FALSE)
  }
  base <- row.m2ll(theta)
  if (any(!is.finite(base))) stop("Row likelihoods are not finite at the estimates.")
  out <- matrix(0, length(base), length(theta),
                dimnames=list(NULL, names(theta)))
  for (j in seq_along(theta)) {
    step <- 1e-5 * max(1, abs(theta[j]))
    plus <- minus <- theta
    plus[j] <- theta[j] + step
    minus[j] <- theta[j] - step
    fp <- tryCatch(row.m2ll(plus), error=function(e) NA_real_)
    fm <- tryCatch(row.m2ll(minus), error=function(e) NA_real_)
    ok.p <- all(is.finite(fp))
    ok.m <- all(is.finite(fm))
    out[, j] <- if (ok.p && ok.m) (fp-fm)/(2*step) else if (ok.p)
      (fp-base)/step else if (ok.m) (base-fm)/step else
      stop("Could not evaluate row likelihoods near the estimates.")
  }
  out
}

## Free parameters estimated close to a lower or upper bound: a logical
## vector with a "side" attribute (-1 lower, 1 upper, 0 neither). Whether
## such a bound binds is decided by .binding.status().
.at.bound <- function(model, theta=coef(model)) {
  bound <- function(fetch) {
    b <- omxGetParameters(model, fetch=fetch)
    b[names(theta)]
  }
  near <- function(b) !is.na(b) & abs(theta - b) <= 1e-6 * pmax(1, abs(b))
  lower <- near(bound("lbound"))
  upper <- near(bound("ubound")) & !lower
  structure(setNames(lower | upper, names(theta)),
            side=ifelse(lower, -1, ifelse(upper, 1, 0)),
            bound=ifelse(lower, bound("lbound"), bound("ubound")))
}

## The mxConstraints of a model and its submodels as g(theta) = lhs - rhs,
## with the path to the (sub)model in which each is evaluated. "sign" is
## 1 for "<" and -1 for ">", so that sign*g <= 0 when the constraint holds,
## and 0 for "==".
.constraint.list <- function(model, path=character(0)) {
  out <- lapply(model@constraints, function(x) {
    op <- as.character(x@formula[[1]])
    list(path=path, expr=call("-", x@formula[[2]], x@formula[[3]]),
         equality=identical(op, "=="),
         sign=switch(op, "<"=1, ">"=-1, 0))
  })
  c(unname(out), unlist(lapply(names(model@submodels), function(k)
    .constraint.list(model@submodels[[k]], c(path, k))), recursive=FALSE))
}

## Values of all constraint functions at "theta"; the "equality" and "sign"
## attributes give, for each element, whether it comes from an equality
## constraint and its sign as in .constraint.list().
.constraint.values <- function(model, constraints, theta) {
  model <- omxSetParameters(model, labels=names(theta), values=theta)
  values <- lapply(constraints, function(x) {
    m <- model
    for (k in x$path) m <- m@submodels[[k]]
    as.vector(eval(bquote(mxEval(.(x$expr), m, compute=TRUE))))
  })
  out <- unlist(values)
  attr(out, "equality") <- rep(vapply(constraints, function(x) x$equality,
                                      logical(1)), lengths(values))
  attr(out, "sign") <- rep(vapply(constraints, function(x) x$sign,
                                  numeric(1)), lengths(values))
  out
}

## The constraints that may hold at the estimates: all equality constraints
## and the inequality constraints close to their boundary. Returns NULL if
## there are none, or a list with their Jacobian (rows: constraints,
## columns: "p"), a function giving their values at a parameter vector
## named as "p", the estimates, their values at the estimates, and their
## "equality" flags and "sign"s. Whether a close inequality binds is
## decided by .binding.status(). OpenMx's own constraintJacobian is not used
## because it reports zero rows for inequality constraints that are
## satisfied, whether or not they bind.
.active.constraints <- function(model, p) {
  constraints <- .constraint.list(model)
  if (!length(constraints)) return(NULL)
  ## Unlabeled free parameters (e.g., in wls() models) cannot be set by name.
  ## They are labeled here; the order of the parameters is unchanged.
  original <- names(omxGetParameters(model))
  model <- omxNameAnonymousParameters(model)[[1]]
  theta <- omxGetParameters(model)
  internal <- names(theta)[match(p, original)]
  g <- .constraint.values(model, constraints, theta)
  if (any(!is.finite(g)))
    stop("Cannot calculate constrained covariance: constraint values are not finite.", call.=FALSE)
  active <- attr(g, "equality") | abs(g) <= 10*.feasibility.tolerance
  if (!any(active)) return(NULL)
  values <- function(x) {
    full <- theta
    full[internal] <- x
    .constraint.values(model, constraints, full)[active]
  }
  x <- setNames(theta[internal], p)
  jacobian <- vapply(seq_along(p), function(j) {
    h <- 1e-6*max(1, abs(x[[j]]))
    up <- down <- x
    up[j] <- up[j] + h
    down[j] <- down[j] - h
    (values(up) - values(down))/(2*h)
  }, numeric(sum(active)))
  list(jacobian=matrix(jacobian, nrow=sum(active), dimnames=list(NULL, p)),
       values=values, estimate=x, value=as.vector(g)[active],
       equality=attr(g, "equality")[active], sign=attr(g, "sign")[active])
}

## Keep the rows "keep" of the constraints from .active.constraints().
.subset.constraints <- function(active, keep) {
  if (is.null(active) || !any(keep)) return(NULL)
  values <- active$values
  list(jacobian=active$jacobian[keep, , drop=FALSE],
       values=function(x) values(x)[keep], estimate=active$estimate,
       value=active$value[keep], equality=active$equality[keep],
       sign=active$sign[keep])
}

.active.constraint.jacobian <- function(model, p) {
  .binding.status(model, p)$active$jacobian
}

## Hessian of the fit function (-2LL) at the estimates, with its gradient as
## the "gradient" attribute. With mxConstraints, OpenMx's reported Hessian
## is not that of the fit function (it can be indefinite), so it is
## recomputed numerically without the constraints.
.fit.hessian <- function(model, p) {
  if (!.has.constraints(model)) {
    hessian <- model@output$hessian
    if (is.null(hessian) || !all(p %in% rownames(hessian)) ||
        !all(p %in% colnames(hessian)))
      stop("Robust covariance with parameters at a bound requires the model's Hessian.")
    gradient <- model@output$gradient
    return(structure(hessian[p, p, drop=FALSE],
                     gradient=if (all(p %in% names(gradient))) gradient[p]))
  }
  drop <- function(m) {
    m@constraints <- list()
    m@submodels <- lapply(m@submodels, drop)
    m
  }
  m <- drop(model)
  m@intervals <- list()
  m@compute <- mxComputeSequence(list(
    mxComputeOnce("fitfunction", "fit"),
    mxComputeNumericDeriv(checkGradient=FALSE),
    mxComputeReportDeriv()))
  m <- mxRun(m, silent=TRUE, suppressWarnings=TRUE, unsafe=TRUE)
  hessian <- m@output$hessian
  if (is.null(hessian) || !all(p %in% rownames(hessian)) ||
      !all(p %in% colnames(hessian)) || any(!is.finite(hessian[p,p,drop=FALSE])))
    stop("Cannot calculate constrained covariance: a finite, named Hessian is unavailable.", call.=FALSE)
  structure(hessian[p, p, drop=FALSE], gradient=m@output$gradient[p])
}

## Which of the bounds and inequality constraints close to the estimates
## bind. Distance alone cannot tell: an estimate can be within a tolerance
## of a bound or an inequality without being restricted by it. Starting
## from the equality constraints, the Newton step of -2LL in the directions
## that the constraints taken as binding leave free predicts where the
## estimates would move without the others; negative curvature is replaced
## by its absolute value so that the step is downhill. A close bound or inequality
## binds if that step would violate it by more than 0.001 of its standard
## error (from the same constrained covariance); the most violated one is
## added and the step recomputed until none is violated. One that the step
## does not violate leaves the estimates where they are, so they have the
## usual covariance, as if it were absent. If the gradient is unavailable or
## the step cannot be computed, every close bound and inequality is taken
## as binding.
## Returns the Hessian (with a "gradient" attribute, or NULL when not
## needed), the binding constraints as from .active.constraints(), and a
## logical vector of parameters at a binding bound.
.binding.status <- function(model, p) {
  theta <- coef(model)[p]
  near <- .at.bound(model, theta)
  constrained <- .has.constraints(model)
  bounded <- setNames(rep(FALSE, length(p)), p)
  if (!any(near) && !constrained)
    return(list(hessian=NULL, active=NULL, bounded=bounded))
  hessian <- .fit.hessian(model, p)
  active <- if (constrained) .active.constraints(model, p)
  close <- if (is.null(active)) logical(0) else !active$equality
  if (!any(near) && !any(close))
    return(list(hessian=hessian, active=active, bounded=bounded))

  side <- attr(near, "side")
  limit <- attr(near, "bound")
  gradient <- attr(hessian, "gradient")
  ## Candidates: close inequalities (rows of the Jacobian), then close bounds
  rows <- which(close)
  params <- which(near)
  n <- length(rows) + length(params)
  coefficient <- rbind(if (length(rows)) active$sign[rows]*active$jacobian[rows, , drop=FALSE],
                       if (length(params)) side[params]*diag(length(p))[params, , drop=FALSE])
  slack <- c(if (length(rows)) active$sign[rows]*active$value[rows],
             if (length(params)) side[params]*(theta[params] - limit[params]))
  equality <- if (is.null(active)) NULL
              else active$jacobian[active$equality, , drop=FALSE]
  binding <- rep(FALSE, n)
  fallback <- function() list(hessian=hessian,
                              active=active,
                              bounded=setNames(near, p))
  if (is.null(gradient) || any(!is.finite(gradient))) {
    warning("The fit gradient is unavailable or non-finite; binding constraints could not be assessed.", call.=FALSE)
    return(fallback())
  }
  h <- matrix(hessian, length(p), length(p))
  for (iteration in seq_len(n)) {
    constraint <- rbind(equality, coefficient[binding, , drop=FALSE])
    z <- if (is.null(constraint) || !nrow(constraint)) diag(length(p)) else {
      decomposition <- svd(t(constraint), nu=length(p))
      rank <- sum(decomposition$d > 1e-8 * max(decomposition$d))
      if (rank >= length(p)) break
      decomposition$u[, seq.int(rank+1L, length(p)), drop=FALSE]
    }
    ## Absolute eigenvalues: where -2LL is not convex (e.g., a variance at
    ## 0), the step still points downhill.
    eigen.h <- eigen(crossprod(z, h %*% z), symmetric=TRUE)
    lambda <- abs(eigen.h$values)
    if (!length(lambda) || !all(is.finite(lambda)) || max(lambda) <= 0)
      return(fallback())
    lambda <- pmax(lambda, 1e-8*max(lambda))
    reduced <- eigen.h$vectors %*% (t(eigen.h$vectors)/lambda)
    step <- -drop(z %*% reduced %*% crossprod(z, gradient))
    covariance <- 2 * z %*% reduced %*% t(z)
    predicted <- slack + drop(coefficient %*% step)
    se <- sqrt(pmax(rowSums((coefficient %*% covariance) * coefficient), 0))
    violation <- ifelse(binding, -Inf, predicted - pmax(1e-3*se, 1e-8))
    if (max(violation) <= 0) break
    binding[which.max(violation/pmax(se, 1e-12))] <- TRUE
  }
  keep <- if (is.null(active)) logical(0) else active$equality
  if (length(rows)) keep[rows] <- binding[seq_along(rows)]
  if (length(params)) bounded[params] <- binding[length(rows) + seq_along(params)]
  list(hessian=hessian, active=.subset.constraints(active, keep),
       bounded=bounded)
}

## Curvature of the active constraints weighted by their Lagrange
## multipliers, sum_k lambda_k d2 c_k / d theta d theta', on the scale of
## -2LL. At a constrained minimum of the fit function F, grad F = C' lambda,
## where C stacks the active constraint Jacobian and the rows e_k of
## parameters at a bound. The Hessian of the Lagrangian, H - this
## curvature, is the Hessian of F with the constraints substituted into the
## model, so its inverse on the tangent space is the covariance of the
## constrained estimates. Linear constraints have no curvature. Second
## derivatives are taken only for parameters on which the constraints
## depend, by central differences of the multiplier-weighted constraints.
.constraint.curvature <- function(active, jacobian, gradient, p) {
  out <- matrix(0, length(p), length(p), dimnames=list(p, p))
  if (is.null(gradient) || anyNA(gradient[p])) {
    warning("The gradient of the fit function is unavailable, so the curvature ",
            "of nonlinear constraints is ignored in the standard errors.",
            call.=FALSE)
    return(out)
  }
  gradient <- gradient[p]
  if (any(!is.finite(gradient)) || any(!is.finite(jacobian)))
    stop("Cannot calculate constrained covariance: the gradient or constraint Jacobian is non-finite.", call.=FALSE)
  decomposition <- svd(t(jacobian))
  keep <- decomposition$d > 1e-8*max(decomposition$d)
  lambda <- decomposition$v[, keep, drop=FALSE] %*%
    (crossprod(decomposition$u[, keep, drop=FALSE], gradient)/decomposition$d[keep])
  residual <- gradient - drop(t(jacobian) %*% lambda)
  if (sqrt(sum(residual^2)) > 0.01*max(1, sqrt(sum(gradient^2))))
    warning("The gradient of the fit function is not orthogonal to the active ",
            "constraints, so the fit may not be at a constrained optimum and ",
            "its standard errors may be unreliable. You are advised to rerun it, ",
            "e.g., with rerun().", call.=FALSE)
  lambda <- lambda[seq_len(nrow(active$jacobian))]
  if (all(abs(lambda) < 1e-12)) return(out)
  x <- active$estimate
  constraint.values <- function(v) {
    out <- active$values(v)
    if (any(!is.finite(out)))
      stop("Cannot calculate constrained covariance: constraint values are non-finite near the estimate.", call.=FALSE)
    out
  }
  weighted <- function(v) sum(lambda*constraint.values(v))
  step <- 1e-4*pmax(1, abs(x))
  ## Parameters on which the constraints depend (checked by a finite change
  ## so that, e.g., a^2 at a = 0 is not missed).
  depends <- vapply(seq_along(x), function(j) {
    v <- x
    v[j] <- v[j] + 1e-3*max(1, abs(x[j]))
    any(abs(constraint.values(v) - constraint.values(x)) > 1e-12)
  }, logical(1))
  idx <- which(depends)
  f0 <- weighted(x)
  shift <- function(i, si, j=NULL, sj=0) {
    v <- x
    v[i] <- v[i] + si*step[i]
    if (!is.null(j)) v[j] <- v[j] + sj*step[j]
    weighted(v)
  }
  for (a in seq_along(idx)) {
    i <- idx[a]
    out[i, i] <- (shift(i, 1) - 2*f0 + shift(i, -1))/step[i]^2
    for (b in seq_len(a - 1L)) {
      j <- idx[b]
      out[i, j] <- out[j, i] <-
        (shift(i, 1, j, 1) - shift(i, 1, j, -1) - shift(i, -1, j, 1) +
           shift(i, -1, j, -1))/(4*step[i]*step[j])
    }
  }
  out
}

## Covariance of the estimates (the inverse information on the vcov() scale)
## under the constraints that bind at the estimates:
## V = 2 Z (Z' H_L Z)^{-1} Z', where H_L is the Hessian of the Lagrangian of
## -2LL (the Hessian of -2LL minus the multiplier-weighted curvature of
## nonlinear constraints) and Z spans the null space of the binding
## constraint Jacobian. Parameters at a binding bound ("bounded") are
## treated as fixed by adding them to the constraints; their rows and
## columns are NA. "status" is from .binding.status().
.constrained.vcov <- function(model, p=names(omxGetParameters(model)),
                              bounded=rep(FALSE, length(p)),
                              status=.binding.status(model, p)) {
  hessian <- status$hessian
  if (is.null(hessian)) hessian <- .fit.hessian(model, p)
  active <- status$active
  jacobian <- rbind(active$jacobian,
                    diag(length(p))[bounded, , drop=FALSE])
  out <- matrix(NA_real_, length(p), length(p), dimnames=list(p, p))
  if (is.null(jacobian) || !nrow(jacobian)) {
    z <- diag(length(p))
  } else {
    decomposition <- svd(t(jacobian), nu=length(p))
    rank <- sum(decomposition$d > 1e-8 * max(decomposition$d))
    if (rank >= length(p)) {
      ## Equality constraints can determine every parameter locally. Their
      ## tangent space is empty, so their sampling covariance is zero;
      ## retain NA only for parameters treated as fixed at binding bounds.
      out[] <- 0
      out[bounded, ] <- NA_real_
      out[, bounded] <- NA_real_
      return(out)
    }
    z <- decomposition$u[, seq.int(rank+1L, length(p)), drop=FALSE]
    if (!is.null(active))
      hessian <- hessian - .constraint.curvature(active, jacobian,
                                                 attr(hessian, "gradient"), p)
  }
  hessian <- matrix(hessian, length(p), length(p))
  reduced <- crossprod(z, hessian %*% z)
  if (any(!is.finite(reduced)))
    stop("Cannot calculate constrained covariance: the reduced Hessian is non-finite.", call.=FALSE)
  inverse <- tryCatch(solve(reduced), error=function(e)
    stop("Cannot calculate constrained covariance: the reduced Hessian is singular; the model may not be identified.", call.=FALSE))
  out[] <- 2 * z %*% inverse %*% t(z)
  out[bounded, ] <- NA_real_
  out[, bounded] <- NA_real_
  out
}

## Model-based covariance of the estimates: vcov(model), except for models
## with mxConstraints, whose OpenMx covariance is replaced as above.
.model.vcov <- function(model) {
  if (!.has.constraints(model)) return(vcov(model))
  .constrained.vcov(model)
}

## Bread of the sandwich: the model-based covariance, with parameters at a
## binding bound treated as fixed (their rows and columns are NA).
.robust.bread <- function(model, p, status=.binding.status(model, p)) {
  if (!any(status$bounded) && !.has.constraints(model))
    return(vcov(model)[p, p, drop=FALSE])
  .constrained.vcov(model, p, status$bounded, status)
}

## Sandwich covariance from log-likelihood scores (independent units x
## parameters), using the constraint-adjusted bread from .robust.bread().
## Parameters at a bound are treated as fixed: their rows and columns are NA.
.score.sandwich <- function(model, scores, factor=1) {
  status <- model@output$status$code
  if (length(status) && !(status %in% c(0, 1)))
    warning("OpenMx status code is ", status, ", so the model may not have ",
            "converged and the robust standard errors may be unreliable. ",
            "You are advised to rerun it, e.g., with rerun(), and check that ",
            "the status code is 0 or 1.", call.=FALSE)
  p <- names(coef(model))
  status <- .binding.status(model, p)
  bounded <- status$bounded
  keep <- !bounded
  bread <- .robust.bread(model, p, status)
  meat <- covariance <- matrix(NA_real_, length(p), length(p),
                               dimnames=list(p, p))
  if (any(keep)) {
    b <- bread[keep, keep, drop=FALSE]
    m <- crossprod(scores[, p[keep], drop=FALSE]) * factor
    v <- b %*% m %*% b
    meat[keep, keep] <- m
    covariance[keep, keep] <- (v + t(v))/2
  }
  if (any(bounded))
    warning("Robust standard errors are unavailable for parameter(s) at a bound (",
            paste(p[bounded], collapse=", "), "); they are NA, and the ",
            "remaining parameters treat them as fixed.", call.=FALSE)
  list(cov=covariance, bread=bread, meat=meat)
}
