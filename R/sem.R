#' Fit a structural equation model using OpenMx
#' 
#' It fits a structural equation model by creating a mxModel from a RAM object.
#' Variable names, raw data column names, and ordinary parameter labels must
#' not contain periods. Rename consistently, e.g., \code{x.1} to \code{x_1},
#' and check for duplicate names. Definition-variable references such as
#' \code{data.x_1} and group or cluster values containing periods are allowed.
#' 
#' 
#' @aliases sem create.mxModel
#' @usage
#' sem(model.name = "sem", RAM = NULL, data = NULL,
#'     Cov = NULL, means = NULL, numObs, intervals.type = c("z",
#'         "LB"), startvalues = NULL, lbound = NULL, ubound = NULL,
#'     replace.constraints = FALSE, mxModel.Args = NULL,
#'     run = TRUE, silent = TRUE, group = NULL, cluster = NULL,
#'     ...)
#' create.mxModel(model.name = "sem", RAM = NULL, data = NULL,
#'     Cov = NULL, means = NULL, numObs, intervals.type = c("z",
#'         "LB"), startvalues = NULL, replace.constraints = FALSE,
#'     mxModel.Args = NULL, run = TRUE, silent = TRUE,
#'     ...)
#' @param model.name A string for the model name in
#' \code{\link[OpenMx]{mxModel}}.
#' @param RAM A RAM object including a list of matrices of the model returned
#' from \code{\link[metaSEM]{lavaan2RAM}}.
#' @param data A data frame or matrix of data.
#' @param Cov If \code{data=NULL}, a covariance matrix for a single-level
#' model, a list of matrices (one per group) for a multiple-group model, or
#' a list of matrices (one per cluster) for a two-level model. For a combined
#' group and two-level model, supply a list of groups, each containing a list
#' of cluster covariance matrices. Matrices must have named variables and
#' use the usual unbiased covariance denominator \code{n-1}. For two-level
#' summary input, covariance reconstruction is checked in each variable's
#' scale; incompatible matrices are rejected. When group or cluster names are
#' supplied on both \code{Cov} and \code{means} or \code{numObs}, entries are
#' matched by name; names must be complete, unique, and consistent. Unnamed
#' entries are matched by position.
#' @param means Means corresponding to \code{Cov}: a named vector for a
#' single-level model, a list of vectors for multiple groups or clusters,
#' or a nested list for combined group and two-level models. Means are
#' required for two-level summary input. For single-level and multiple-group
#' input they are optional; omitted means remove the mean structure.
#' @param numObs Sample sizes corresponding to \code{Cov}: a scalar for a
#' single-level model, a vector for multiple groups or clusters, or a list
#' of vectors for combined group and two-level models.
#' @param intervals.type Either \code{z} (default if missing) or \code{LB}. If
#' it is \code{z}, it calculates the 95\% confidence intervals (CIs) based on
#' the estimated standard error. If it is \code{LB}, it calculates the 95\%
#' likelihood-based CIs on the parameter estimates.
#' @param startvalues A list of named starting values of the free parameters,
#' e.g., list(a=1, b=2)
#' @param lbound A list of lower bound of the free parameters. If it is not
#' provided, all free parameters are assumed \code{NA}. Bounds from
#' inequalities in the lavaan syntax, recorded by
#' \code{\link[metaSEM]{lavaan2RAM}}, are also applied; a bound given here
#' for the same parameter takes precedence. Bounds must be named and refer
#' to model parameters. Bounds on substituted expressions remain constraints;
#' definition-variable substitutions with bounds are unsupported.
#' @param ubound A list of upper bound of the free parameters. If it is not
#' provided, all free parameters are assumed \code{NA}. Bounds from the
#' lavaan syntax are handled as for \code{lbound}.
#' @param replace.constraints Logical. If \code{TRUE}, the parameters on the
#' left hand side will be replaced by the constraints on the right hand side.
#' That is, the parameters on the left hand side are no longer parameters in
#' the model.
#' @param mxModel.Args A list of arguments passed to
#' \code{\link[OpenMx]{mxModel}}.
#' @param run Logical. If \code{FALSE}, only return the mx model without
#' running the analysis.
#' @param silent Logical. Passed to \code{\link[OpenMx]{mxRun}}. If
#' \code{FALSE}, warnings and errors during the \code{mxTryHard()} retries
#' are also printed as messages. A fit that ends with an OpenMx status code
#' other than 0 or 1 is retried once with \code{\link[OpenMx]{mxTryHard}}
#' from its estimates; a warning is given if the status code is still
#' neither 0 nor 1.
#' Combined group and cluster summary inputs require distinct column names.
#' @param group A character string naming the column in \code{data}
#' identifying group membership. Required when \code{RAM} is
#' a multiple-group RAM object, i.e., created by \code{\link[metaSEM]{lavaan2RAM}}
#' with \code{ngroups > 1}. \code{data} must be a single data frame in this
#' case; it is split internally by the unique values of \code{data[[group]]}
#' in \emph{first-appearance order in the data} (not sorted, and not a
#' factor's own \code{levels()} order, even for a factor column -- this
#' matches how \code{lavaan} itself numbers groups), which must correspond,
#' in order, to the groups used when the model was specified (e.g.,
#' \code{lavaan}'s \code{c(...)} syntax numbers groups the same way, and so
#' do \code{group:} blocks, in their order in the syntax rather than by
#' their labels). A warning is given when the labels of \code{group:} blocks
#' are the values of \code{data[[group]]} in a different order. A
#' group value that is not itself a legal \code{\link[OpenMx]{mxModel}} name
#' (e.g. purely numeric, or containing a character such as \code{.} or
#' \code{-}) is sanitized into one internally; the resulting submodel names
#' are available via \code{names(fit$mx.fit@submodels)}, and a mapping back
#' to the original data values via \code{fit$group.map}, so
#' \code{\link{plot.mxsem}}'s own \code{group} argument can still be called
#' with the original, unsanitized group value. An unused factor level in
#' \code{data[[group]]} (one that never actually occurs in the data) is
#' simply never treated as a group in the first place (this is what
#' first-appearance order, above, naturally does); if that leaves fewer
#' groups than \code{RAM} expects, it is reported as a group-count
#' mismatch. \code{Cov}/\code{numObs} input is supported with one covariance
#' matrix and sample size per group. \code{replace.constraints=TRUE} is supported
#' (substituting a constraint into whichever group's \code{A}/\code{S}/
#' \code{M} matrices contain it, including a constraint that spans
#' multiple groups, e.g. constraining one group's variance to be a
#' multiple of another's); a definition variable (\code{data.*}) inside
#' such a constraint resolves against whichever group's own data it is
#' used in, since (unlike two-level's between level) every group already
#' carries its own full data slice.
#' @param cluster A character string naming the column in \code{data}
#' identifying cluster membership. It is required when \code{RAM} is a
#' two-level RAM object, i.e., created by \code{\link[metaSEM]{lavaan2RAM}}
#' from \code{lavaan}'s two-level (\code{level: 1} / \code{level: 2})
#' syntax; \code{data} must then be a single data frame of individual-level
#' (level-1) rows. In a single-level or flat multiple-group model, it records
#' the independent units for cluster-robust standard errors in
#' \code{summary(fit, robust=TRUE)} and \code{vcov(fit, robust=TRUE)}; the
#' fitted likelihood and estimates do not change, and raw data are required.
#' See Details.
#' @param \dots Further arguments will be passed to either
#' \code{\link[OpenMx]{mxRun}} or \code{\link[OpenMx]{mxTryHard}}
#' @details
#' \strong{Two-level models.} A manifest variable appearing at both levels
#' is internally split into within and between components. Consequently, a
#' level-1 regression uses only the predictor's within component, as in
#' \code{lavaan}. Two-level summary input requires separate means,
#' covariance matrices, and sample sizes for each cluster; a pooled
#' covariance matrix cannot identify within- and between-cluster variation.
#' For complete Gaussian data these summaries reproduce the raw-data
#' likelihood. Combined group and two-level models use nested
#' \code{group:}/\code{level:} lavaan syntax and both \code{group} and
#' \code{cluster} arguments. \code{replace.constraints=TRUE} is supported
#' for two-level models (substituting a constraint into whichever of the
#' within-/between-level \code{A}/\code{S}/\code{M} matrices contains it,
#' including a constraint that spans both levels); a definition variable
#' (\code{data.*}) may be used in such a constraint at either level if its
#' between-level values are constant within each cluster. Cluster-level
#' observed covariates likewise require a constant value within each
#' cluster. The reported number of observations and the information
#' criteria use the number of level-1 rows with at least one observed
#' level-1 variable, as in \code{lavaan}. A cluster with a measured
#' between-level variable still contributes, even if all its level-1 data
#' are missing.
#'
#' \strong{Constraints.} Models with constraints between parameters
#' (\code{==}, \code{<}, or \code{>}) are fitted with an OpenMx
#' feasibility tolerance of 1e-6 unless one has been set with
#' \code{mxOption()}, and their model-based standard errors are computed
#' from the Hessian of the Lagrangian of -2LL and the constraints that hold
#' at the estimates instead of being taken from OpenMx. They equal those of
#' the same model with the constraints substituted
#' (\code{replace.constraints=TRUE}); for nonlinear constraints they can
#' differ from \code{lavaan}'s, which omit the curvature of the constraints.
#'
#' \strong{Robust standard errors.} With a single-level model and no
#' \code{cluster}, robust covariance uses the observation-level sandwich of
#' OpenMx's \code{imxRobustSE()}, whose meat is multiplied by \eqn{N/(N-1)}.
#' With \code{cluster}, it sums the row scores within clusters before
#' forming the sandwich and multiplies the meat by \eqn{G/(G-1)}, where
#' \eqn{G} is the number of clusters, as \code{lavaan} does for
#' single-group models; this allows within-cluster dependence. With one row
#' per cluster it equals the robust covariance without \code{cluster}. For
#' flat multiple-group models, the same cluster ID in different groups
#' denotes the same independent unit; use distinct IDs when groups contain
#' distinct clusters. The bread uses the observed information, so it
#' matches \code{lavaan} with \code{information="observed"};
#' \code{lavaan}'s default for \code{cluster=}, the expected information,
#' can give much smaller standard errors with data missing at random (see
#' \code{\link{clusterRobustSE}}). The optional \code{smallSample=TRUE} in
#' \code{summary()} or \code{vcov()} also multiplies the meat by
#' \eqn{(N-1)/(N-p)}, where \eqn{N} is the number of rows and \eqn{p} the
#' number of free parameters (including those determined by constraints or
#' at a bound), for a total factor of \eqn{G/(G-1) \cdot (N-1)/(N-p)}. With
#' few independent clusters, \code{summary(fit, robust=TRUE, t.df=)} uses a
#' \eqn{t} reference distribution instead of \eqn{z}; see
#' \code{\link[metaSEM]{summary}} for how to choose its degrees of freedom.
#' Row scores are computed by numerical differentiation of the row
#' likelihoods, so their cost hardly depends on the number of rows.
#' Parameters estimated at a binding bound have no valid robust covariance:
#' their robust standard errors are \code{NA}, with a warning, and the other
#' parameters treat them as fixed. A bound binds if the estimate would move
#' past it without the bound; an estimate merely close to a bound has the
#' usual robust standard error. A defined (\code{:=}) parameter has an
#' \code{NA} robust standard error if it depends on any bounded parameter.
#' Robust covariance also warns when the OpenMx status code is neither 0
#' nor 1, because the fit may not have converged.
#'
#' For two-level fits, robust covariance uses each cluster's complete
#' relational-model likelihood score, without a \eqn{G/(G-1)} factor, as in
#' \code{lavaan}'s two-level models. Each score needs one OpenMx run that
#' differentiates the cluster's likelihood numerically, which takes several
#' likelihood evaluations per free parameter, so computing time grows with
#' both the number of clusters and the number of parameters. In combined
#' group and two-level models, clusters are nested within groups: each
#' group's clusters are formed from that group's rows only, so a cluster ID
#' used in two groups denotes two different clusters, and a unit spanning
#' groups cannot be represented. \code{sem()} warns when an ID appears in
#' more than one group, and robust covariance stops with an error in that
#' case. Cluster IDs generated from summary input are prefixed by the
#' group's position (e.g., \code{"2.5"} for cluster \code{5} in the second
#' group) when they would otherwise repeat. Robust summaries require
#' \code{intervals.type="z"}; likelihood-based intervals are not robust.
#' @return An object of class \code{mxsem}
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @keywords utilities
#' @examples
#' 
#' \donttest{
#' ## Generate data
#' set.seed(100)
#' n <- 100
#' x <- rnorm(n)
#' y <- 0.5*x + rnorm(n, mean=0, sd=sqrt(1-0.5^2))
#' my.df <- data.frame(y=y, x=x)
#' 
#' ## A regression model
#' model <- "y ~ x   # Regress y on x
#'           y ~ 1   # Intercept of y
#'           x ~ 1   # Mean of x"
#' plot(model)
#' 
#' RAM <- lavaan2RAM(model, obs.variables=c("y", "x"))
#' 
#' my.fit <- sem(RAM=RAM, data=my.df)
#' summary(my.fit)
#' 
#' ## A meta-analysis
#' model <- "f =~ 1*yi
#'           f ~ mu*1          ## Average effect
#'           f ~~ tau2*f       ## Heterogeneity variance
#'           yi ~~ data.vi*yi  ## Known sampling variance"
#' plot(model)
#' 
#' ## Do not standardize the latent variable (f): std.lv=FALSE 
#' RAM <- lavaan2RAM(model, obs.variables="yi", std.lv=FALSE)
#' 
#' ## Use likelihood-based CI
#' my.fit <- sem(RAM=RAM, data=Hox02, intervals="LB")
#' summary(my.fit)
#'
#' ## A 2-group regression with a shared slope but group-specific intercepts
#' x1 <- rnorm(100); y1 <- 2 + 0.5*x1 + rnorm(100)
#' x2 <- rnorm(100); y2 <- -1 + 0.5*x2 + rnorm(100)
#' my.df2 <- data.frame(y=c(y1, y2), x=c(x1, x2),
#'                      g=rep(c("g1", "g2"), each=100))
#'
#' model <- "y ~ c(b,b)*x + c(a1,a2)*1"
#' RAM <- lavaan2RAM(model, ngroups=2)
#' my.fit <- sem(RAM=RAM, data=my.df2, group="g")
#' summary(my.fit)
#'
#' ## The same 2-group regression, extended with replace.constraints=TRUE:
#' ## group 2's residual variance is constrained to be exactly TWICE group
#' ## 1's, and group 1's own residual variance is reparameterised through
#' ## exp() (guaranteeing positivity). "resid1" never appears as a free
#' ## parameter itself -- lavaanify() would otherwise reject "a0" as
#' ## undeclared, and an ordinary (non-replaced) "==" constraint cannot
#' ## introduce a new free parameter this way; only replace.constraints=
#' ## TRUE can. The constraint chain is resolved across groups, so
#' ## "resid2" ends up expressed directly in terms of "a0" (not "resid1").
#' model <- "y ~ c(b,b)*x + c(a1,a2)*1
#'           y ~~ c(resid1,resid2)*y
#'           resid1 == exp(a0)
#'           resid2 == 2*resid1"
#' RAM <- lavaan2RAM(model, ngroups=2)
#' my.fit <- sem(RAM=RAM, data=my.df2, group="g", replace.constraints=TRUE)
#' summary(my.fit)
#'
#' ## A two-level (random-intercept) regression model
#' set.seed(123)
#' nClusters <- 60; nPerCluster <- 6
#' clusterID <- rep(seq_len(nClusters), each=nPerCluster)
#' u <- rep(rnorm(nClusters, 0, 0.7), each=nPerCluster)
#' x <- rnorm(nClusters*nPerCluster)
#' y <- 1.5 + 0.6*x + u + rnorm(nClusters*nPerCluster)
#' my.df3 <- data.frame(clusterID=clusterID, x=x, y=y)
#'
#' model <- "level: 1
#'             y ~ b1*x
#'             y ~~ residW*y
#'           level: 2
#'             y ~ 1
#'             y ~~ residB*y"
#' RAM <- lavaan2RAM(model)
#' my.fit <- sem(RAM=RAM, data=my.df3, cluster="clusterID")
#' summary(my.fit)
#'
#' ## A 2-group CFA: test whether factor loadings are invariant across groups
#' set.seed(103)
#' simGroup <- function(n, loadings) {
#'   f <- rnorm(n)
#'   Y <- sapply(loadings, function(l) l*f + rnorm(n, sd=0.5))
#'   colnames(Y) <- paste0("y", seq_along(loadings))
#'   data.frame(Y)
#' }
#' loadings <- c(1, 0.8, 1.3, 0.9)
#' my.df4 <- rbind(cbind(simGroup(250, loadings), g="g1"),
#'                cbind(simGroup(250, loadings), g="g2"))
#'
#' ## Metric invariance: loadings constrained equal across groups.
#' ## std.lv=FALSE: the factor's scale is already identified by fixing
#' ## y1's loading at 1 (c(1,1)*y1), so its variance ("f ~~ c(vf1,vf2)*f")
#' ## is left free and group-specific instead of also being fixed at 1 --
#' ## lavaan2RAM()'s default std.lv=TRUE would otherwise silently discard
#' ## that "f ~~ ..." line, fixing the variance at 1 in both groups with
#' ## no error (see lavaan2RAM()'s "std.lv" argument for details).
#' model.invariant <- "f =~ c(1,1)*y1 + c(l2,l2)*y2 + c(l3,l3)*y3 + c(l4,l4)*y4
#'                     f ~~ c(vf1,vf2)*f
#'                     y1 ~ c(m1a,m1b)*1
#'                     y2 ~ c(m2a,m2b)*1
#'                     y3 ~ c(m3a,m3b)*1
#'                     y4 ~ c(m4a,m4b)*1"
#' fit.invariant <- sem(RAM=lavaan2RAM(model.invariant, ngroups=2, std.lv=FALSE),
#'                      data=my.df4, group="g")
#'
#' ## Configural model: loadings free in each group
#' model.free <- "f =~ c(1,1)*y1 + c(l2a,l2b)*y2 + c(l3a,l3b)*y3 + c(l4a,l4b)*y4
#'                f ~~ c(vf1,vf2)*f
#'                y1 ~ c(m1a,m1b)*1
#'                y2 ~ c(m2a,m2b)*1
#'                y3 ~ c(m3a,m3b)*1
#'                y4 ~ c(m4a,m4b)*1"
#' fit.free <- sem(RAM=lavaan2RAM(model.free, ngroups=2, std.lv=FALSE),
#'                 data=my.df4, group="g")
#'
#' summary(fit.invariant)
#' anova(fit.free, fit.invariant)   # likelihood-ratio test of metric invariance
#'
#' ## Group difference in the factor variance, as a derived quantity
#' model.diff <- paste(model.free, "\ndiff_v := vf2 - vf1")
#' fit.diff <- sem(RAM=lavaan2RAM(model.diff, ngroups=2, std.lv=FALSE),
#'                data=my.df4, group="g")
#' summary(fit.diff)   # see "Mxalgebras" in the printed output for diff_v
#'
#' ## A two-level CFA: a factor at the within level and another at the
#' ## between level, with the same three indicators loading on both
#' set.seed(104)
#' nClusters <- 100; nPerCluster <- 8; N <- nClusters*nPerCluster
#' clusterID <- rep(seq_len(nClusters), each=nPerCluster)
#' fb <- rep(rnorm(nClusters, 0, 1), each=nPerCluster)
#' fw <- rnorm(N, 0, 1)
#' ly2 <- 0.8; ly3 <- 1.2
#' y1 <- fb + fw + rnorm(N, 0, 0.5)
#' y2 <- ly2*fb + ly2*fw + rnorm(N, 0, 0.5)
#' y3 <- ly3*fb + ly3*fw + rnorm(N, 0, 0.5)
#' my.df5 <- data.frame(clusterID=clusterID, y1=y1, y2=y2, y3=y3)
#'
#' model <- "level: 1
#'             fw =~ y1 + y2 + y3
#'           level: 2
#'             fb =~ y1 + y2 + y3"
#' RAM <- lavaan2RAM(model)
#' my.fit <- sem(RAM=RAM, data=my.df5, cluster="clusterID")
#' summary(my.fit)
#'
#' ## A two-level (multilevel) meta-analysis: studies contributing several
#' ## dependent effect sizes each, separating within- from between-study
#' ## heterogeneity while accounting for each effect size's own known
#' ## sampling variance
#' set.seed(106)
#' nStudies <- 60
#' esPerStudy <- sample(2:4, nStudies, replace=TRUE)
#' studyID <- rep(seq_len(nStudies), esPerStudy)
#' N <- length(studyID)
#' mu <- 0.35; tau2.between <- 0.03; tau2.within <- 0.015
#' studyEffect <- rep(rnorm(nStudies, mu, sqrt(tau2.between)), esPerStudy)
#' vi <- runif(N, 0.02, 0.06)
#' yi <- rnorm(N, studyEffect + rnorm(N, 0, sqrt(tau2.within)), sqrt(vi))
#' my.df6 <- data.frame(studyID=studyID, yi=yi, vi=vi)
#'
#' ## fw's own variance is the WITHIN-study heterogeneity; yi's own residual
#' ## (after removing fw) is fixed at the KNOWN sampling variance; the
#' ## BETWEEN-study heterogeneity is yi's own variance at level 2
#' model <- "level: 1
#'             fw =~ 1*yi
#'             fw ~~ tauW*fw
#'             yi ~~ data.vi*yi
#'           level: 2
#'             yi ~ mu*1
#'             yi ~~ tauB*yi"
#' RAM <- lavaan2RAM(model, std.lv=FALSE)
#' my.fit <- sem(RAM=RAM, data=my.df6, cluster="studyID")
#' summary(my.fit)
#'
#' ## A two-level location-scale meta-analysis: replace.constraints=TRUE
#' ## lets the WITHIN-study heterogeneity (tauW, on fw) be a function of an
#' ## effect-size-level moderator xi, while yi's own residual stays fixed
#' ## at its KNOWN sampling variance (data.vi), exactly as in the previous
#' ## example. data.xi is a definition variable, so its per-row value is
#' ## looked up from the data instead of estimated. "tauW" never appears
#' ## as a free parameter itself here -- lavaanify() would otherwise reject
#' ## a0/a1 as undeclared, and an ordinary (non-replaced) "==" constraint
#' ## cannot introduce new free parameters this way; only
#' ## replace.constraints=TRUE can.
#' set.seed(107)
#' xi <- rnorm(N)
#' a0 <- -4; a1 <- 0.8
#' tauWi <- exp(a0 + a1*xi)
#' yi2 <- rnorm(N, studyEffect + rnorm(N, 0, sqrt(tauWi)), sqrt(vi))
#' my.df7 <- data.frame(studyID=studyID, yi=yi2, vi=vi, xi=xi)
#'
#' model <- "level: 1
#'             fw =~ 1*yi
#'             fw ~~ tauW*fw
#'             yi ~~ data.vi*yi
#'           level: 2
#'             yi ~ mu*1
#'             yi ~~ tauB*yi
#'           tauW == exp(a0 + a1*data.xi)"
#' RAM <- lavaan2RAM(model, std.lv=FALSE)
#' my.fit <- sem(RAM=RAM, data=my.df7, cluster="studyID",
#'               replace.constraints=TRUE)
#' summary(my.fit)
#' }
#'
sem <- function(model.name="sem", RAM=NULL, data=NULL, Cov=NULL,
                means=NULL, numObs, intervals.type=c("z", "LB"),
                startvalues=NULL, lbound=NULL, ubound=NULL,
                replace.constraints=FALSE,
                mxModel.Args=NULL, run=TRUE, silent=TRUE, group=NULL,
                cluster=NULL, ...) {

  .check.ram.names(RAM)
  .check.no.period(colnames(data), "Data column names")
  .check.no.period(c(group, cluster), "Group and cluster column names")
  .check.cov.names(Cov)
  .check.mean.names(means)
  .check.parameter.names(c(names(startvalues), names(lbound), names(ubound)))

  intervals.type <- match.arg(intervals.type)

  ## Bounds from inequalities in the lavaan syntax (see lavaan2RAM());
  ## bounds given in 'lbound'/'ubound' take precedence.
  lbound <- .merge.bounds(attr(RAM, "lbound"), lbound)
  ubound <- .merge.bounds(attr(RAM, "ubound"), ubound)
  .check.parameter.bounds(RAM, mxModel.Args, lbound, ubound)
  .check.replacement.parameters(RAM, replace.constraints, mxModel.Args)

  if (inherits(RAM, "RAM_multigroup_twolevel")) {
    if (is.null(data)) {
      if (!is.list(Cov) || length(Cov) != length(RAM) ||
          !is.list(means) || length(means) != length(RAM) ||
          missing(numObs) || !is.list(numObs) ||
          length(numObs) != length(RAM))
        stop("Combined group and level summary input requires nested 'Cov', 'means', and 'numObs' lists, one per group and then one per cluster.")
      if (is.null(group) || is.null(cluster))
        stop("'group' and 'cluster' must name the generated ID columns for combined summary input.")
      if (identical(group, cluster))
        stop("'group' and 'cluster' must name different generated columns for combined summary input.", call.=FALSE)
      means <- .align.summary.entries(Cov, means, "means", "group")
      numObs <- .align.summary.entries(Cov, numObs, "numObs", "group")
      values <- names(Cov)
      if (is.null(values) || any(!nzchar(values))) values <- as.character(seq_along(Cov))
      pieces <- lapply(seq_along(Cov), function(i) {
        x <- .expand.cluster.summaries(Cov[[i]], means[[i]], numObs[[i]], cluster)
        x[[group]] <- values[i]
        x
      })
      ## Clusters are nested within groups; prefix generated IDs with the
      ## group's position when they would otherwise repeat across groups.
      ## A position has no period, so "<i>.<id>" cannot collide across
      ## groups, unlike a group name such as "a.b".
      ids <- lapply(pieces, function(x) unique(x[[cluster]]))
      if (anyDuplicated(unlist(ids)))
        pieces <- lapply(seq_along(pieces), function(i) {
          x <- pieces[[i]]
          x[[cluster]] <- paste(i, x[[cluster]], sep=".")
          x
        })
      data <- do.call(rbind, pieces)
    }
    return(.sem.multigroup.twolevel(model.name=model.name, RAM=RAM,
                                    data=data, group=group, cluster=cluster,
                                    intervals.type=intervals.type,
                                    startvalues=startvalues, lbound=lbound,
                                    ubound=ubound,
                                    replace.constraints=replace.constraints,
                                    mxModel.Args=mxModel.Args, run=run,
                                    silent=silent, ...))
  }

  ## Multiple-group RAM (from lavaan2RAM(..., ngroups>1)) and two-level RAM
  ## (from lavaan2RAM() on "level: 1"/"level: 2" syntax): dispatch each to
  ## a dedicated builder rather than threading extra code paths through the
  ## single-group/single-level logic below (kept untouched to avoid
  ## regressions).
  if (.is.RAM.multigroup(RAM)) {
    ## A plain multiple-group RAM has no between-level component.
    return(.sem.multigroup(model.name=model.name, RAM=RAM, data=data,
                           group=group, cluster=cluster, Cov=Cov, means=means,
                           numObs=if (missing(numObs)) NULL else numObs,
                           intervals.type=intervals.type,
                           startvalues=startvalues, lbound=lbound, ubound=ubound,
                           replace.constraints=replace.constraints,
                           mxModel.Args=mxModel.Args, run=run, silent=silent, ...))
  }
  if (.is.RAM.twolevel(RAM)) {
    if (!is.null(group)) {
      stop("'group' has no effect for a two-level RAM object (class ",
          "\"RAM_twolevel\"); combined two-level and multiple-group ",
          "models require nested group:/level: RAM syntax.")
    }
    if (is.null(data) && !is.null(Cov)) {
      if (missing(numObs))
        stop("Two-level summary input requires per-cluster 'numObs'.")
      data <- .expand.cluster.summaries(Cov, means, numObs, cluster)
    }
    return(.sem.twolevel(model.name=model.name, RAM=RAM, data=data,
                         cluster=cluster, intervals.type=intervals.type,
                         startvalues=startvalues, lbound=lbound, ubound=ubound,
                         replace.constraints=replace.constraints,
                         mxModel.Args=mxModel.Args, run=run, silent=silent, ...))
  }

  if (!is.null(group))
    stop("'group' requires a multiple-group RAM object.")
  if (!is.null(cluster)) {
    if (is.null(data) || !is.data.frame(data) ||
        !is.character(cluster) || length(cluster) != 1L ||
        !(cluster %in% names(data)) || anyNA(data[[cluster]]))
      stop("For a single-level model, 'cluster' must name a complete column in raw data.")
  }

  Amatrix <- as.mxMatrix(RAM$A, name="Amatrix")
  Smatrix <- as.mxMatrix(RAM$S, name="Smatrix")
  Fmatrix <- as.mxMatrix(RAM$F, name="Fmatrix")
  Mmatrix <- as.mxMatrix(RAM$M, name="Mmatrix")

  ## Some basic checking in RAM
  checkRAM(Amatrix, Smatrix, cor.analysis=FALSE)

  ## Extract all observed and latent variable names
  var.names <- colnames(Fmatrix$values)

  ## Without raw data
  if (is.null(data))  {
    ## Without means
    if (is.null(means)) {
      mx.data <- mxData(observed=Cov, type="cov", numObs=numObs)
      expFun <- mxExpectationRAM(A="Amatrix", S="Smatrix",
                                 F="Fmatrix", dimnames=var.names)
    } else {
      ## With means
      mx.data <- mxData(observed=Cov, type="cov", means=means, numObs=numObs)
      expFun <- mxExpectationRAM(A="Amatrix", S="Smatrix", F="Fmatrix",
                                 M="Mmatrix", dimnames=var.names)
    }
  } else {
    ## With raw data
    mx.data <- mxData(observed=data, type="raw")
    expFun <- mxExpectationRAM(A="Amatrix", S="Smatrix", F="Fmatrix",
                               M="Mmatrix", dimnames=var.names)
  }    

  ## Create an incomplete model, which will be used to store other mx objects.
  mx.model <- mxModel(model.name, mx.data, expFun, mxFitFunctionML())

  ## Note. startvalues may overwrite the starting values in RAM
  ## Collate the starting values from RAM and add them to startvalues
  para.labels <- c(Amatrix$labels[Amatrix$free], Smatrix$labels[Smatrix$free],
                   Mmatrix$labels[Mmatrix$free])
  para.values <- c(Amatrix$values[Amatrix$free], Smatrix$values[Smatrix$free],
                   Mmatrix$values[Mmatrix$free])
  ## Name the starting values with names, which is consistent with the startvalues
  names(para.values) <- para.labels
  para.values <- as.list(para.values)

  ## Prepare startvalues
  if (is.null(startvalues)) {## Note. startvalues may overwrite the starting values in RAM
    startvalues <- para.values
  } else {
    ## Remove starting values from para.values if they are overlapped with startvalues    
    para.values[names(para.values) %in% names(startvalues)] <- NULL
    startvalues <- c(startvalues, para.values)

    ## Replace the startvalues in Amatrix, Smatrix, and Mmatrix
    for (i in seq_along(startvalues)) {
      Amatrix$values[Amatrix$labels==names(startvalues)[i]] <- startvalues[[i]]
      Smatrix$values[Smatrix$labels==names(startvalues)[i]] <- startvalues[[i]]
      Mmatrix$values[Mmatrix$labels==names(startvalues)[i]] <- startvalues[[i]]
    }
  }

  ## Add lbound and ubound
  Amatrix <- .addbound(Amatrix, lbound=lbound, ubound=ubound)
  Smatrix <- .addbound(Smatrix, lbound=lbound, ubound=ubound)
  Mmatrix <- .addbound(Mmatrix, lbound=lbound, ubound=ubound)
  
  ## Extract a local copy for ease of reference
  ## Remove starting values for ease of matching
  A <- as.symMatrix(RAM$A)
  S <- as.symMatrix(RAM$S)
  M <- as.symMatrix(RAM$M)
  mxalgebras <- RAM$mxalgebras
  
  ## Replace the parameters on the left-hand side of "==" constraints by
  ## their right-hand sides if replace.constraints=TRUE (see
  ## .select.substitutions()).
  if (isTRUE(replace.constraints)) {
    selection <- .select.substitutions(mxalgebras, para.labels, lbound, ubound)
    A <- .apply.replacements(A, selection$replacements)
    S <- .apply.replacements(S, selection$replacements)
    M <- .apply.replacements(M, selection$replacements)
    mxalgebras <- selection$remaining
  }

  ## Check whether there are replacements
  ## Remove the starting values before comparisons
  if (all(A==as.symMatrix(RAM$A))) {
    mx.model <- mxModel(mx.model, Amatrix)
  } else if (length(all.vars(parse(text=A))) == 0) {
    ## All free parameters in A were replaced with numeric constants: use
    ## a fixed MxMatrix instead of an mxAlgebra (as.mxAlgebra fails on
    ## zero free variables due to paste0("0*", character(0)) yielding
    ## "0*"). as.mxMatrix() cannot evaluate a non-trivial constant
    ## EXPRESSION -- it only recognises bare numbers or "value*label"
    ## strings, silently treating anything else (e.g. "exp(0.4)", which
    ## contains no "*") as an UNLABELLED FREE parameter with no starting
    ## value instead of a fixed one. Confirmed by testing: without
    ## evaluating by hand first, "sigma == exp(0.4)" left sigma free and
    ## estimated at the data's own MLE, not fixed at exp(0.4) -- same fix
    ## already applied to .twolevel.replace.one() and .build.group.mxmodel().
    values <- matrix(vapply(A, function(z) eval(parse(text=z)), numeric(1)),
                     nrow=nrow(A), ncol=ncol(A), dimnames=dimnames(A))
    mx.model <- mxModel(mx.model, as.mxMatrix(values, name="Amatrix"))
  } else {
    A <- as.mxAlgebra(A, startvalues=startvalues, lbound=lbound, ubound=ubound,
                      name="Amatrix")
    mx.model <- mxModel(mx.model, A$mxalgebra, A$parameters, A$list)
  }

  if (all(S==as.symMatrix(RAM$S))) {
    mx.model <- mxModel(mx.model, Smatrix)
  } else if (length(all.vars(parse(text=S))) == 0) {
    ## Same guard as for A above -- evaluate any constant EXPRESSION by
    ## hand first, rather than handing raw expression text to
    ## as.mxMatrix().
    values <- matrix(vapply(S, function(z) eval(parse(text=z)), numeric(1)),
                     nrow=nrow(S), ncol=ncol(S), dimnames=dimnames(S))
    mx.model <- mxModel(mx.model, as.mxMatrix(values, name="Smatrix"))
  } else {
    S <- as.mxAlgebra(S, startvalues=startvalues, lbound=lbound, ubound=ubound,
                      name="Smatrix")
    mx.model <- mxModel(mx.model, S$mxalgebra, S$parameters, S$list)
  }
         
  ## Create an identity matrix from the no. of columens of Fmatrix,
  ## including all latent and observed variables
  Id <- as.mxMatrix(diag(ncol(Fmatrix$values)), name="Id")
  Id_A <- mxAlgebra(solve(Id - Amatrix), name="Id_A")
  
  ## Note. expCov and expMean are NOT used in the fit function.
  ## They are included so the implied structures include the latent variables.
  ## It may be useful for future applications.
  expCov <- mxAlgebra(Id_A %&% Smatrix, name="expCov")   
  expMean <- mxAlgebra(Mmatrix %*% t(Id_A), name="expMean")
  
  ## Add the mean structure only if there are means
  if (!is.null(data) | !is.null(means)) {

    if (all(M==as.symMatrix(RAM$M))) {
      mx.model <- mxModel(mx.model, Mmatrix)
    } else if (length(all.vars(parse(text=M))) == 0) {
      ## Same guard as for A and S above -- evaluate any constant
      ## EXPRESSION by hand first, rather than handing raw expression text
      ## to as.mxMatrix().
      values <- matrix(vapply(M, function(z) eval(parse(text=z)), numeric(1)),
                       nrow=nrow(M), ncol=ncol(M), dimnames=dimnames(M))
      mx.model <- mxModel(mx.model, as.mxMatrix(values, name="Mmatrix"))
    } else {
      M <- as.mxAlgebra(M, startvalues=startvalues, lbound=lbound,
                        ubound=ubound, name="Mmatrix")
      mx.model <- mxModel(mx.model, M$mxalgebra, M$parameters, M$list)
    }
          
    mx.model <- mxModel(mx.model, Fmatrix, Id, Id_A, expCov, expMean)
  } else {
    ## No mean structure
    mx.model <- mxModel(mx.model, Fmatrix, Id, Id_A, expCov)
  }
  
  ## Add additional arguments to mxModel
  if (!is.null(mxModel.Args)) {
    for (i in seq_along(mxModel.Args)) {
      mx.model <- mxModel(mx.model, mxModel.Args[[i]])
    }
  }

  ## New parameter labels including those in constraints
  new.para.labels <- unique(c(A, S, if (!is.null(data) || !is.null(means)) M))
  ## Get the variable names
  new.para.labels <- all.vars(parse(text=new.para.labels))
  ## Drop the definition variables. "." in grepl()'s regex means "any
  ## character" -- a literal-looking label such as "database" (data + b +
  ## ase) also matches "data.", wrongly dropping an ordinary free
  ## parameter from mxCI() here. Match only the literal "data." prefix.
  new.para.labels <- new.para.labels[!startsWith(new.para.labels, "data.")]
  ## mxCI() rejects an empty reference vector outright ("'reference'
  ## argument must be a character vector") rather than treating it as
  ## "no CIs requested" -- confirmed by testing. Reachable whenever every
  ## free parameter has been fixed (e.g. replace.constraints=TRUE fixing
  ## every "==" target at a literal constant leaves nothing free at all),
  ## which is a valid, if unusual, model to fit; skip the call entirely
  ## rather than let it fail before the model is even run.
  if (length(new.para.labels) > 0) {
    mx.model <- mxModel(mx.model, mxCI(new.para.labels))
  }

  ## A list of mxalgebras required SE or CI
  mxalgebras.ci <- NULL
  
  ## Add mxAlgebra and mxConstraint from RAM$mxalgebra
  if (!is.null(mxalgebras)) {
    for (i in seq_along(mxalgebras)) {
      mx.model <- mxModel(mx.model, mxalgebras[[i]])
      ## Name of the mxalgebra
      name.mxalgebra <- names(mxalgebras)[i]
      ## Include the mxalgebras, but not the mxConstraints, in mxCI. The
      ## class is checked: a defined parameter may be named, e.g.,
      ## "constraint1_effect".
      if (!inherits(mxalgebras[[i]], "MxConstraint")) {
        mx.model <- mxModel(mx.model, mxCI(c(name.mxalgebra)))
        mxalgebras.ci <- c(mxalgebras.ci, name.mxalgebra)
      }
    }
  }
    
  ## Initialise here so they are available regardless of whether run=TRUE
  warning.msg <- error.msg <- NULL

  if (run) {
    fallback <- .run.mxTryHard.fallback(mx.model, intervals.type, silent=silent, ...)
    mx.fit <- fallback$mx.fit
    warning.msg <- fallback$warning.msg
    error.msg <- fallback$error.msg
  } else {
    mx.fit <- .tighten.feasibility(mx.model)
  }

  out <- list(mx.fit=mx.fit, RAM=RAM, data=data, cluster=cluster,
              mxalgebras=mxalgebras.ci,
              intervals.type=intervals.type, warning=warning.msg,
              error=error.msg)
  class(out) <- "mxsem"
  out
}

.sem.multigroup.twolevel <- function(model.name, RAM, data, group, cluster,
                                     intervals.type, startvalues, lbound,
                                     ubound, replace.constraints, mxModel.Args,
                                     run, silent, ...) {
  if (!is.data.frame(data) || length(group) != 1L ||
      !(group %in% names(data)) || length(cluster) != 1L ||
      !(cluster %in% names(data)))
    stop("Combined group and level models require 'data' and columns named by 'group' and 'cluster'.")
  if (anyNA(data[[group]]))
    stop("Group membership must be known for every row.")
  values <- unique(as.character(data[[group]]))
  ## Each group's clusters are formed from that group's rows only, so a
  ## reused ID denotes different clusters; a unit spanning groups cannot be
  ## represented. Make this explicit rather than silent.
  ids <- split(as.character(data[[cluster]]), as.character(data[[group]]))
  ids <- unlist(lapply(ids, unique), use.names=FALSE)
  shared <- unique(ids[duplicated(ids) & !is.na(ids)])
  if (length(shared))
    warning("Cluster ID(s) ", paste(utils::head(shared, 5), collapse=", "),
            if (length(shared) > 5) ", ...",
            " appear in more than one group and are treated as distinct ",
            "clusters within each group. Use distinct IDs across groups ",
            "(required for robust standard errors).", call.=FALSE)
  if (length(values) != length(RAM))
    stop("Number of groups in 'data' does not match 'RAM'.")
  .check.group.labels(RAM, values)
  global.algebras <- attr(RAM, "mxalgebras")
  remaining.algebras <- global.algebras
  replacements <- NULL
  if (isTRUE(replace.constraints) && length(global.algebras)) {
    free.labels <- unique(unlist(lapply(RAM, function(r) {
      mats <- list(r$within$A, r$within$S, r$within$M,
                   r$between$A, r$between$S, r$between$M)
      unlist(lapply(mats, function(x) {
        m <- as.mxMatrix(x, name="tmp")
        m$labels[m$free]
      }))
    })))
    selection <- .select.substitutions(global.algebras, free.labels, lbound, ubound)
    replacements <- selection$replacements
    remaining.algebras <- selection$remaining
  }
  startvalues2 <- startvalues
  if (isTRUE(replace.constraints)) {
    global.starts <- .ram.starts(unlist(lapply(RAM, function(r)
      list(r$within$A, r$within$S, r$within$M, r$between$A, r$between$S,
           r$between$M)), recursive=FALSE))
    startvalues2 <- c(global.starts, startvalues)
    startvalues2 <- startvalues2[!duplicated(names(startvalues2), fromLast=TRUE)]
  }
  ## Internal group and between-model names must not collide with user
  ## parameters, variables, parent entities, or auxiliary model names.
  reserved <- c(model.name, "data", "fitfunction", "expectation", "compute",
                names(global.algebras), names(remaining.algebras), unlist(lapply(global.algebras, function(x)
                  all.vars(x$formula)), use.names=FALSE))
  for (r in RAM) for (level in c("within", "between")) {
    for (nm in c("A", "S", "M")) {
      mat <- as.mxMatrix(r[[level]][[nm]], name="reserved")
      reserved <- c(reserved, mat$labels, colnames(r[[level]][[nm]]))
    }
  }
  extra.names <- function(x) {
    if (inherits(x, "MxModel"))
      return(c(x$name, unlist(lapply(x@matrices, extra.names), use.names=FALSE),
               unlist(lapply(x@algebras, extra.names), use.names=FALSE),
               unlist(lapply(x@submodels, extra.names), use.names=FALSE)))
    if (inherits(x, "MxMatrix")) return(c(x$name, x$labels))
    if (inherits(x, c("MxAlgebra", "MxConstraint"))) return(x$name)
    character()
  }
  reserved <- unique(na.omit(c(reserved, unlist(lapply(mxModel.Args, extra.names),
                                                 use.names=FALSE))))
  group.names <- character(length(values))
  for (i in seq_along(values)) {
    base <- paste0("group_", i)
    candidate <- base
    suffix <- 0L
    while (candidate %in% reserved || paste0("between_", candidate) %in% reserved) {
      suffix <- suffix + 1L
      candidate <- paste0(base, "_", suffix)
    }
    group.names[i] <- candidate
    reserved <- c(reserved, candidate, paste0("between_", candidate))
  }
  child <- lapply(seq_along(values), function(i) {
    rows <- as.character(data[[group]]) == values[i]
    ## The global algebras and constraints are added to the parent model
    ## below; the substitutions selected above are applied in every group.
    .sem.twolevel(model.name=group.names[i], RAM=RAM[[i]],
                  data=data[rows, , drop=FALSE], cluster=cluster,
                  between.name=paste0("between_", group.names[i]),
                  intervals.type=intervals.type,
                  startvalues=startvalues2, lbound=lbound, ubound=ubound,
                  replace.constraints=replace.constraints,
                  replacements=replacements,
                  mxModel.Args=NULL, run=FALSE, silent=silent)$mx.fit
  })
  model <- do.call(mxModel, c(list(model.name), child,
                              list(mxFitFunctionMultigroup(group.names))))
  algebra.ci <- NULL
  if (length(remaining.algebras)) {
    for (i in seq_along(remaining.algebras)) {
      model <- mxModel(model, remaining.algebras[[i]])
      name <- names(remaining.algebras)[i]
      if (!inherits(remaining.algebras[[i]], "MxConstraint")) {
        model <- mxModel(model, mxCI(name))
        algebra.ci <- c(algebra.ci, name)
      }
    }
  }
  if (!is.null(mxModel.Args))
    for (arg in mxModel.Args) model <- mxModel(model, arg)
  warning.msg <- error.msg <- NULL
  if (run) {
    fallback <- .run.mxTryHard.fallback(model, intervals.type, silent=silent, ...)
    mx.fit <- fallback$mx.fit
    warning.msg <- fallback$warning.msg
    error.msg <- fallback$error.msg
  } else {
    mx.fit <- .tighten.feasibility(model)
  }
  out <- list(mx.fit=mx.fit, RAM=RAM, data=data,
              group.map=setNames(group.names, values),
              cluster=cluster,
              mxalgebras=algebra.ci,
              intervals.type=intervals.type, warning=warning.msg,
              error=error.msg)
  class(out) <- "mxsem"
  out
}

## Match named summary inputs to the order of their covariance matrices.
## Unnamed inputs retain the historical positional convention.
.align.summary.entries <- function(Cov, value, argument, unit) {
  if (is.null(value)) return(NULL)
  reference <- names(Cov)
  supplied <- names(value)
  valid.names <- function(x, label) {
    if (is.null(x)) return(NULL)
    if (anyNA(x) || any(!nzchar(x)) || anyDuplicated(x))
      stop("Names of ", label, " must be complete and unique when supplied.", call.=FALSE)
    x
  }
  reference <- valid.names(reference, paste0("'Cov' ", unit, "s"))
  supplied <- valid.names(supplied, paste0("'", argument, "' ", unit, "s"))
  if (!is.null(reference) && !is.null(supplied)) {
    if (!setequal(reference, supplied))
      stop("Names of '", argument, "' must match the ", unit,
           " names in 'Cov'.", call.=FALSE)
    value <- value[match(reference, supplied)]
  } else if (is.null(reference) && !is.null(supplied)) {
    warning("'Cov' has no ", unit, " names; '", argument,
            "' is matched by position.", call.=FALSE)
  }
  value
}

## Reconstruct deterministic Gaussian pseudo-observations with the supplied
## per-cluster means and unbiased covariance matrices. For complete-data
## Gaussian ML these are sufficient statistics, so the resulting likelihood
## and estimates match those of the original observations. A pooled covariance
## cannot recover the within/between decomposition and is rejected.
.expand.cluster.summaries <- function(Cov, means, numObs, cluster) {
  if (!is.character(cluster) || length(cluster) != 1L || !nzchar(cluster))
    stop("'cluster' must name the generated cluster ID column.")
  if (!is.list(Cov) || !length(Cov) || !is.list(means) ||
      length(means) != length(Cov) || length(numObs) != length(Cov))
    stop("Two-level summaries require one covariance matrix, mean vector, and sample size per cluster.")
  means <- .align.summary.entries(Cov, means, "means", "cluster")
  numObs <- .align.summary.entries(Cov, numObs, "numObs", "cluster")
  ids <- names(Cov)
  if (is.null(ids) || any(!nzchar(ids))) ids <- as.character(seq_along(Cov))
  if (anyDuplicated(ids)) stop("Cluster IDs must be unique.")
  variables <- colnames(Cov[[1]])
  if (is.null(variables) || cluster %in% variables)
    stop("Per-cluster covariance matrices need variable names distinct from 'cluster'.")
  pieces <- lapply(seq_along(Cov), function(i) {
    s <- Cov[[i]]
    n <- numObs[i]
    mu <- means[[i]]
    p <- length(variables)
    if (!is.matrix(s) || !identical(dim(s), c(p, p)) ||
        !identical(colnames(s), variables) ||
        !identical(rownames(s), variables) ||
        !isSymmetric(s) || anyNA(s))
      stop("Each cluster covariance matrix must be complete, symmetric, and use the same variable order.")
    if (length(n) != 1L || !is.finite(n) || n < 1 || n != as.integer(n))
      stop("Each per-cluster sample size must be a positive integer.")
    if (!is.numeric(mu) || length(mu) != p || any(!is.finite(mu)))
      stop("Each cluster mean vector must match the covariance variable order.")
    mu <- .align.summary.entries(setNames(as.list(variables), variables), mu,
                                 "means", "variable")
    ## Work in correlation units to avoid discarding small variances merely
    ## because another variable is measured on a larger scale.
    if (any(!is.finite(s)) || any(diag(s) < 0))
      stop("Per-cluster covariance matrices must have finite nonnegative variances.")
    sd <- sqrt(diag(s))
    active <- which(sd > 0)
    inactive <- which(sd == 0)
    if (length(inactive) && any(s[inactive, , drop=FALSE] != 0))
      stop("A zero-variance variable must have zero covariances.")
    centered <- matrix(0, n, p)
    if (length(active)) {
      scale <- outer(sd[active], sd[active])
      corr <- s[active, active, drop=FALSE] / scale
      eig <- eigen(corr, symmetric=TRUE)
      tolerance <- 16 * .Machine$double.eps * length(active) * max(abs(eig$values))
      if (min(eig$values) < -tolerance)
        stop("Per-cluster covariance matrices must be positive semidefinite.")
      eig$values[eig$values < 0] <- 0
      ## Retain every positive direction when the sample has room. If the
      ## covariance is rank-deficient, allow numerical roundoff at zero.
      keep <- which(eig$values > 0)
      if (length(keep) > n-1L) {
        keep <- which(eig$values > tolerance)
        if (length(keep) > n-1L)
          stop("A cluster covariance matrix has rank greater than numObs - 1.")
      }
      if (length(keep)) {
        q <- matrix(0, n, length(keep))
        for (j in seq_along(keep)) {
          q[seq_len(j), j] <- 1/sqrt(j*(j+1))
          q[j+1L, j] <- -j/sqrt(j*(j+1))
        }
        z <- q %*% (diag(sqrt((n-1)*eig$values[keep]),
                          nrow=length(keep)) %*%
                     t(eig$vectors[, keep, drop=FALSE]))
        centered[, active] <- sweep(z, 2L, sd[active], "*")
      }
      if (n > 1L) {
        reproduced <- crossprod(centered[, active, drop=FALSE])/(n-1)
        error <- max(abs((reproduced - s[active, active, drop=FALSE]) / scale))
        if (!is.finite(error) || error > 1e-10)
          stop("The covariance matrix of cluster '", ids[i],
               "' cannot be reproduced by ", n, " observations.")
      }
    }
    out <- as.data.frame(sweep(centered, 2, unname(mu), "+"),
                         check.names=FALSE)
    names(out) <- variables
    ## Adding a large mean can itself erase small deviations in floating
    ## point arithmetic. Check the data actually passed to OpenMx too.
    if (n > 1L && length(active)) {
      reconstructed <- as.matrix(out[, active, drop=FALSE])
      reconstructed <- sweep(reconstructed, 2L, colMeans(reconstructed), "-")
      reconstructed <- crossprod(reconstructed)/(n-1)
      error <- max(abs((reconstructed - s[active, active, drop=FALSE]) / scale))
      if (!is.finite(error) || error > 1e-8)
        stop("The covariance matrix of cluster '", ids[i],
             "' cannot be reproduced accurately after adding its means.")
    }
    out[[cluster]] <- ids[i]
    out
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

## Will be depreciated in the future
create.mxModel <- function(model.name="sem", RAM=NULL, data=NULL,
                           Cov=NULL, means=NULL, numObs,
                           intervals.type=c("z", "LB"), startvalues=NULL,
                           replace.constraints=FALSE, mxModel.Args=NULL,
                           run=TRUE, silent=TRUE, ...) {
  .Deprecated("sem")
  sem(model.name=model.name, RAM=RAM, data=data, Cov=Cov, means=means,
      numObs=numObs, intervals.type=intervals.type, startvalues=startvalues,
      replace.constraints=replace.constraints, mxModel.Args=mxModel.Args,
      run=run, silent=silent, ...)
}


summary.mxsem <- function(object, robust=FALSE, smallSample=FALSE, t.df=NULL,
                          ...) {
  if (!is.element("mxsem", class(object)))
    stop("\"object\" must be an object of class \"mxsem\".")
  .check.mx.fit(object)

  ## Reference distribution of the Wald tests and intervals: z by default,
  ## or t with user-supplied degrees of freedom (e.g., K - p for
  ## meta-analytic fixed effects with few studies).
  if (!is.null(t.df)) {
    if (!is.numeric(t.df) || length(t.df) != 1L || is.na(t.df) ||
        !is.finite(t.df) || t.df <= 0)
      stop("'t.df' must be NULL or a single positive, finite number.")
    if (object$intervals.type != "z")
      stop("'t.df' requires intervals.type='z'; likelihood-based intervals do not use a reference distribution.")
  }
  crit <- if (is.null(t.df)) qnorm(.975) else qt(.975, t.df)

  if (isTRUE(robust) && object$intervals.type != "z")
    stop("Robust standard errors require intervals.type='z'.")
  if (isTRUE(smallSample) && !isTRUE(robust))
    stop("'smallSample' requires robust=TRUE.")
  robust.cov <- if (isTRUE(robust))
    vcov(object, robust=TRUE, smallSample=smallSample) else NULL
  se.method <- if (!isTRUE(robust)) "model-based" else if (
    .is.RAM.twolevel(object$RAM) ||
    inherits(object$RAM, "RAM_multigroup_twolevel") ||
    !is.null(object$cluster)) "cluster sandwich" else "observation sandwich"
  if (isTRUE(smallSample) && se.method == "cluster sandwich")
    se.method <- "cluster sandwich (small-sample correction)"

  # calculate coefficients
  ## OpenMx counts the rows of the between-level data in a two-level model
  ## as observations. Use the number of level-1 rows, as lavaan does, for the
  ## reported sample size and the information criteria.
  if (.is.RAM.twolevel(object$RAM) ||
      inherits(object$RAM, "RAM_multigroup_twolevel")) {
    info <- .twolevel.information(object)
    my.mx <- summary(object$mx.fit, numObs=sum(vapply(info, `[[`, integer(1), "nobs")))
  } else {
    my.mx <- summary(object$mx.fit)
  }
  ## Exclude lbound ubound etc
  my.para <- my.mx$parameters[, 1:6, drop=FALSE]   

  # Determine if CIs on parameter estimates are present
  if (object$intervals.type=="z") {

    ## Replace the SEs with robust SEs. With mxConstraints, OpenMx's SEs are
    ## replaced by those from .model.vcov().
    if (robust) {
      my.robust <- sqrt(diag(robust.cov))
      my.para[, "Std.Error"] <- my.robust[my.para$name]
    } else if (.has.constraints(object$mx.fit)) {
      model.cov <- .model.vcov(object$mx.fit)
      my.para[, "Std.Error"] <- sqrt(diag(model.cov))[my.para$name]
    }
        
    my.para$lbound <- my.para$Estimate - crit * my.para$Std.Error
    my.para$ubound <- my.para$Estimate + crit * my.para$Std.Error
    coefficients <- my.para[, -c(1:4), drop=FALSE]
    dimnames(coefficients)[[1]] <- my.para$name
      
  } else {
 
    ## Convert a data frame with length of 0 in my.mx$CI and remove the last column "note"
    my.ci <- my.mx$CI
    if (length(my.ci)==0) my.ci <- NULL else my.ci <- my.ci[, 1:3, drop=FALSE]        
    
    ## Select the elements matched my.para (excluded I2)  
    my.ci <- my.ci[row.names(my.ci) %in% my.para$name, ]
      
    my.ci <- data.frame(name=row.names(my.ci), my.ci)
    my.para <- merge(my.para, my.ci, by=c("name"), sort=FALSE)
    ## Keep the order of the parameters in the model, as with z intervals
    my.para <- my.para[order(match(my.para$name, my.mx$parameters$name)), ,
                       drop=FALSE]
    coefficients <- my.para[, -c(1:4,8)]
    dimnames(coefficients)[[1]] <- my.para$name
    # NA for LBCI
    coefficients$Std.Error <- NA
  }
  
  if (is.null(t.df)) {
    coefficients$"z value" <- coefficients$Estimate/coefficients$Std.Error
    coefficients$"Pr(>|z|)" <- 2*pnorm(-abs(coefficients$"z value"))
  } else {
    coefficients$"t value" <- coefficients$Estimate/coefficients$Std.Error
    coefficients$"Pr(>|t|)" <- 2*pt(-abs(coefficients$"t value"), t.df)
  }

  informationCriteria <- my.mx$informationCriteria
  ## Better column names
  colnames(informationCriteria) <- c("df Penalty", "Parameters Penalty",
                                     "Sample-Size Adjusted")

  ## Get the mxalgebras
  if (!is.null(object$mxalgebras)) {
    if (object$intervals.type=="z") {
      algebra.expr <- str2lang(paste0("rbind(", paste(object$mxalgebras, collapse=","), ")"))
      estimate <- eval(bquote(mxEval(.(algebra.expr), object$mx.fit)))
      ## mxSE()'s own default (omitting "cov") is the model's ordinary,
      ## non-robust vcov(), regardless of "robust" above -- so a defined
      ## parameter (":=", e.g. an indirect effect a*b) silently kept a
      ## non-robust delta-method SE/CI even when robust=TRUE correctly
      ## replaced every ordinary coefficient's SE with a robust one just
      ## above. Pass the same robust sandwich covariance (mxSE()'s "cov"
      ## argument accepts any free-parameter covariance matrix in place of
      ## vcov(model); imxRobustSE(..., details=TRUE)$cov is exactly that,
      ## with matching dimnames) so the two stay consistent.
      if (robust) {
        ## Parameters at a bound have NA robust covariance, which mxSE()
        ## would propagate to every defined parameter. They are treated as
        ## fixed (zero variance), and NA is reported only for defined
        ## parameters that depend on them.
        bounded <- rownames(robust.cov)[is.na(diag(robust.cov))]
        fixed.cov <- robust.cov
        fixed.cov[bounded, ] <- 0
        fixed.cov[, bounded] <- 0
        SE <- eval(bquote(mxSE(.(algebra.expr), model=object$mx.fit, cov=fixed.cov, silent=TRUE)))
        ## Shift each bounded parameter separately so that their changes
        ## cannot cancel (e.g., a - b with both a and b at a bound).
        if (length(bounded)) {
          theta <- coef(object$mx.fit)
          for (k in bounded) {
            shifted <- omxSetParameters(object$mx.fit, labels=k,
                                        values=theta[k] + 1e-4*max(1, abs(theta[k])))
            moved <- vapply(object$mxalgebras, mxEvalByName, numeric(1),
                             model=shifted, compute=TRUE)
            SE[abs(moved - estimate) > 1e-10 * pmax(1, abs(estimate))] <- NA
          }
        }
      } else if (.has.constraints(object$mx.fit)) {
        model.cov <- .model.vcov(object$mx.fit)
        SE <- eval(bquote(mxSE(.(algebra.expr), model=object$mx.fit, cov=model.cov, silent=TRUE)))
      } else {
        SE <- eval(bquote(mxSE(.(algebra.expr), model=object$mx.fit, silent=TRUE)))
      }
      mxalgebras <- cbind(lbound=estimate - crit*SE,
                          estimate=estimate,
                          ubound=estimate + crit*SE)
      dimnames(mxalgebras) <- list(object$mxalgebras,
                                   c("lbound", "estimate", "ubound"))
    } else {
      ## The CI of a 1 x 1 algebra is named "<model>.<algebra>[1,1]". Match
      ## the whole name: a partial match of "sem.ab" would also select
      ## "sem.ab2[1,1]".
      index <- match(paste0(object$mx.fit$name, ".", object$mxalgebras, "[1,1]"),
                     rownames(my.mx$CI))
      mxalgebras <- my.mx$CI[index, c("lbound", "estimate", "ubound"), drop=FALSE]
      dimnames(mxalgebras) <- list(object$mxalgebras,
                                   c("lbound", "estimate", "ubound"))
    }
  } else {
    mxalgebras <- NULL
  }

  out <- list(coefficients=coefficients, mxalgebras=mxalgebras,
              intervals.type=object$intervals.type,
              robust=robust, se.method=se.method, t.df=t.df,
              no.studies=my.mx$numObs, no.clusters=.n.clusters(object),
              obsStat=my.mx$observedStatistics,
              estPara=my.mx$estimatedParameters, df=my.mx$degreesOfFreedom,
              Minus2LL=my.mx$Minus2LogLikelihood,
              Mx.status1=object$mx.fit@output$status[[1]],
              informationCriteria=informationCriteria)
    class(out) <- "summary.mxsem"
    out
}

## Select actual group models, excluding auxiliary OpenMx submodels.
.twolevel.group.models <- function(object) {
  fit <- object$mx.fit
  if (!inherits(object$RAM, "RAM_multigroup_twolevel")) return(list(fit))
  groups <- if (!is.null(object$group.map)) unname(object$group.map) else
    sub("\\.fitfunction$", "", fit@fitfunction$groups)
  if (length(groups) != length(object$RAM) || !all(groups %in% names(fit@submodels)))
    stop("Cannot identify the fitted two-level group models.", call.=FALSE)
  fit@submodels[groups]
}

## Level-1 sample size counts rows with observed level-1 variables. A cluster
## also contributes when its between-level observed covariates are present,
## even if every level-1 observation in the cluster is missing.
.twolevel.information <- function(object, models=.twolevel.group.models(object)) {
  grouped <- inherits(object$RAM, "RAM_multigroup_twolevel")
  rams <- if (grouped) object$RAM else list(object$RAM)
  key <- object$cluster
  lapply(seq_along(models), function(i) {
    model <- models[[i]]
    ram <- rams[[i]]
    d <- model@data@observed
    within <- rownames(ram$within$F)
    if (is.null(d) || !key %in% names(d) || !all(within %in% names(d)))
      stop("Cannot identify the data used by the two-level model.", call.=FALSE)
    observed <- if (length(within)) rowSums(!is.na(d[, within, drop=FALSE])) > 0 else
      rep(FALSE, nrow(d))
    ids <- unique(as.character(d[[key]][observed]))
    between.vars <- rownames(ram$between$F)
    if (length(between.vars)) {
      between <- grep("^between($|_)", names(model@submodels), value=TRUE)
      if (length(between) != 1L)
        stop("Cannot identify the between-level model.", call.=FALSE)
      bd <- model@submodels[[between]]@data@observed
      present <- rowSums(!is.na(bd[, between.vars, drop=FALSE])) > 0
      ids <- union(ids, as.character(bd[[key]][present]))
    }
    list(nobs=as.integer(sum(observed)), clusters=ids)
  })
}

## Number of independent clusters of a fit with a cluster argument, as used
## by the cluster-robust covariance: the distinct IDs over all groups of a
## flat model (an ID shared by two groups is one unit), or the clusters of a
## two-level model summed over groups. NULL without a cluster argument.
.n.clusters <- function(object) {
  key <- object$cluster
  if (is.null(key)) return(NULL)
  fit <- object$mx.fit
  if (.is.RAM.twolevel(object$RAM) ||
      inherits(object$RAM, "RAM_multigroup_twolevel")) {
    info <- .twolevel.information(object)
    return(sum(vapply(info, function(x) length(x$clusters), integer(1))))
  }
  ids <- if (length(fit@submodels)) {
    groups <- fit@submodels[sub("\\.fitfunction$", "", fit@fitfunction$groups)]
    unlist(lapply(groups, function(m) as.character(m@data@observed[[key]])),
           use.names=FALSE)
  } else as.character(fit@data@observed[[key]])
  if (is.null(ids)) return(NULL)
  length(unique(ids))
}

print.summary.mxsem <- function(x, ...) {
    if (!is.element("summary.mxsem", class(x))) {
      stop("\"x\" must be an object of class \"summary.mxsem\".")
    }
    
    cat("95% confidence intervals: ")
    switch(x$intervals.type,
           z = if (is.null(x$t.df))
                 cat("z statistic approximation (robust=", x$robust, ")", sep="")
               else
                 cat("t distribution with ", format(x$t.df), " df (robust=",
                     x$robust, ")", sep=""),
           LB = cat("Likelihood-based statistic") )

    cat("\nStandard errors: ", x$se.method, "\nCoefficients:\n", sep="")
    printCoefmat(x$coefficients, P.values=TRUE, ...)

    if (!is.null(x$mxalgebras)) {
      cat("\nMxalgebras:\n")
      print(x$mxalgebras)
    }
    
    cat("\nInformation Criteria:\n")
    print(x$informationCriteria)
   
    cat("\nNumber of subjects (or studies):", x$no.studies)
    if (!is.null(x$no.clusters))
      cat("\nNumber of clusters:", x$no.clusters)
    cat("\nNumber of observed statistics:", x$obsStat)
    cat("\nNumber of estimated parameters:", x$estPara)
    cat("\nDegrees of freedom:", x$df)
    cat("\n-2 log likelihood:", x$Minus2LL, "\n")        
    cat("OpenMx status1:", x$Mx.status1, "(\"0\" or \"1\": The optimization is considered fine.\nOther values may indicate problems.)\n")

    if (!(x$Mx.status1 %in% c(0,1))) warning("OpenMx status1 is neither 0 or 1. You are advised to 'rerun' it again.\n")
}

coef.mxsem <- function(object, ...) {
  if (!is.element("mxsem", class(object)))
    stop("\"object\" must be an object of class \"mxsem\".")
  .check.mx.fit(object)
  coef(object$mx.fit)
}

vcov.mxsem <- function(object, robust=FALSE, smallSample=FALSE, ...) {
  if (!is.element("mxsem", class(object)))
    stop("\"object\" must be an object of class \"mxsem\".")
  .check.mx.fit(object)
  if (isTRUE(smallSample) &&
      (!isTRUE(robust) || .is.RAM.twolevel(object$RAM) ||
       inherits(object$RAM, "RAM_multigroup_twolevel") ||
       is.null(object$cluster)))
    stop("'smallSample' is available only for single-level cluster-robust covariance.")

  if (robust) {
    if (.is.RAM.twolevel(object$RAM) ||
        inherits(object$RAM, "RAM_multigroup_twolevel"))
      .robust.twolevel.vcov(object)
    else if (!is.null(object$cluster))
      clusterRobustSE(object$mx.fit, cluster=object$cluster, details=TRUE,
                      smallSample=smallSample)$cov
    else
      .robust.observation.vcov(object$mx.fit)
  } else {
    .model.vcov(object$mx.fit)
  } 
}

## Observation-level sandwich covariance for a flat raw-data ML model. As in
## OpenMx::imxRobustSE(), the meat is multiplied by N/(N-1).
.robust.observation.vcov <- function(model) {
  n <- .check.flat.ml(model)
  scores <- .row.m2ll.gradients(model) / -2
  .score.sandwich(model, scores, n/(n-1))$cov
}

## Cluster-level sandwich covariance for relational two-level OpenMx models.
## The Gaussian likelihood factors by independent clusters. OpenMx does not
## expose row-wise derivatives through a primaryKey/joinKey submodel, so each
## cluster score is the gradient of that cluster's exact -2 log likelihood at
## the full-data estimate. One mxRun() per cluster computes the whole gradient
## by finite differences in the backend (mxComputeNumericDeriv); the central
## difference is used, or a one-sided one if the other side cannot be
## evaluated (e.g., near a bound). Almost all of the cost is OpenMx's
## frontend processing in mxRun(), so each cluster's model is formed by
## replacing the data slots of a copy without intervals rather than by
## mxModel(). The sandwich is formed by .score.sandwich(), which also handles
## parameters at a bound.
.robust.twolevel.vcov <- function(object) {
  fit <- object$mx.fit
  theta <- coef(fit)
  if (!length(theta))
    stop("Robust covariance requires free parameters.")

  grouped <- inherits(object$RAM, "RAM_multigroup_twolevel")
  models <- .twolevel.group.models(object)
  key <- object$cluster
  if (is.null(key) || !all(vapply(models, function(m) key %in% names(m@data@observed),
                                  logical(1))))
    stop("Cluster identifier is unavailable for robust covariance.")

  ## Each group's clusters are formed from that group's rows only, so an ID
  ## used in two groups would be scored as two independent clusters.
  if (grouped) {
    ids <- lapply(.twolevel.information(object, models), `[[`, "clusters")
    shared <- unique(unlist(ids)[duplicated(unlist(ids))])
    if (length(shared))
      stop("Cluster-robust covariance requires clusters nested within groups, ",
           "but cluster ID(s) ", paste(utils::head(shared, 5), collapse=", "),
           if (length(shared) > 5) ", ...", " appear in more than one group. ",
           "Give each group's clusters distinct IDs.")
  }

  derivative <- mxComputeSequence(list(
    mxComputeOnce("fitfunction", "fit"),
    mxComputeNumericDeriv(hessian=FALSE, checkGradient=FALSE)))
  scores <- list()
  information <- .twolevel.information(object, models)
  for (i in seq_along(models)) {
    model <- models[[i]]
    template <- model
    template@intervals <- list()
    template@output <- list()
    template@compute <- derivative
    between.name <- grep("^between($|_)", names(model@submodels), value=TRUE)
    if (length(between.name) != 1L)
      stop("Cannot identify the between-level submodel for robust covariance.")
    between <- model@submodels[[between.name]]
    within.data <- model@data@observed
    between.data <- between@data@observed
    within.rows <- split(seq_len(nrow(within.data)),
                         as.character(within.data[[key]]))
    between.rows <- split(seq_len(nrow(between.data)),
                          as.character(between.data[[key]]))
    for (id in information[[i]]$clusters) {
      one <- template
      one@data <- mxData(within.data[within.rows[[id]], , drop=FALSE],
                         type="raw")
      one@submodels[[between.name]]@data <-
        mxData(between.data[between.rows[[id]], , drop=FALSE], type="raw",
               primaryKey=key)
      result <- tryCatch(mxRun(one, silent=TRUE, suppressWarnings=TRUE,
                               unsafe=TRUE),
                         error=function(e) e)
      if (inherits(result, "error"))
        stop("Could not evaluate the score of cluster '", id, "': ",
             conditionMessage(result))
      steps <- result@compute$steps
      d <- steps[[length(steps)]]$output$gradient
      g <- ifelse(is.finite(d$central), d$central,
                  ifelse(is.finite(d$forward), d$forward, d$backward))
      if (any(!is.finite(g)))
        stop("Could not evaluate the score of cluster '", id,
             "' near the fitted parameter values.")
      gradient <- setNames(numeric(length(theta)), names(theta))
      gradient[rownames(d)] <- g
      scores[[length(scores)+1L]] <- gradient
    }
  }
  if (length(scores) < 2L)
    stop("Robust covariance requires at least two independent clusters.")
  .score.sandwich(fit, do.call(rbind, scores) / -2)$cov
}

anova.mxsem <- function(object, ..., all=FALSE) {
  lapply(c(list(object), list(...)), .check.mx.fit)
  base <- lapply(list(object), function(x) x$mx.fit)
  comparison <- lapply(list(...), function(x) x$mx.fit)
  mxCompare(base=base, comparison=comparison, all=all)
}

plot.mxsem <- function(x, manNames=NULL, latNames=NULL,
                       labels=c("labels", "RAM"), what="est", nCharNodes=0,
                       nCharEdges=0, layout=c("tree", "circle", "spring",
                                              "tree2", "circle2"),
                       sizeMan=8, sizeLat=8, edge.label.cex=1.3,
                       color="white", weighted=FALSE,
                       group=NULL, level=NULL,
                       combine=missing(group) && missing(level),
                       main.line=2.5, main.adj=NULL, ...) {

  if (!requireNamespace("semPlot", quietly=TRUE))
    stop("\"semPlot\" package is required for this function.")

  if (!inherits(x, "mxsem"))
    stop("'mxsem' object is required.\n")
  .check.mx.fit(x)

  RAM <- x$RAM
  layout <- match.arg(layout)
  labels <- match.arg(labels)
  if (!is.null(level)) level <- match.arg(level, c("within", "between"))

  ## Raw level-1 data of a fitted (sub)model, used to evaluate the paths set
  ## by definition variables at their means; NULL for summary input.
  raw.data <- function(m)
    if (!is.null(m@data) && m@data$type == "raw") m@data@observed else NULL

  ## Arguments shared by every diagram
  draw <- function(mx.sub, RAM, shadow.vars=character(0), main=NULL,
                   within.aliases=NULL, data=raw.data(mx.sub))
    .plot.mxsem.one(mx.sub, RAM, data, shadow.vars,
                    within.aliases=within.aliases, what=what,
                    manNames=manNames, latNames=latNames, labels=labels,
                    nCharNodes=nCharNodes, nCharEdges=nCharEdges,
                    layout=layout, sizeMan=sizeMan, sizeLat=sizeLat,
                    edge.label.cex=edge.label.cex, color=color,
                    weighted=weighted, main=main, main.line=main.line,
                    main.adj=main.adj, ...)

  ## One level of a two-level (sub)model: "mx.sub" is the fitted within
  ## model, which contains its between model.
  draw.level <- function(mx.sub, RAM, level, main=NULL) {
    if (level == "within") {
      draw(mx.sub, RAM$within, main=main,
           within.aliases=.twolevel.within.aliases(RAM))
    } else {
      between.name <- grep("^between($|_)", names(mx.sub@submodels), value=TRUE)
      ## A between model may contain no manifest variables of its own, and
      ## semPlot::ramModel() errors for an all-latent model. Draw the
      ## between components of level-1 variables and the cluster-level
      ## covariates as manifest; level-2 factors stay latent.
      manifest <- intersect(colnames(RAM$between$A),
                            c(rownames(RAM$within$F), rownames(RAM$between$F)))
      ## Definition variables are averaged over the level-1 rows, as at
      ## level 1.
      draw(mx.sub@submodels[[between.name]], RAM$between, manifest, main=main,
           data=raw.data(mx.sub))
    }
  }

  ## Titles of the groups: the user's own group values (e.g. "site-2"),
  ## not the sanitized submodel names (e.g. "site_2").
  group.titles <- function(sub.names) {
    if (is.null(x$group.map)) return(sub.names)
    orig <- names(x$group.map)[match(sub.names, x$group.map)]
    ifelse(is.na(orig), sub.names, orig)
  }

  ## combine=TRUE plots every submodel (all groups, or within and between)
  ## as panels of one figure via par(mfrow=), the same manual approach
  ## semPlot's own examples use for multi-group lavaan fits (semPaths()'s
  ## panelGroups=TRUE only recognises objects it builds internally, not
  ## semPlotModel objects built here via ramModel()).
  if (inherits(RAM, "RAM_multigroup_twolevel")) {
    sub.names <- unname(x$group.map)
    if (isTRUE(combine)) {
      titles <- group.titles(sub.names)
      op <- graphics::par(mfrow=c(length(sub.names), 2))
      on.exit(graphics::par(op))
      for (i in seq_along(sub.names))
        for (lev in c("within", "between"))
          draw.level(x$mx.fit[[sub.names[i]]], RAM[[i]], lev,
                     main=paste0(titles[i], ": ", lev))
      return(invisible(NULL))
    }
    idx <- .resolve.plot.group(group, sub.names, x$group.map)
    return(draw.level(x$mx.fit[[sub.names[idx]]], RAM[[idx]],
                      if (is.null(level)) "within" else level))
  }

  if (.is.RAM.multigroup(RAM)) {
    sub.names <- if (!is.null(x$group.map)) unname(x$group.map) else
      names(x$mx.fit@submodels)[seq_along(RAM)]
    if (isTRUE(combine)) {
      titles <- group.titles(sub.names)
      op <- graphics::par(mfrow=c(1, length(sub.names)))
      on.exit(graphics::par(op))
      for (i in seq_along(sub.names))
        draw(x$mx.fit[[sub.names[i]]], RAM[[i]], main=titles[i])
      return(invisible(NULL))
    }
    ## Plot one group, matched by its original data value first (see
    ## .resolve.plot.group()).
    idx <- .resolve.plot.group(group, sub.names, x$group.map)
    return(draw(x$mx.fit[[sub.names[idx]]], RAM[[idx]]))
  }

  if (.is.RAM.twolevel(RAM)) {
    if (isTRUE(combine)) {
      op <- graphics::par(mfrow=c(1, 2))
      on.exit(graphics::par(op))
      draw.level(x$mx.fit, RAM, "within", main="within")
      draw.level(x$mx.fit, RAM, "between", main="between")
      return(invisible(NULL))
    }
    return(draw.level(x$mx.fit, RAM, if (is.null(level)) "within" else level))
  }

  draw(x$mx.fit, RAM)
}

## Names of the internal within components of the level-1 variables that
## also have a between component (see .decompose.twolevel.within()), or
## NULL if there are none.
.twolevel.within.aliases <- function(RAM) {
  shadow.vars <- intersect(rownames(RAM$within$F), colnames(RAM$between$A))
  attr(.decompose.twolevel.within(RAM$within, shadow.vars), "shadow.aliases")
}

## Resolve plot.mxsem()'s "group" argument to an index into sub.names
## (names(x$mx.fit@submodels)). Factored out of plot.mxsem() so the
## resolution logic -- easy to get subtly wrong, see the comments below --
## can be unit-tested directly instead of only indirectly through
## side-effecting plot() calls.
##
## "group.map" (fit$group.map) maps original data values to their
## (possibly sanitized) submodel name; consulted FIRST, not as a fallback
## after checking sub.names, because a sanitized name can collide with a
## different group's own literal value: original groups "a.b" and "a_b"
## sanitize/de-duplicate to submodels "a_b" and "a_b_1" -- matching
## sub.names first would let group="a_b" silently resolve to "a.b"'s
## submodel instead of the group literally named "a_b", since "a_b" is
## *also* a valid submodel name in its own right.
##
## The one exception is the no-argument default (group=NULL): sub.names[1]
## is already the right submodel name directly and must NOT be run through
## group.map, because in the exact collision scenario above, sub.names[1]
## ("a_b") can itself coincidentally be a valid *original* group value too
## (mapping to "a_b_1"), which would wrongly redirect the default to a
## different group than "the first one".
.resolve.plot.group <- function(group, sub.names, group.map) {
  if (is.null(group)) {
    group <- sub.names[1]
  } else {
    group <- as.character(group)
    if (!is.null(group.map) && group %in% names(group.map)) {
      group <- group.map[[group]]
    }
  }
  idx <- match(group, sub.names)
  if (is.na(idx)) {
    choices <- sub.names
    if (!is.null(group.map) && !setequal(names(group.map), sub.names)) {
      choices <- paste0(sub.names, " (data value: ",
                        names(group.map)[match(sub.names, group.map)], ")")
    }
    stop("'group' must be one of: ", paste(choices, collapse=", "), ".\n")
  }
  idx
}

## Build and draw one submodel's path diagram (one group, or one level).
## Factored out of plot.mxsem() so both the single-diagram path and the
## combine=TRUE multi-panel path share exactly one implementation.
.plot.mxsem.one <- function(mx.sub, RAM, data, shadow.vars, what, nCharNodes,
                            nCharEdges, layout, sizeMan, sizeLat,
                            edge.label.cex, color, weighted, main=NULL,
                            main.line=2.5, main.adj=NULL,
                            within.aliases=NULL, manNames=NULL, latNames=NULL,
                            labels=c("labels", "RAM"), ...) {

  ## Single-group/multi-group submodels name their matrices "Amatrix" etc.
  ## (as.mxMatrix(RAM$A, name="Amatrix")); two-level submodels are built via
  ## type="RAM" + mxPath(), which names them "A" etc. instead. Checked
  ## against EITHER @matrices or @algebras: replace.constraints=TRUE can
  ## turn any of A/S/M into an mxAlgebra (not a plain mxMatrix), so
  ## checking @matrices alone would misidentify the naming convention
  ## entirely whenever "Amatrix" itself happens to be the one replaced.
  has.name <- function(name) name %in% names(mx.sub@matrices) || name %in% names(mx.sub@algebras)
  if (has.name("Amatrix")) {
    mat.names <- c(A="Amatrix", S="Smatrix", F="Fmatrix", M="Mmatrix")
  } else {
    mat.names <- c(A="A", S="S", F="F", M="M")
  }

  ## With replace.constraints=TRUE, A, S, or M can be an mxAlgebra, which
  ## may contain definition variables. OpenMx evaluates them at the first
  ## data row; evaluate them at the means of the definition variables in
  ## "data" instead, as for the paths set by definition variables below.
  at.means <- mx.sub
  defvars <- unique(unlist(lapply(mx.sub@matrices, function(m) c(m$labels))))
  defvars <- sub("^data\\.", "", grep("^data\\.", defvars, value=TRUE))
  defvars <- intersect(defvars, names(data))
  if (length(defvars) && length(mx.sub@algebras) && !is.null(mx.sub@data) &&
      mx.sub@data$type == "raw") {
    observed <- mx.sub@data@observed[1, , drop=FALSE]
    for (v in defvars) observed[[v]] <- mean(data[[v]], na.rm=TRUE)
    at.means@data <- mxData(observed, type="raw",
                            primaryKey=mx.sub@data@primaryKey)
  }

  ## A plain mxMatrix's $values already carries dimnames; an mxAlgebra's
  ## result never does (confirmed by testing) -- fall back to "RAM"'s own
  ## dimnames for that specific matrix, which are unaffected by
  ## replace.constraints and already match "mx.sub" here (RAM is the
  ## caller's per-group/per-level slice, same as mx.sub).
  get.mat <- function(name, dimnames.fallback) {
    if (name %in% names(mx.sub@matrices)) return(mx.sub@matrices[[name]]$values)
    if (!name %in% names(mx.sub@algebras)) return(NULL)
    values <- if (identical(at.means, mx.sub)) mx.sub@algebras[[name]]$result
              else mxEvalByName(name, at.means, compute=TRUE)
    if (is.null(values) || !length(values))
      values <- mxEvalByName(name, at.means, compute=TRUE)
    if (!is.null(values) && is.null(dimnames(values))) {
      dims <- dimnames(mx.sub@algebras[[name]])
      dimnames(values) <- if (is.null(dims)) dimnames.fallback else dims
    }
    values
  }
  A <- get.mat(mat.names[["A"]], dimnames(RAM$A))
  S <- get.mat(mat.names[["S"]], dimnames(RAM$S))
  F <- get.mat(mat.names[["F"]], dimnames(RAM$F))
  M <- get.mat(mat.names[["M"]], dimnames(RAM$M))
  mean.ref <- mx.sub@expectation$M
  has.mean <- !is.null(M) && length(mean.ref) == 1L && !is.na(mean.ref)
  ## Fixed when M is NULL, i.e., no mean structure
  if (is.null(M)) M <- matrix(0, nrow=1, ncol=ncol(A))
  if (!has.mean) M[,] <- 0

  ## The fitted two-level within model has extra latent components for
  ## shared manifest variables. Show the original level-1 paths by reading
  ## their fitted cells from those components and hiding the fixed unit
  ## paths that only implement the within/between decomposition.
  if (length(within.aliases) > 0L) {
    original <- colnames(RAM$A)
    mapped <- original
    hit <- original %in% names(within.aliases)
    mapped[hit] <- within.aliases[original[hit]]
    A <- A[mapped, mapped, drop=FALSE]
    S <- S[mapped, mapped, drop=FALSE]
    M <- M[, mapped, drop=FALSE]
    dimnames(A) <- dimnames(RAM$A)
    dimnames(S) <- dimnames(RAM$S)
    dimnames(M) <- dimnames(RAM$M)
    F <- RAM$F
  }

  ## A path set by a definition variable is shown at the variable's mean in
  ## "data", the level-1 data of the group plotted. This is display only and
  ## does not affect the fitted model.
  for (i in seq_len(nrow(S)))
    for (j in seq_len(ncol(S))) {
      ## fixed=TRUE: matches the literal "data." substring, consistent
      ## with the strsplit() below it -- without this, "." is a regex
      ## wildcard and a literal-looking label such as "database" (data +
      ## b + ase) would also match, then fail the strsplit(fixed=TRUE)
      ## right after it (no literal "data." to split on, so tmp is NA).
      if (grepl("data.", RAM$S[i, j], fixed=TRUE)) {
        tmp <- strsplit(RAM$S[i, j], "data.", fixed=TRUE)[[1]][2]
        S[i, j] <- mean(data[[tmp]], na.rm=TRUE)
      }
    }

  for (i in seq_len(nrow(A)))
    for (j in seq_len(ncol(A))) {
      if (grepl("data.", RAM$A[i, j], fixed=TRUE)) {
        tmp <- strsplit(RAM$A[i, j], "data.", fixed=TRUE)[[1]][2]
        A[i, j] <- mean(data[[tmp]], na.rm=TRUE)
      }
    }

  for (j in seq_len(ncol(M))) {
    if (grepl("data.", RAM$M[1, j], fixed=TRUE)) {
      tmp <- strsplit(RAM$M[1, j], "data.", fixed=TRUE)[[1]][2]
      M[1, j] <- mean(data[[tmp]], na.rm=TRUE)
    }
  }

  ## F rows identify the order of observed variables; its columns identify
  ## their positions among all RAM variables. These orders can differ.
  allNames <- colnames(A)
  if (length(shadow.vars) > 0) {
    F <- create.Fmatrix(allNames %in% shadow.vars, as.mxMatrix=FALSE)
    dimnames(F) <- list(allNames[allNames %in% shadow.vars], allNames)
  }
  selector <- all(is.finite(F)) && all(F %in% c(0, 1)) &&
              all(rowSums(F) == 1) && all(colSums(F) <= 1)
  if (!selector)
    stop("Plotting requires F to select distinct observed variables.", call.=FALSE)
  obs.index <- if (nrow(F)) max.col(F, ties.method="first") else integer()
  index_obs <- seq_along(allNames) %in% obs.index
  if (is.null(manNames)) manNames <- allNames[obs.index]
  if (is.null(latNames)) latNames <- allNames[!index_obs]
  if (length(manNames) != length(obs.index) || length(latNames) != sum(!index_obs))
    stop("Lengths of 'manNames' and 'latNames' must match the observed and latent variables.")
  node.names <- allNames
  node.names[obs.index] <- manNames
  node.names[!index_obs] <- latNames

  ## ramModel() drops numeric zero entries and treats every retained path as
  ## free. Supply tiny placeholders for zero-valued specified paths, then
  ## restore their estimates and set the free/fixed and equality metadata.
  metadata <- lapply(c("A", "S", "M"), function(nm) {
    original <- as.mxMatrix(RAM[[nm]], name=paste0("plot", nm))
    free <- original$free
    labs <- original$labels
    specified <- original$free | !is.na(original$labels)
    derived <- matrix(FALSE, nrow=nrow(free), ncol=ncol(free))
    actual <- mx.sub@matrices[[mat.names[[nm]]]]
    if (!is.null(actual)) {
      target <- if (length(within.aliases)) mapped else allNames
      if (nm == "M") {
        if (!is.null(colnames(actual$free)) && all(target %in% colnames(actual$free))) {
          free <- actual$free[, target, drop=FALSE]
          labs <- actual$labels[, target, drop=FALSE]
        } else if (identical(dim(actual$free), dim(free))) {
          free <- actual$free
          labs <- actual$labels
        }
      } else if (!is.null(rownames(actual$free)) && all(target %in% rownames(actual$free))) {
        free <- actual$free[target, target, drop=FALSE]
        labs <- actual$labels[target, target, drop=FALSE]
      } else if (identical(dim(actual$free), dim(free))) {
        free <- actual$free
        labs <- actual$labels
      }
    } else if (mat.names[[nm]] %in% names(mx.sub@algebras)) {
      ## as.mxAlgebra() represents each cell as a separate named algebra.
      ## A direct label replacement must share the same plotting parameter
      ## ID as that label elsewhere; literal constants are fixed.
      for (row in seq_len(nrow(free))) for (col in seq_len(ncol(free))) {
        if (length(within.aliases)) {
          expanded <- c(colnames(RAM$A), unname(within.aliases))
          cell.row <- if (nm == "M") 1L else match(mapped[row], expanded)
          cell.col <- match(mapped[col], expanded)
        } else {
          cell.row <- row
          cell.col <- col
        }
        cell <- mx.sub@algebras[[paste0(mat.names[[nm]], cell.row, "_", cell.col)]]
        if (is.null(cell)) next
        expr <- cell$formula
        if (is.symbol(expr)) {
          label <- as.character(expr)
          free[row, col] <- !startsWith(label, "data.")
          labs[row, col] <- label
        } else if (is.numeric(expr) && length(expr) == 1L) {
          free[row, col] <- FALSE
          labs[row, col] <- NA_character_
        } else if (is.call(expr)) {
          variables <- all.vars(expr)
          free[row, col] <- any(!startsWith(variables, "data."))
          derived[row, col] <- free[row, col]
        }
      }
    }
    list(free=free, labels=labs, specified=specified, derived=derived)
  })
  names(metadata) <- c("A", "S", "M")
  shown <- list(A=A, S=S, M=M)
  for (nm in names(shown)) {
    mask <- metadata[[nm]]$specified | metadata[[nm]]$free | !is.na(metadata[[nm]]$labels)
    if (nm == "M" && !has.mean) mask[,] <- FALSE
    zeros <- mask & shown[[nm]] == 0
    shown[[nm]][zeros] <- 1e-100
  }
  sem.plot <- semPlot::ramModel(A=shown$A, S=shown$S, F=F, M=shown$M,
                                manNames=manNames, latNames=latNames,
                                Names=node.names)
  equal.ids <- new.env(parent=emptyenv())
  next.id <- 0L
  for (i in seq_len(nrow(sem.plot@Pars))) {
    path <- sem.plot@Pars[i, ]
    from <- match(path$lhs, node.names)
    to <- match(path$rhs, node.names)
    nm <- switch(path$edge, "->"="A", "<->"="S", "int"="M", NULL)
    if (is.null(nm)) next
    row <- if (nm == "M") 1L else to
    col <- if (nm == "M") to else from
    actual.value <- switch(nm, A=A[row, col], S=S[row, col], M=M[row, col])
    free <- metadata[[nm]]$free[row, col]
    label <- metadata[[nm]]$labels[row, col]
    sem.plot@Pars$est[i] <- actual.value
    if (actual.value == 0 && nm != "M") sem.plot@Pars$std[i] <- 0
    sem.plot@Pars$fixed[i] <- !free
    if (free && !metadata[[nm]]$derived[row, col]) {
      key <- if (is.na(label) || !nzchar(label)) paste0(".unlabelled_", i) else label
      if (!exists(key, envir=equal.ids, inherits=FALSE)) {
        next.id <- next.id + 1L
        assign(key, next.id, envir=equal.ids)
      }
      sem.plot@Pars$par[i] <- get(key, envir=equal.ids, inherits=FALSE)
    } else sem.plot@Pars$par[i] <- 0
    if (match.arg(labels) == "labels")
      sem.plot@Pars$label[i] <- if (is.na(label)) "" else label
  }

  ## qgraph drops an edge whose weight is exactly zero, even when it exists
  ## in the semPlotModel. A negligible weight retains the edge in estimate
  ## and standardized diagrams and is displayed as zero at normal precision.
  render.plot <- sem.plot
  if (grepl("est|par", what, ignore.case=TRUE))
    render.plot@Pars$est[!is.na(render.plot@Pars$est) & render.plot@Pars$est == 0] <- 1e-100
  if (grepl("stand|std", what, ignore.case=TRUE))
    render.plot@Pars$std[!is.na(render.plot@Pars$std) & render.plot@Pars$std == 0] <- 1e-100
  out <- semPlot::semPaths(render.plot, what=what, nCharNodes=nCharNodes,
                           nCharEdges=nCharEdges, layout=layout,
                           sizeMan=sizeMan, sizeLat=sizeLat,
                           edge.label.cex=edge.label.cex, color=color,
                           weighted=weighted, ...)
  ## Also preserve exact values for callers that intercept semPaths() and
  ## return its semPlotModel (e.g. for programmatic inspection).
  if (inherits(out, "semPlotModel")) out <- sem.plot
  ## Default line=2.5: semPaths()/qgraph draw close to the top margin (e.g.
  ## a variance edge label there), so the default title() position collides
  ## with it; push the title up and out of the way. main.adj (NA = par("adj"),
  ## i.e. centered, same as title()'s own default) additionally lets the
  ## caller nudge the title horizontally (e.g. main.adj=1 for right-aligned)
  ## when it still overlaps an edge label centered near the top of the
  ## diagram -- which edge that is depends on the model, so no single fixed
  ## adj avoids every collision.
  if (!is.null(main)) graphics::title(main=main, line=main.line, adj=main.adj)
  invisible(out)
}

## Groups are matched to the data by position (order of first appearance),
## as in lavaan, not by the labels of "group:" blocks. Warn when those labels
## are the data's group values in a different order, which is likely a
## mistake.
.check.group.labels <- function(RAM, values) {
  labels <- attr(RAM, "group.labels")
  values <- as.character(values)
  if (is.null(labels) || length(labels) != length(values) ||
      !setequal(labels, values) || identical(labels, values))
    return(invisible(NULL))
  differ <- labels != values
  warning("Groups are matched to the data by position, as in lavaan, not by ",
          "their labels: ", paste0("model group '", labels[differ],
                                   "' is fitted to data group '", values[differ],
                                   "'", collapse=", "),
          ". Reorder the 'group:' blocks or the data if this is not intended.",
          call.=FALSE)
  invisible(NULL)
}

## Run "mx.model" via mxRun(). If that errors, retry with mxTryHard() and a
## final mxRun(); if it returns a status code other than 0 or 1, retry with
## mxTryHard() from its estimates (see below). Warnings during the retries
## are recorded rather than allowed to abort them, and are printed only if
## silent=FALSE. A model that cannot be fitted gives a NULL fit with a
## warning; .check.mx.fit() reports it when the fit is used. Shared by all
## sem() fitting paths.
## A model without free parameters has status code NA.
.status.ok <- function(fit) inherits(fit, "MxModel") &&
  (fit@output$status$code %in% c(0, 1) || !length(coef(fit)))

## A fit returned with a status code other than 0 or 1 is retried with
## mxTryHard() from its estimates. The retry is kept if it has status 0 or 1
## or a lower -2LL; otherwise the original fit is kept. 'record' wraps each
## call to collect its warnings and errors, and 'rerun', if given, reruns a
## kept retry (e.g., to compute likelihood-based CIs).
.retry.bad.status <- function(mx.fit, record, rerun=NULL) {
  tried <- record(.quiet(mxTryHard(mx.fit, extraTries=10, intervals=FALSE,
                                   silent=TRUE)))
  if (inherits(tried, "MxModel") &&
      (.status.ok(tried) || isTRUE(tried@output$fit < mx.fit@output$fit - 1e-8))) {
    if (!is.null(rerun)) tried <- record(rerun(tried))
    if (inherits(tried, "MxModel")) mx.fit <- tried
  }
  mx.fit
}

## Warn if a fit that was retried (status = its status code before the
## retry; NULL if not retried) still does not have status 0 or 1.
.warn.bad.status <- function(mx.fit, status) {
  if (!is.null(status) && !.status.ok(mx.fit))
    warning("OpenMx status code is ", mx.fit@output$status$code,
            if (!identical(mx.fit@output$status$code, status))
              paste0(" (", status, " before mxTryHard())"),
            ", so the model may not have converged. You are advised to ",
            "rerun it, e.g., with rerun(), and check that the status code ",
            "is 0 or 1.", call.=FALSE)
  invisible(NULL)
}

.run.mxTryHard.fallback <- function(mx.model, intervals.type, silent=TRUE, ...) {
  warning.msg <- error.msg <- NULL
  intervals <- intervals.type == "LB"
  mx.model <- .tighten.feasibility(mx.model)
  ## Record warnings and errors; print them only if silent=FALSE.
  record <- function(expr) withCallingHandlers(
    tryCatch(expr, error=function(e) {
      error.msg <<- conditionMessage(e)
      if (!silent) message(error.msg)
      e
    }),
    warning=function(w) {
      warning.msg <<- c(warning.msg, conditionMessage(w))
      if (!silent) message(conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  status <- NULL

  mx.fit <- tryCatch(mxRun(mx.model, intervals=intervals,
                           suppressWarnings=TRUE, silent=silent, ...),
                     error=function(e) e)

  if (inherits(mx.fit, "error")) {
    mx.fit <- record(.quiet(mxTryHard(mx.model, extraTries=50, intervals=FALSE,
                                      silent=TRUE)))
    if (!inherits(mx.fit, "error"))
      mx.fit <- record(mxRun(mx.fit, intervals=intervals,
                             suppressWarnings=TRUE, silent=silent, ...))
    if (inherits(mx.fit, "error")) {
      warning("Error in running mxModel.\n")
      mx.fit <- NULL
    }
  } else if (!.status.ok(mx.fit)) {
    status <- mx.fit@output$status$code
    mx.fit <- .retry.bad.status(mx.fit, record=record,
                                rerun=if (intervals) function(fit)
                                  mxRun(fit, intervals=TRUE, suppressWarnings=TRUE,
                                        silent=silent, ...))
  }
  mx.fit <- .check.constraints(mx.fit, intervals=intervals)
  .warn.bad.status(mx.fit, status)

  list(mx.fit=mx.fit, warning.msg=warning.msg, error.msg=error.msg)
}

## Stop with a clear, actionable message if "object" (an "mxsem" object)
## failed to fit -- i.e. its mx.fit is not a valid MxModel, the scenario
## .run.mxTryHard.fallback() above now correctly signals (see its comment)
## instead of leaving mx.fit as an untyped NULL that downstream OpenMx/base
## R calls in summary.mxsem()/coef.mxsem()/vcov.mxsem()/plot.mxsem() would
## otherwise choke on with unrelated-looking internal errors.
.check.mx.fit <- function(object) {
  if (!inherits(object$mx.fit, "MxModel")) {
    reason <- if (!is.null(object$error)) object$error
             else if (!is.null(object$warning))
               paste(object$warning, collapse="; ")
             else "unknown reason"
    stop("The model was not fitted successfully, so there is no fitted ",
        "model to use here. Reported reason: ", reason, call.=FALSE)
  }
}

## Combine bounds recorded on a RAM object with those given as arguments
## (named lists keyed by label); the arguments take precedence.
.check.parameter.bounds <- function(RAM, extras, lbound, ubound) {
  for (b in list(lbound, ubound)) {
    if (!length(b)) next
    if (is.null(names(b)) || anyNA(names(b)) || any(!nzchar(names(b))) || anyDuplicated(names(b)))
      stop("Bounds must have unique, nonempty parameter names.", call.=FALSE)
    if (any(!vapply(b, function(x) is.numeric(x) && length(x)==1L &&
                    (is.na(x) || is.finite(x)), logical(1))))
      stop("Each bound must be a finite number or NA.", call.=FALSE)
  }
  collect <- function(x) {
    if (is.null(x)) return(character())
    if (inherits(x, "MxModel"))
      return(c(unlist(lapply(x@matrices, collect), use.names=FALSE),
               unlist(lapply(x@submodels, collect), use.names=FALSE),
               unlist(lapply(x@algebras, collect), use.names=FALSE),
               unlist(lapply(x@constraints, collect), use.names=FALSE)))
    if (inherits(x, c("MxAlgebra", "MxConstraint"))) return(all.vars(x$formula))
    if (inherits(x, "MxMatrix")) return(as.character(x$labels[!is.na(x$labels)]))
    if (is.matrix(x)) {
      if (!is.character(x)) return(character())
      entries <- x[!is.na(x) & is.na(suppressWarnings(as.numeric(x)))]
      return(sub("^.*[*@]", "", entries))
    }
    if (is.list(x)) return(unlist(lapply(x, collect), use.names=FALSE))
    if (methods::is(x, "MxPath")) return(as.character(x@labels))
    character()
  }
  known <- unique(c(collect(RAM), collect(extras)))
  missing <- setdiff(c(names(lbound), names(ubound)), known)
  if (length(missing))
    stop("Bound parameter name(s) not found in the model: ", paste(unique(missing), collapse=", "), ".", call.=FALSE)
  invisible(NULL)
}

.merge.bounds <- function(from.RAM, given) {
  if (is.null(from.RAM)) return(given)
  given <- as.list(given)
  c(from.RAM[setdiff(names(from.RAM), names(given))], given)
}

## Add lbound and ubound to Amatrix, Smatrix, and Mmatrix
.addbound <- function(x, lbound=NULL, ubound=NULL) {
  if (!is.null(lbound)) {
    for (i in seq_along(lbound)) {
      index <- x$labels==names(lbound)[i]
      ## NAs are treated as FALSE
      index[is.na(index)] <- FALSE
      x$lbound[index] <- lbound[[i]]
    }
  }
  
  if (!is.null(ubound)) {
    for (i in seq_along(ubound)) {
      index <- x$labels==names(ubound)[i]
      index[is.na(index)] <- FALSE
      x$ubound[index] <- ubound[[i]]
    }
  }
  x
}

## Detect a multiple-group RAM object, i.e., one produced by
## lavaan2RAM(..., ngroups > 1): a list of per-group RAM lists, each with
## (at least) A/S/F/M elements. Class-tagged objects (the normal case, set
## by lavaan2RAM()) are recognised immediately; a structural fallback covers
## RAM lists assembled by hand (e.g. via lapply(list_of_models, lavaan2RAM))
## that were not tagged with class "RAM_multigroup" explicitly.
.is.RAM.multigroup <- function(RAM) {
  if (inherits(RAM, "RAM_multigroup")) return(TRUE)
  ## A two-level RAM's "within"/"between" elements each structurally look
  ## like a valid per-group RAM list too (both are lists with A/S/F/M), so
  ## it must be excluded explicitly here, not just left to rely on
  ## dispatch order in sem().
  if (.is.RAM.twolevel(RAM)) return(FALSE)
  is.list(RAM) && length(RAM) > 0 &&
    all(vapply(RAM, function(x) is.list(x) && all(c("A","S","F","M") %in% names(x)),
               FUN.VALUE=logical(1)))
}

## Build one group's self-contained mxModel (its own mxData, expectation,
## and fit function) from a single-group RAM list. Mirrors the matrix-
## building portion of sem() for the single-group case, including its
## replace.constraints=TRUE branching (unchanged / fully-constant /
## algebra) -- "replacements" is the GLOBAL, already chain-resolved
## label->rhs map computed once by .sem.multigroup() (NULL when
## replace.constraints=FALSE, or when there is nothing to replace); it is
## simply not found in this group's own A/S/M text when irrelevant to it.
##
## Unlike .sem.twolevel()'s equivalent (.twolevel.replace.one()), no
## remove-and-readd/dimname patch is needed here: this function builds
## Amatrix/Smatrix/Mmatrix as ordinary named mxMatrix objects up front and
## passes dimnames=var.names directly to mxExpectationRAM() (the same
## construction style the single-group path already uses, where its own
## replace.constraints branch needs no such patch either).
.build.group.mxmodel <- function(model.name, RAM, data, startvalues=NULL,
                                 lbound=NULL, ubound=NULL, replacements=NULL) {
  Amatrix <- as.mxMatrix(RAM$A, name="Amatrix")
  Smatrix <- as.mxMatrix(RAM$S, name="Smatrix")
  Fmatrix <- as.mxMatrix(RAM$F, name="Fmatrix")
  Mmatrix <- as.mxMatrix(RAM$M, name="Mmatrix")

  no.means <- inherits(data, "metaSEM_group_summary") && is.null(data$means)
  if (no.means) {
    Mmatrix$free[,] <- FALSE
    Mmatrix$labels[,] <- NA_character_
  }

  checkRAM(Amatrix, Smatrix, cor.analysis=FALSE)

  var.names <- colnames(Fmatrix$values)

  if (inherits(data, "metaSEM_group_summary")) {
    mx.data <- mxData(observed=data$Cov, type="cov", numObs=data$numObs,
                      means=if (is.null(data$means)) NA else data$means)
    if (is.null(data$means)) {
      expFun <- mxExpectationRAM(A="Amatrix", S="Smatrix", F="Fmatrix",
                                 dimnames=var.names)
    } else {
      expFun <- mxExpectationRAM(A="Amatrix", S="Smatrix", F="Fmatrix",
                                 M="Mmatrix", dimnames=var.names)
    }
  } else {
    mx.data <- mxData(observed=data, type="raw")
    expFun <- mxExpectationRAM(A="Amatrix", S="Smatrix", F="Fmatrix",
                               M="Mmatrix", dimnames=var.names)
  }
  mx.model <- mxModel(model.name, mx.data, expFun, mxFitFunctionML())

  ## Starting values: RAM's own defaults, then overridden by any matching
  ## labels in the user-supplied "startvalues" (shared across all groups --
  ## a label shared across groups via lavaan's c(l,l) syntax picks up the
  ## same override in every group it appears in, which is the correct
  ## behaviour for an equality-constrained parameter). When
  ## replace.constraints applied, "startvalues" here also carries the
  ## existing-parameter fallback values .sem.multigroup() collected (see
  ## its own comment on the starting-value collision this avoids).
  para.labels <- c(Amatrix$labels[Amatrix$free], Smatrix$labels[Smatrix$free],
                   Mmatrix$labels[Mmatrix$free])
  if (!is.null(startvalues)) {
    to.apply <- startvalues[names(startvalues) %in% para.labels]
    for (i in seq_along(to.apply)) {
      Amatrix$values[Amatrix$labels==names(to.apply)[i]] <- to.apply[[i]]
      Smatrix$values[Smatrix$labels==names(to.apply)[i]] <- to.apply[[i]]
      Mmatrix$values[Mmatrix$labels==names(to.apply)[i]] <- to.apply[[i]]
    }
  }

  Amatrix <- .addbound(Amatrix, lbound=lbound, ubound=ubound)
  Smatrix <- .addbound(Smatrix, lbound=lbound, ubound=ubound)
  Mmatrix <- .addbound(Mmatrix, lbound=lbound, ubound=ubound)

  ## Substitute "replacements" into local symbolic copies of A/S/M (bare
  ## labels/constants, no "value*" prefix -- same as.symMatrix() single-
  ## group's own replace.constraints branch uses). A no-op (all(orig==new)
  ## stays TRUE for every matrix) when replacements is NULL or irrelevant
  ## to this group, so the branch below reduces to exactly today's
  ## behaviour in that case.
  A <- as.symMatrix(RAM$A); S <- as.symMatrix(RAM$S); M <- as.symMatrix(RAM$M)
  if (no.means) M[] <- "0"
  A0 <- A; S0 <- S; M0 <- M
  A <- .apply.replacements(A, replacements)
  S <- .apply.replacements(S, replacements)
  M <- .apply.replacements(M, replacements)

  ## Same unchanged/fully-constant/algebra branch single-group sem() uses.
  add.one <- function(mx.model, orig, new, mat, name) {
    if (all(orig==new)) {
      mxModel(mx.model, mat)
    } else if (length(all.vars(parse(text=new))) == 0) {
      ## Fully constant after substitution (e.g. a constraint fixing a
      ## parameter at a literal number, or at a numeric expression like
      ## "exp(0.4)"). as.mxMatrix() cannot evaluate a non-trivial constant
      ## EXPRESSION -- it only recognises bare numbers or "value*label"
      ## strings, silently treating anything else (e.g. "exp(0.4)", which
      ## contains no "*") as an UNLABELLED FREE parameter with no starting
      ## value instead of a fixed one. Confirmed by testing: without
      ## evaluating by hand first, "sigma1 == exp(0.4)" left sigma1 free
      ## and estimated at the data's own MLE, not fixed at exp(0.4) --
      ## same fix already applied to .twolevel.replace.one().
      values <- matrix(vapply(new, function(z) eval(parse(text=z)), numeric(1)),
                       nrow=nrow(new), ncol=ncol(new), dimnames=dimnames(new))
      mxModel(mx.model, as.mxMatrix(values, name=name))
    } else {
      alg <- as.mxAlgebra(new, startvalues=startvalues, lbound=lbound,
                          ubound=ubound, name=name)
      mxModel(mx.model, alg$mxalgebra, alg$parameters, alg$list)
    }
  }
  mx.model <- add.one(mx.model, A0, A, Amatrix, "Amatrix")
  mx.model <- add.one(mx.model, S0, S, Smatrix, "Smatrix")
  mx.model <- add.one(mx.model, M0, M, Mmatrix, "Mmatrix")

  ## Same expCov/expMean bookkeeping as the single-group path -- not used by
  ## the fit function, but kept for parity/future use (e.g. implied stats).
  ## References "Amatrix"/"Smatrix"/"Mmatrix" by NAME, so it resolves
  ## correctly regardless of whether add.one() above added a plain matrix
  ## or an algebra under that name.
  Id <- as.mxMatrix(diag(ncol(Fmatrix$values)), name="Id")
  Id_A <- mxAlgebra(solve(Id - Amatrix), name="Id_A")
  expCov <- mxAlgebra(Id_A %&% Smatrix, name="expCov")
  expMean <- mxAlgebra(Mmatrix %*% t(Id_A), name="expMean")

  mx.model <- mxModel(mx.model, Fmatrix, Id, Id_A, expCov, expMean)

  ## Post-substitution free-parameter labels (identical to "para.labels"
  ## above when nothing was substituted for this group -- as.symMatrix()
  ## strips to bare labels either way, and all.vars() naturally excludes
  ## numeric constants -- so this single formula covers both cases,
  ## mirroring single-group sem()'s own "new.para.labels").
  new.para.labels <- unique(c(A, S, M))
  new.para.labels <- all.vars(parse(text=new.para.labels))
  new.para.labels <- new.para.labels[!startsWith(new.para.labels, "data.")]

  list(model=mx.model, para.labels=new.para.labels)
}

## Fit a multiple-group RAM object: one mxModel(type="RAM")-style submodel
## per group (built by .build.group.mxmodel()), combined with OpenMx's
## mxFitFunctionMultigroup(). Cross-group equality constraints (lavaan's
## c(l,l) syntax) fall out automatically because equal labels across
## submodels are the same free parameter to OpenMx.
.sem.multigroup <- function(model.name="sem", RAM, data=NULL, group=NULL,
                            cluster=NULL, Cov=NULL, means=NULL, numObs=NULL,
                            intervals.type=c("z","LB"), startvalues=NULL,
                            lbound=NULL, ubound=NULL, replace.constraints=FALSE,
                            mxModel.Args=NULL, run=TRUE, silent=TRUE, ...) {

  intervals.type <- match.arg(intervals.type)

  if (is.null(data)) {
    if (!is.null(cluster))
      stop("Cluster-robust covariance requires raw individual data.")
    if (!is.list(Cov) || length(Cov) != length(RAM) ||
        is.null(numObs) || length(numObs) != length(RAM) ||
        any(!is.finite(numObs) | numObs <= 0))
      stop("Summary input for multiple-group models requires one covariance matrix and positive sample size per group.")
    if (!is.null(means) && (!is.list(means) || length(means) != length(RAM)))
      stop("'means' must be a list with one mean vector per group.")
    means <- .align.summary.entries(Cov, means, "means", "group")
    numObs <- .align.summary.entries(Cov, numObs, "numObs", "group")
    group.values.chr <- names(Cov)
    if (is.null(group.values.chr) || any(!nzchar(group.values.chr)))
      group.values.chr <- as.character(seq_along(Cov))
    data.list <- lapply(seq_along(Cov), function(i) {
      if (!is.matrix(Cov[[i]]) || nrow(Cov[[i]]) != ncol(Cov[[i]]))
        stop("Each element of 'Cov' must be a square covariance matrix.")
      group.means <- if (is.null(means)) NULL else means[[i]]
      if (!is.null(group.means) && !is.null(colnames(Cov[[i]])))
        group.means <- .align.summary.entries(
          setNames(as.list(colnames(Cov[[i]])), colnames(Cov[[i]])),
          group.means, "means", "variable")
      structure(list(Cov=Cov[[i]], means=group.means,
                     numObs=numObs[i]), class="metaSEM_group_summary")
    })
    names(data.list) <- group.values.chr
  } else {
    if (!is.data.frame(data))
      stop("Multiple-group raw data must be a data frame.")
    if (is.null(group) || !(group %in% colnames(data)))
      stop("'group' must name a column in 'data' identifying group membership.")
    if (!is.null(cluster) &&
        (!is.character(cluster) || length(cluster) != 1L ||
         !(cluster %in% names(data))))
      stop("'cluster' must name a column in raw data.")
  ## Group values in FIRST-APPEARANCE order in the data -- confirmed
  ## empirically against lavaan (not just assumed): lavaan numbers groups
  ## by first appearance in the data for BOTH character and factor group
  ## columns, ignoring a factor's own levels() order entirely. This
  ## function originally used split(data, data[[group]]), which sorts
  ## alphabetically for a character column (or follows levels() for a
  ## factor) -- silently mismatching a lavaan-authored model's group order
  ## whenever the first-appearing group in the data isn't also the
  ## alphabetically-first one. Example: group values appear as "b" then
  ## "a" in the data, with model labels c(int_b, int_a) (group 1 = "b",
  ## group 2 = "a", matching where each label was written); split() would
  ## instead put "a" first, silently assigning int_b's estimate to group
  ## "a" and int_a's to group "b". Caught on review (with lavaan itself as
  ## the reference, not assumed) -- every test up to that point happened
  ## to use group values already in alphabetical/first-appearance
  ## agreement, so this never surfaced.
  g <- data[[group]]
  if (anyNA(g)) {
    warning(sum(is.na(g)), " row(s) with a missing '", group,
           "' value were dropped (not assigned to any group).")
  }
  group.values <- unique(g)
  group.values <- group.values[!is.na(group.values)]
  group.values.chr <- as.character(group.values)

  data.list <- lapply(group.values.chr, function(v) {
    data[!is.na(g) & as.character(g)==v, , drop=FALSE]
  })
  names(data.list) <- group.values.chr
  }

  if (!is.null(cluster) &&
      any(vapply(data.list, function(x) anyNA(x[[cluster]]), logical(1))))
    stop("Cluster IDs must be known for every retained row.")

  if (length(data.list) != length(RAM)) {
    stop("Number of groups implied by the first-appearance-order unique ",
        "values of 'data[[\"", group, "\"]]' (", length(data.list), ") ",
        "does not match the number of groups in 'RAM' (", length(RAM), "). ",
        "These values must correspond, in order, to the groups used when ",
        "the model was specified (lavaan numbers groups by first ",
        "appearance in the data too).")
  }
  ## Originally added because an unused factor level in the group column
  ## produced a 0-row entry that could pass the count check above
  ## undetected (whenever it happened to leave the right number of
  ## groups) and silently fit a degenerate submodel with no data --
  ## mxRun() did not even error on that, it returned a catastrophic-
  ## failure status with every parameter frozen at its starting value.
  ## That was against the old split(data, data[[group]]) implementation,
  ## which enumerates every level of a factor column whether or not it
  ## actually occurs in the data. group.values above is now built from
  ## unique(data[[group]]) instead (first-appearance order, see the
  ## comment above it), which by construction only ever contains values
  ## that occur at least once -- an unused factor level is excluded before
  ## data.list is even built, and would instead surface as a group-count
  ## mismatch on the check above. This 0-row check should therefore be
  ## unreachable now; kept as a defensive fallback (e.g. against a future
  ## refactor of how group.values is computed) rather than removed.
  empty <- if (is.null(data)) rep(FALSE, length(data.list)) else
    vapply(data.list, nrow, integer(1)) == 0
  if (any(empty)) {
    stop("Group(s) ", paste(names(data.list)[empty], collapse=", "),
        " have 0 rows in 'data'.")
  }

  ## Every group's own genuine free-parameter labels (needed below to keep
  ## a sanitized group name from colliding with one).
  harvest.labels <- function(x) {
    m <- as.mxMatrix(x)
    m$labels[m$free]
  }
  para.labels <- unique(unlist(lapply(RAM, function(r) {
    c(harvest.labels(r$A), harvest.labels(r$S), harvest.labels(r$M))
  })))

  ## ---- replace.constraints=TRUE: substitute "==" constraints whose LHS
  ## is a genuine free parameter in ANY group directly into that group's
  ## (and every other qualifying group's) A/S/M matrices, mirroring sem()'s
  ## single-group substitution logic and .sem.twolevel()'s generalization
  ## of it -- applied here across groups instead of levels. A constraint
  ## may reference a label that lives in a DIFFERENT group entirely (e.g. a
  ## ratio constraint across groups' variances); cross-group label
  ## resolution needs no special handling, since OpenMx labels are already
  ## global across a model tree, exactly as for two-level.
  ##
  ## mxAlgebra/mxConstraint from lavaan2RAM's group==0 rows (fn1 := ...,
  ## a == b, etc.) are global, not group-specific, and lavaan2RAM() only
  ## ever attaches them to the first group's RAM element.
  ##
  ## Deliberately computed BEFORE group names are sanitized below (not
  ## just before .build.group.mxmodel() is called): a constraint's RHS
  ## can introduce a genuinely NEW free-parameter label that never
  ## appears anywhere in the original RAM matrices at all (e.g. "a0" in
  ## "sigma1 == exp(a0)") -- "para.labels" above cannot know about these,
  ## so group-name reservation must wait until they are known too, or a
  ## group value equal to one of them (e.g. "a0") hits the exact same
  ## named-entity/free-parameter collision as an ordinary label (see
  ## "reserved" below) -- confirmed by testing.
  mxalgebras <- RAM[[1]]$mxalgebras
  replacements <- NULL
  startvalues2 <- startvalues
  introduced.labels <- character(0)

  if (isTRUE(replace.constraints) && !is.null(mxalgebras)) {
    selection <- .select.substitutions(mxalgebras, para.labels, lbound, ubound)
    replacements <- selection$replacements
    mxalgebras <- selection$remaining
    introduced.labels <- selection$introduced

    if (length(replacements)) {
      ## Starting values of all groups, for labels used across groups (see
      ## .ram.starts()); those in "startvalues" take precedence.
      existing.starts <- .ram.starts(unlist(lapply(RAM, function(r)
        list(r$A, r$S, r$M)), recursive=FALSE))
      startvalues2 <- c(existing.starts, startvalues)
      startvalues2 <- startvalues2[!duplicated(names(startvalues2), fromLast=TRUE)]
    }
  }

  ## Sanitize group values into legal, unique OpenMx model names. mxModel()
  ## rejects far more than just numeric-looking names -- e.g. "." and "-"
  ## are illegal too ("site.1"/"site-2" both fail), while some other
  ## punctuation is fine (checked directly rather than assumed), so
  ## anything not in [A-Za-z0-9_] is replaced rather than trying to
  ## enumerate a safe/unsafe character list. A mapping from the original
  ## group value back to its (possibly rewritten) submodel name is kept in
  ## out$group.map below so plot.mxsem(group=) can still be called with the
  ## user's own group labels, not just the sanitized name.
  sanitize.name <- function(x) {
    x <- gsub("[^A-Za-z0-9_]", "_", x)
    if (x=="" || !is.na(suppressWarnings(as.numeric(x)))) x <- paste0("group_", x)
    x
  }
  ## A sanitized-but-legal name can still collide with a name OpenMx (or
  ## .build.group.mxmodel()) treats specially: "data", "fitfunction",
  ## "expectation" and "compute" are GLOBALLY reserved by mxModel() itself
  ## (confirmed directly -- mxModel("data") errors with "is illegal
  ## because it is a reserved name", regardless of context); "sem" (this
  ## call's own model.name) collides because a model cannot contain a
  ## child entity sharing its own name; "Amatrix"/"Smatrix"/etc. collide
  ## the same way once .build.group.mxmodel() tries to name a matrix
  ## "Amatrix" INSIDE a submodel that is ITSELF already named "Amatrix";
  ## and a group's own FREE-PARAMETER LABELS collide because OpenMx
  ## rejects a name used as both a named entity (the submodel) and a free
  ## parameter within it (confirmed by testing: "In model 'b' the
  ## following are used both as named entities and free parameters: 'b'",
  ## reproducible with or without replace.constraints -- e.g. a group
  ## value of "b" together with a shared slope labelled "b" via lavaan's
  ## c(b,b) syntax). "introduced.labels" (computed above) covers the same
  ## risk for a label that only comes into existence via a
  ## replace.constraints substitution (e.g. a group literally named "a0"
  ## colliding with "sigma1 == exp(a0)") -- para.labels alone cannot see
  ## these, since they never appear in the original RAM matrices at all.
  ## Reserve all of them up front, before make.unique(), so a group value
  ## equal to one of them gets a "_1" suffix like any other collision
  ## instead of surfacing as an opaque OpenMx error.
  reserved <- c(model.name, "data", "fitfunction", "expectation", "compute",
               "Amatrix", "Smatrix", "Fmatrix", "Mmatrix",
               "Id", "Id_A", "expCov", "expMean", para.labels, introduced.labels,
               names(mxalgebras))
  .check.group.labels(RAM, group.values.chr)
  group.names <- vapply(group.values.chr, sanitize.name, character(1))
  group.names <- make.unique(c(reserved, unname(group.names)), sep="_")[-seq_along(reserved)]
  group.map <- setNames(group.names, group.values.chr)

  built <- lapply(seq_along(RAM), function(i) {
    .build.group.mxmodel(model.name=group.names[i], RAM=RAM[[i]],
                         data=data.list[[i]], startvalues=startvalues2,
                         lbound=lbound, ubound=ubound, replacements=replacements)
  })

  submodels <- lapply(built, function(x) x$model)
  mx.model <- do.call(mxModel, c(list(model.name), submodels,
                                 list(mxFitFunctionMultigroup(group.names))))

  ## mxCI over the union of free-parameter labels across all groups (shared
  ## labels are OpenMx parameters known by name regardless of which
  ## submodel(s) declare them, so a single mxCI() at the container level is
  ## sufficient -- no need to qualify by submodel name). Each built$para.
  ## labels already reflects post-substitution labels when
  ## replace.constraints applied (see .build.group.mxmodel()).
  all.para.labels <- unique(unlist(lapply(built, function(x) x$para.labels)))
  ## mxCI() rejects an empty reference vector outright rather than treating
  ## it as "no CIs requested" (confirmed by testing) -- reachable whenever
  ## every group's free parameters have all been fixed via
  ## replace.constraints=TRUE, a valid, if unusual, model to fit.
  if (length(all.para.labels) > 0) {
    mx.model <- mxModel(mx.model, mxCI(all.para.labels))
  }

  ## Constraints already consumed by the replace.constraints substitution
  ## above have already been removed from "mxalgebras" at this point.
  mxalgebras.ci <- NULL
  if (!is.null(mxalgebras)) {
    for (i in seq_along(mxalgebras)) {
      mx.model <- mxModel(mx.model, mxalgebras[[i]])
      name.mxalgebra <- names(mxalgebras)[i]
      if (!inherits(mxalgebras[[i]], "MxConstraint")) {
        mx.model <- mxModel(mx.model, mxCI(c(name.mxalgebra)))
        mxalgebras.ci <- c(mxalgebras.ci, name.mxalgebra)
      }
    }
  }

  if (!is.null(mxModel.Args)) {
    for (i in seq_along(mxModel.Args)) {
      mx.model <- mxModel(mx.model, mxModel.Args[[i]])
    }
  }

  warning.msg <- error.msg <- NULL

  if (run) {
    fallback <- .run.mxTryHard.fallback(mx.model, intervals.type, silent=silent, ...)
    mx.fit <- fallback$mx.fit
    warning.msg <- fallback$warning.msg
    error.msg <- fallback$error.msg
  } else {
    mx.fit <- .tighten.feasibility(mx.model)
  }

  out <- list(mx.fit=mx.fit, RAM=RAM, data=data, cluster=cluster,
              mxalgebras=mxalgebras.ci,
              intervals.type=intervals.type, warning=warning.msg,
              error=error.msg, group.map=group.map)
  class(out) <- "mxsem"
  out
}

## Detect a two-level RAM object, i.e., one produced by lavaan2RAM() from
## "level: 1"/"level: 2" syntax: a list with (at least) "within" and
## "between" elements, each themselves a list with A/S/F/M. Class-tagged
## objects (the normal case, set by lavaan2RAM()) are recognised
## immediately; a structural fallback covers hand-assembled two-level RAM
## objects that were not tagged explicitly.
.is.RAM.twolevel <- function(RAM) {
  if (inherits(RAM, "RAM_twolevel")) return(TRUE)
  is.list(RAM) && all(c("within","between") %in% names(RAM)) &&
    is.list(RAM$within) && all(c("A","S","F","M") %in% names(RAM$within)) &&
    is.list(RAM$between) && all(c("A","S","F","M") %in% names(RAM$between))
}

## A level-1 manifest that also has a between component must be represented
## internally as the sum of two variables. In particular, a regression from
## x to y at level 1 must use x's WITHIN component: using the observed x as
## the source also carries x's between component through that path and fits
## a different model. Keep the public RAM unchanged; this expansion is only
## for OpenMx's relational model built by .sem.twolevel().
.decompose.twolevel.within <- function(RAM, shadow.vars) {
  if (length(shadow.vars) == 0L) return(RAM)

  old.names <- colnames(RAM$A)
  shadow.vars <- intersect(old.names, shadow.vars)
  aliases <- make.unique(c(old.names, paste0(shadow.vars, "_within")), sep="_")
  aliases <- utils::tail(aliases, length(shadow.vars))
  names(aliases) <- shadow.vars
  mapped <- old.names
  mapped[old.names %in% shadow.vars] <- aliases[old.names[old.names %in% shadow.vars]]
  new.names <- c(old.names, unname(aliases))

  A <- S <- matrix("0", length(new.names), length(new.names),
                    dimnames=list(new.names, new.names))
  M <- matrix("0", 1L, length(new.names),
              dimnames=list(rownames(RAM$M), new.names))
  A[mapped, mapped] <- RAM$A
  S[mapped, mapped] <- RAM$S
  M[1L, mapped] <- RAM$M[1L, ]
  for (i in seq_along(shadow.vars)) {
    A[shadow.vars[i], aliases[i]] <- "1"
  }

  F <- matrix(0, nrow(RAM$F), length(new.names),
              dimnames=list(rownames(RAM$F), new.names))
  F[, old.names] <- RAM$F
  out <- list(A=A, S=S, F=F, M=M)
  attr(out, "shadow.aliases") <- aliases
  out
}

## Translate a RAM-style character A/S/M matrix trio into a list of
## mxPath() objects (the "local" structure of one level -- no cross-level
## joinKey paths, which are added separately by the caller).
##
## This exists because OpenMx's type="RAM" model construction only reliably
## wires up an mxExpectationRAM() (in particular its mean vector, and any
## cross-model joinKey linkage) when the ENTIRE local structure is expressed
## via mxPath() -- see dev/technical.md, "OpenMx construction".
## Directly supplying named "A"/"S"/"M" mxMatrix objects instead (the style
## used elsewhere in this package, e.g. in sem()'s single-group path) was
## tried and silently breaks the join: the between-level parameters simply
## never enter the optimizer's gradient. So, for the two-level path only,
## every free/fixed cell is expressed as its own mxPath() call instead.
.RAM2paths <- function(A, S, M, startvalues=NULL, lbound=NULL, ubound=NULL) {
  var.names <- colnames(A)
  nr <- length(var.names)

  ## as.mxMatrix() does not preserve dimnames, so variable names are taken
  ## from the original character matrices ("A" etc.), not from Amx/Smx/Mmx.
  Amx <- .addbound(as.mxMatrix(A, name="tmpA"), lbound=lbound, ubound=ubound)
  Smx <- .addbound(as.mxMatrix(S, name="tmpS"), lbound=lbound, ubound=ubound)
  Mmx <- .addbound(as.mxMatrix(M, name="tmpM"), lbound=lbound, ubound=ubound)

  if (!is.null(startvalues)) {
    for (i in seq_along(startvalues)) {
      nm <- names(startvalues)[i]
      Amx$values[Amx$labels==nm] <- startvalues[[i]]
      Smx$values[Smx$labels==nm] <- startvalues[[i]]
      Mmx$values[Mmx$labels==nm] <- startvalues[[i]]
    }
  }

  paths <- list()

  ## A: asymmetric (regressions/loadings), full scan. Amatrix[row=DV,
  ## col=IV] -- a path FROM the column variable TO the row variable.
  ## Structurally-zero, unlabelled cells are skipped: type="RAM" defaults
  ## any untouched pair to 0/fixed on its own, and skipping keeps the path
  ## list a manageable size. The skip is keyed on "no label" rather than
  ## "value==0 and not free", because a fixed cell CAN have a nonzero value
  ## (0) and still need to survive: definition variables (e.g. "0*data.mod")
  ## and "@"-tied fixed constants (see as.mxMatrix()) are both free=FALSE
  ## with a real label attached, and a 0 numeric prefix is a common,
  ## legitimate convention for them (the starting-value slot is irrelevant
  ## when the cell is data-driven or tied to another label). An earlier
  ## version of this check used "!free && value==0" and silently dropped
  ## exactly this case.
  for (i in seq_len(nr)) {
    for (j in seq_len(nr)) {
      if (!Amx$free[i,j] && Amx$values[i,j]==0 && is.na(Amx$labels[i,j])) next
      paths[[length(paths)+1]] <- mxPath(from=var.names[j], to=var.names[i], arrows=1,
                                         free=Amx$free[i,j], values=Amx$values[i,j],
                                         labels=Amx$labels[i,j], lbound=Amx$lbound[i,j],
                                         ubound=Amx$ubound[i,j])
    }
  }

  ## S: symmetric (variances/covariances); lower triangle incl. diagonal
  ## only, to avoid emitting each covariance path twice. Same "no label"
  ## skip rationale as A above.
  for (i in seq_len(nr)) {
    for (j in seq_len(i)) {
      if (!Smx$free[i,j] && Smx$values[i,j]==0 && is.na(Smx$labels[i,j])) next
      paths[[length(paths)+1]] <- mxPath(from=var.names[j], to=var.names[i], arrows=2,
                                         free=Smx$free[i,j], values=Smx$values[i,j],
                                         labels=Smx$labels[i,j], lbound=Smx$lbound[i,j],
                                         ubound=Smx$ubound[i,j])
    }
  }

  ## M: means. Unlike A/S, EVERY variable gets an explicit "one" path, even
  ## fixed-at-zero ones -- type="RAM" only builds a mean vector into the
  ## expectation at all when at least one mxPath(from="one", ...) is
  ## present; a directly-supplied M matrix with no such path is otherwise
  ## silently ignored (raising "no expected means vector was provided" at
  ## mxRun() time).
  for (j in seq_len(nr)) {
    paths[[length(paths)+1]] <- mxPath(from="one", to=var.names[j],
                                       free=Mmx$free[1,j], values=Mmx$values[1,j],
                                       labels=Mmx$labels[1,j], lbound=Mmx$lbound[1,j],
                                       ubound=Mmx$ubound[1,j])
  }

  list(paths=paths,
       para.labels=c(Amx$labels[Amx$free], Smx$labels[Smx$free], Mmx$labels[Mmx$free]))
}

## Substitute every bare-symbol occurrence of a name in "subs" appearing
## within language object "e" with the corresponding (already-resolved)
## replacement language object in "subs". A plain recursive AST walk, not
## string regex -- so a label can't be mismatched against a substring of
## an unrelated identifier, and operator precedence in the surrounding
## expression is automatically preserved (e.g. "2*residW" with residW's
## replacement "a+b" correctly becomes "2 * (a + b)", never "2*a+b",
## because the substitution happens on the parsed call tree, where
## "residW" is a single argument node, not by textually splicing "a+b" in
## place of a matched substring).
.substitute.symbols <- function(e, subs) {
  if (is.symbol(e)) {
    nm <- as.character(e)
    if (nm %in% names(subs)) return(subs[[nm]])
    return(e)
  } else if (is.call(e)) {
    args <- as.list(e)
    e[] <- c(args[1], lapply(args[-1], .substitute.symbols, subs=subs))
    return(e)
  } else {
    return(e)
  }
}

## Resolve CHAINED replace.constraints=TRUE substitutions to a fixed point
## before they are applied to the RAM matrices. Given the labels and RHS
## text of every "label == RHS" constraint about to be substituted, expand
## any RHS that itself references another constraint's own LHS label, so
## every RETURNED rhs is expressed purely in terms of parameters that are
## NOT themselves being eliminated in this same batch.
##
## Without this, applying each constraint independently -- each LHS label
## replaced by its OWN, unexpanded RHS text, in one pass over the RAM
## matrices, which is what both sem() and .sem.twolevel() used to do --
## can silently reintroduce an eliminated label as a brand-new,
## disconnected free parameter wherever ANOTHER constraint's RHS still
## mentions it by name. Confirmed by testing: "residW == exp(a0); residB
## == 2*residW" left "residW" as an independent free parameter (unrelated
## to "exp(a0)") in whichever matrix the second constraint touched,
## silently fitting a different model than the one written.
##
## "labels" and "rhs" are parallel character vectors, one per qualifying
## constraint; order-independent (a chain may be written in either
## dependency direction). Raises a clear error, naming the offending
## labels, for a constraint set that can never resolve to a fixed point
## (e.g. "a == b; b == a", or a parameter referencing its own label, even
## indirectly through another constraint) rather than silently leaving an
## unresolved self-reference for the matrix-substitution step to mishandle
## the same way. The per-symbol loop below is bounded (at most
## length(labels)+1 passes) regardless of whether a cycle is present, so
## it cannot hang; the definitive cycle/self-reference check happens once,
## after the loop, by testing whether a label still appears within its own
## (as far as possible) resolved expression.
.resolve.chained.constraints <- function(labels, rhs) {
  exprs <- lapply(rhs, function(r) parse(text=r)[[1]])
  names(exprs) <- labels

  for (iter in seq_len(length(labels) + 1L)) {
    changed <- FALSE
    for (i in seq_along(exprs)) {
      hit <- intersect(all.vars(exprs[[i]]), labels[-i])
      if (length(hit) > 0) {
        exprs[[i]] <- .substitute.symbols(exprs[[i]], exprs[hit])
        changed <- TRUE
      }
    }
    if (!changed) break
  }

  self.ref <- vapply(seq_along(exprs), function(i) labels[i] %in% all.vars(exprs[[i]]),
                     logical(1))
  if (any(self.ref)) {
    stop("Cyclic (or self-referential) 'replace.constraints=TRUE' ",
        "substitution detected among: ", paste(labels[self.ref], collapse=", "),
        ". Each parameter's constraint must not reference its own label, ",
        "even indirectly through another constraint.", call.=FALSE)
  }
  ## deparse() (unlike as.character() on a whole call, used to extract
  ## "rhs" in the first place) wraps source text onto multiple lines once
  ## it exceeds its default width.cutoff -- returning a character vector
  ## of length > 1 for a sufficiently long resolved expression, which
  ## vapply(..., character(1)) then rejects outright ("result is length
  ## N"), failing a perfectly valid long constraint before the model is
  ## even built. deparse1() (R >= 4.0.0; this package requires >= 4.2.0)
  ## exists exactly for this -- deparse, then collapse onto one line.
  vapply(exprs, deparse1, character(1))
}

## Replace each label in a symbolic RAM matrix (see as.symMatrix()) by its
## right-hand side in "replacements" (label -> RHS text).
.apply.replacements <- function(x, replacements) {
  for (lab in names(replacements)) x[x == lab] <- replacements[[lab]]
  x
}

## Starting values of the labeled free parameters ("value*label") in RAM
## matrices, as a list named by label. With replace.constraints=TRUE, they
## are passed to as.mxAlgebra(), which would otherwise start every label it
## finds at 0 and give a parameter two different starting values.
.ram.starts <- function(matrices) {
  starts <- lapply(matrices, function(x) {
    if (is.null(x)) return(NULL)
    lab <- gsub("^[^*]*\\*", "", x)
    val <- suppressWarnings(as.numeric(gsub("\\*.*$", "", x)))
    ok <- !is.na(val) & lab != x & !startsWith(lab, "data.")
    setNames(as.list(val[ok]), lab[ok])
  })
  do.call(c, unname(starts))
}

## Select the "==" constraints that replace.constraints=TRUE substitutes:
## those whose left-hand side is a free parameter label in "para.labels".
## Only the first constraint for each label is substituted; any further
## constraint on the same label is kept and rewritten through the
## substitution (e.g., "b == 2*a; b == c + 0.3" gives b = 2*a and the
## constraint 2*a == c + 0.3). Returns the resolved replacements (label ->
## RHS text), the remaining mxalgebras with the substitution applied, and
## the labels introduced by the replacements.
.select.substitutions <- function(mxalgebras, para.labels, lbound=NULL, ubound=NULL) {
  lhs <- vapply(mxalgebras, function(x) {
    form <- as.character(x$formula)
    if (length(form) >= 3L && form[1] == "==" && form[2] %in% para.labels)
      form[2] else NA_character_
  }, character(1))
  index <- !is.na(lhs) & !duplicated(lhs, incomparables=NA)
  if (!any(index))
    return(list(replacements=NULL, remaining=mxalgebras,
                introduced=character(0)))
  labels <- unname(lhs[index])
  rhs <- vapply(mxalgebras[index], function(x) as.character(x$formula)[3],
                character(1))
  rhs <- .resolve.chained.constraints(labels, rhs)
  subs <- setNames(lapply(as.list(rhs), function(r) parse(text=r)[[1]]), labels)
  remaining <- .substitute.remaining.mxalgebras(mxalgebras[!index], subs)
  ## A bound belongs to the substituted parameter, not to the remaining
  ## free variables. Enforce it on its replacement expression.
  for (lab in labels) for (side in c("lower", "upper")) {
    bounds <- if (side == "lower") lbound else ubound
    if (is.null(names(bounds)) || !(lab %in% names(bounds))) next
    bound <- bounds[[lab]]
    if (is.null(bound) || is.na(bound)) next
    expr <- subs[[lab]]
    vars <- all.vars(expr)
    if (any(startsWith(vars, "data.")))
      stop("Cannot preserve a bound on substituted parameter '", lab,
           "' with definition variables; use replace.constraints=FALSE.", call.=FALSE)
    if (!length(vars)) {
      value <- tryCatch(eval(expr, envir=baseenv()), error=function(e) NA_real_)
      if (length(value) != 1L || !is.finite(value) ||
          (side == "lower" && value < bound) || (side == "upper" && value > bound))
        stop("Substitution for '", lab, "' violates its ", side, " bound.", call.=FALSE)
      next
    }
    name <- paste0("substituted_", side, "_", lab)
    if (name %in% names(remaining))
      stop("Cannot create a unique bound constraint for '", lab, "'.", call.=FALSE)
    comparison <- call(if (side == "lower") ">" else "<", expr, as.numeric(bound))
    remaining[[name]] <- eval(call("mxConstraint", comparison, name=name))
  }
  if (!length(remaining)) remaining <- NULL
  introduced <- all.vars(parse(text=rhs))
  list(replacements=setNames(unname(rhs), labels), remaining=remaining,
       introduced=introduced[!startsWith(introduced, "data.")])
}

## After replace.constraints=TRUE eliminates a label (folding it into an
## algebra computed from other, "real" free parameters), any REMAINING
## mxAlgebra (a lavaan ":=" defined parameter) or mxConstraint (an
## unconsumed "=="/"<"/">" row) that still references that label by name
## needs the SAME substitution applied to ITS OWN formula, or it fails at
## mxRun() with an "Unknown reference" error the moment the model tries to
## resolve it -- confirmed by testing: "foo := sigma1 + sigma2; sigma1 ==
## exp(a0)" errors this way for "foo" once "sigma1" is fully replaced.
## ":=" rows were never part of the "==" substitution set to begin with,
## in any of sem()'s three RAM dispatch paths (single-group, two-level,
## multiple-group) -- this is called from all three, right after each
## removes its own consumed "==" entries from "mxalgebras".
##
## "mxalgebras" is the REMAINING (already-consumed-entries-removed)
## constraint list; "subs" is the resolved label->expression map (as
## PARSED language objects, not text) from the "==" constraints that WERE
## just replaced. An entry that references none of "subs"'s labels is
## returned completely untouched (same object, not just an equivalent
## rebuild) -- cheap to check via all.vars(), and avoids needlessly
## rebuilding the vast majority of models' constraints, which reference
## no replaced label at all.
##
## mxAlgebra/mxConstraint objects are rebuilt via deparse-then-eval,
## mirroring lavaan2RAM()'s own construction of these exact objects
## (Amatrix/Smatrix's ":="/"==" rows) -- both take their first argument as
## a single whole-expression formula (for mxConstraint, INCLUDING the
## comparison operator, e.g. "sigma1 == exp(a0)", not separate lhs/op/rhs
## arguments), confirmed directly against mxAlgebra()/mxConstraint()'s own
## $formula slot. A ":="-defined name is distinguished from an unconsumed
## "=="/"<"/">" row by its class (MxConstraint), not by its name: a defined
## parameter may itself be named like "constraint1_effect".
.substitute.remaining.mxalgebras <- function(mxalgebras, subs) {
  if (is.null(mxalgebras) || length(subs) == 0) return(mxalgebras)

  for (nm in names(mxalgebras)) {
    formula <- mxalgebras[[nm]]$formula
    if (length(intersect(all.vars(formula), names(subs))) == 0) next

    formula <- .substitute.symbols(formula, subs)
    new.text <- deparse1(formula)
    if (inherits(mxalgebras[[nm]], "MxConstraint")) {
      ## A constraint that the substitution makes identical on both sides
      ## (e.g., a repeated "b == 2*a") always holds and is dropped.
      if (length(formula) == 3L && identical(formula[[2]], formula[[3]])) {
        mxalgebras[[nm]] <- NULL
        next
      }
      mxalgebras[[nm]] <- eval(parse(text=paste0(
        "mxConstraint(", new.text, ", name='", nm, "')")))
    } else {
      mxalgebras[[nm]] <- eval(parse(text=paste0(
        "mxAlgebra(", new.text, ", name='", nm, "')")))
    }
  }
  mxalgebras
}

## Apply one A/S/M-matrix "==" constraint substitution (replace.constraints
## =TRUE) to "model", which is either the top ("within") model or the
## "between" submodel of a two-level fit: remove the auto-built matrix
## named "mat.name" and re-add either a plain fixed mxMatrix (if every
## free symbol was substituted away, leaving a constant) or an mxAlgebra
## decomposition (if free parameters remain) -- the same unchanged/
## constant/algebra branching sem()'s single-group path already uses.
## Returns "model" untouched if "new.sym" is identical to "orig.sym"
## (nothing to replace in this particular matrix).
##
## "model" here is always a type="RAM" model whose A/S/F/M matrices were
## auto-built by mxModel() from mxPath() calls (see .RAM2paths() above).
## Unlike sem()'s single-group path -- which passes dimnames=var.names
## directly to mxExpectationRAM() -- that auto-built expectation relies on
## the matrix OBJECTS themselves carrying dimnames, which neither
## as.mxAlgebra()'s result nor a plain numeric matrix has by default;
## both are given "new.sym"'s dimnames explicitly below (confirmed
## necessary by testing: an M-matrix replacement lacking them fails with
## "does not contain dimnames").
.twolevel.replace.one <- function(model, mat.name, orig.sym, new.sym,
                                  startvalues, lbound, ubound) {
  if (all(orig.sym==new.sym)) return(model)

  model <- mxModel(model, mat.name, remove=TRUE)

  if (length(all.vars(parse(text=new.sym))) == 0) {
    ## Fully constant after substitution (e.g. a constraint fixing a
    ## parameter at a literal number, or at a numeric expression like
    ## "exp(0.4)"). as.mxMatrix() cannot evaluate a non-trivial constant
    ## EXPRESSION -- it only recognises bare numbers or "value*label"
    ## strings -- so the expression is evaluated by hand first and a
    ## plain numeric matrix is passed instead.
    values <- matrix(vapply(new.sym, function(z) eval(parse(text=z)), numeric(1)),
                     nrow=nrow(new.sym), ncol=ncol(new.sym), dimnames=dimnames(new.sym))
    mat <- as.mxMatrix(values, name=mat.name)
    ## as.mxMatrix() silently SKIPS setting dimnames on the mxMatrix it
    ## builds whenever rownames(x)[1]=="1" (its own convention for
    ## identifying an M/mean-vector row -- see as.mxMatrix.R) -- which is
    ## ALWAYS true for an M-matrix replacement, even though "values" above
    ## was itself built with dimnames. Confirmed by testing: without this,
    ## a constant M replacement (e.g. "ymean_between == 1") reaches
    ## mxRun() as "the M matrix ... does not contain dimnames", exactly the
    ## failure the algebra branch below already guards against. Set
    ## unconditionally here too, for consistency with that branch (a no-op
    ## for A/S, which as.mxMatrix() already dimnames correctly).
    dimnames(mat) <- dimnames(new.sym)
    mxModel(model, mat)
  } else {
    alg <- as.mxAlgebra(new.sym, startvalues=startvalues, lbound=lbound,
                        ubound=ubound, name=mat.name)
    dimnames(alg$mxalgebra) <- dimnames(new.sym)
    mxModel(model, alg$mxalgebra, alg$parameters, alg$list)
  }
}

## Fit a two-level RAM object using OpenMx's relational SEM mechanism: a
## level-2 ("between") mxModel(type="RAM") with mxData(primaryKey=cluster),
## nested inside a level-1 ("within") mxModel(type="RAM") that links to it
## via fixed unit-loading mxPath(joinKey=cluster) paths -- one per "shadow"
## variable (a between-level variable that is also a level-1 manifest
## variable). See dev/technical.md for the current design and
## dev/debug-two-level-multigroup/03-04 for the original comparison.
.sem.twolevel <- function(model.name="sem", RAM, data=NULL, cluster=NULL,
                          intervals.type=c("z","LB"), startvalues=NULL,
                          lbound=NULL, ubound=NULL, replace.constraints=FALSE,
                          mxModel.Args=NULL, run=TRUE, silent=TRUE,
                          between.name="between", replacements=NULL, ...) {

  intervals.type <- match.arg(intervals.type)

  if (is.null(data) || !is.data.frame(data)) {
    stop("Two-level models require either a single data frame of ",
        "individual-level (level-1) rows passed via 'data', or per-cluster ",
        "'Cov', 'means', and 'numObs', plus 'cluster' naming the column ",
        "that identifies level-2 (cluster) membership.")
  }
  if (is.null(cluster) || !(cluster %in% colnames(data))) {
    stop("'cluster' must name a column in 'data' identifying level-2 ",
        "(cluster) membership for two-level models.")
  }
  ## Left unguarded, unique(data[[cluster]]) below would include NA as a
  ## cluster in its own right, creating a spurious "phantom" between-level
  ## row (with primaryKey NA) that rows with a missing cluster ID silently
  ## join to -- found by testing, not something OpenMx itself rejects.
  if (anyNA(data[[cluster]])) {
    stop("'data[[\"", cluster, "\"]]' contains missing values; two-level ",
        "models require every row to have a known cluster membership. ",
        "Remove or impute rows with a missing cluster ID before calling sem().")
  }
  ## OpenMx's primaryKey/joinKey mechanism requires the key column to be
  ## an integer or a factor -- a character ID (e.g. "school-1") or a plain
  ## numeric/double ID (e.g. as.numeric(1:10)), both very ordinary ways to
  ## code a cluster ID in real data, are rejected deep inside mxRun() with
  ## "primary key must be an integer or factor column in raw observed
  ## data". Because that error is a hard R error (not a soft non-zero
  ## status), sem()'s own mxRun-then-mxTryHard fallback treats it as
  ## "maybe just needs a better start" and burns through 50 retries before
  ## finally re-raising it -- confirmed empirically, a real usability trap,
  ## not just a slower error. Converting to a factor up front (used for
  ## both the between-level primaryKey data and the within-level data
  ## below, so the two stay consistent) avoids ever reaching mxRun() with
  ## an unsupported key type in the first place.
  if (!is.integer(data[[cluster]]) && !is.factor(data[[cluster]])) {
    data[[cluster]] <- factor(data[[cluster]], levels=unique(data[[cluster]]))
  }

  checkRAM(as.mxMatrix(RAM$within$A), as.mxMatrix(RAM$within$S), cor.analysis=FALSE)
  checkRAM(as.mxMatrix(RAM$between$A), as.mxMatrix(RAM$between$S), cor.analysis=FALSE)

  manifest.vars <- rownames(RAM$within$F)
  shadow.vars <- intersect(manifest.vars, colnames(RAM$between$A))
  within.RAM <- .decompose.twolevel.within(RAM$within, shadow.vars)

  ## ---- Level 2 ("between"): one row per cluster ----
  between.build <- .RAM2paths(RAM$between$A, RAM$between$S, RAM$between$M,
                              startvalues=startvalues, lbound=lbound, ubound=ubound)
  between.manifest <- rownames(RAM$between$F)
  missing.between <- setdiff(between.manifest, colnames(data))
  if (length(missing.between))
    stop("Between-level observed variable(s) missing from 'data': ",
         paste(missing.between, collapse=", "))
  for (v in between.manifest) {
    if (!is.numeric(data[[v]]))
      stop("Between-level observed variable '", v, "' must be numeric.")
    by.cluster <- split(data[[v]], data[[cluster]], drop=TRUE)
    if (any(vapply(by.cluster, function(x) length(unique(x)) > 1L, logical(1))))
      stop("Between-level observed variable '", v,
           "' must be constant within each cluster.")
  }
  between.data <- data[!duplicated(data[[cluster]]),
                       c(cluster, between.manifest), drop=FALSE]
  between.model <- mxModel(between.name, type="RAM",
                           manifestVars=between.manifest,
                           latentVars=setdiff(colnames(RAM$between$A), between.manifest),
                           mxData(between.data, type="raw", primaryKey=cluster),
                           between.build$paths)

  ## ---- Level 1 ("within"): individual rows, nests "between" ----
  within.build <- .RAM2paths(within.RAM$A, within.RAM$S, within.RAM$M,
                             startvalues=startvalues, lbound=lbound, ubound=ubound)
  latent.vars <- setdiff(colnames(within.RAM$A), manifest.vars)

  join.paths <- NULL
  if (length(shadow.vars) > 0) {
    join.paths <- mxPath(from=paste0(between.name, ".", shadow.vars), to=shadow.vars,
                         free=FALSE, values=1, joinKey=cluster)
  }

  ## ---- replace.constraints=TRUE: substitute "==" constraints whose LHS
  ## is a genuine free parameter at either level directly into the
  ## within/between A/S/M matrices, mirroring sem()'s single-group
  ## substitution logic (see as.symMatrix()/as.mxAlgebra() there), but
  ## applied across both levels -- including a constraint that spans them
  ## (e.g. equating a between-level parameter to a function of a
  ## within-level one; label resolution across the two submodels needs no
  ## special handling, since OpenMx labels are already global across a
  ## model tree). mxalgebras is mutated here (consumed constraints
  ## removed) so the ordinary "==" / ":=" loop further below only adds
  ## whatever was NOT replaced.
  mxalgebras <- RAM$mxalgebras
  startvalues2 <- startvalues
  Aw <- Aw0 <- as.symMatrix(within.RAM$A)
  Sw <- Sw0 <- as.symMatrix(within.RAM$S)
  Mw <- Mw0 <- as.symMatrix(within.RAM$M)
  Ab <- Ab0 <- as.symMatrix(RAM$between$A)
  Sb <- Sb0 <- as.symMatrix(RAM$between$S)
  Mb <- Mb0 <- as.symMatrix(RAM$between$M)

  ## "replacements" (label -> resolved RHS text) is given by
  ## .sem.multigroup.twolevel(), which selects the substitutions once for
  ## all groups; otherwise they are selected here.
  if (isTRUE(replace.constraints) && is.null(replacements)) {
    selection <- .select.substitutions(
      mxalgebras, unique(c(within.build$para.labels, between.build$para.labels)), lbound, ubound)
    replacements <- selection$replacements
    mxalgebras <- selection$remaining
  }
  if (!isTRUE(replace.constraints)) replacements <- NULL

  if (length(replacements)) {
    Aw <- .apply.replacements(Aw, replacements)
    Sw <- .apply.replacements(Sw, replacements)
    Mw <- .apply.replacements(Mw, replacements)
    Ab <- .apply.replacements(Ab, replacements)
    Sb <- .apply.replacements(Sb, replacements)
    Mb <- .apply.replacements(Mb, replacements)

    ## Starting values of both levels, for labels used across levels (see
    ## .ram.starts()); those in "startvalues" take precedence.
    existing.starts <- .ram.starts(list(RAM$within$A, RAM$within$S,
                                        RAM$within$M, RAM$between$A,
                                        RAM$between$S, RAM$between$M))
    startvalues2 <- c(existing.starts, startvalues)
    startvalues2 <- startvalues2[!duplicated(names(startvalues2), fromLast=TRUE)]

    ## Between-level replacements must land on "between.model" BEFORE it
    ## is nested into "mx.model" below -- mxModel(model, "S",
    ## remove=TRUE) only ever reaches into "model"'s own direct children,
    ## not a submodel's.
    between.model <- .twolevel.replace.one(between.model, "A", Ab0, Ab,
                                           startvalues2, lbound, ubound)
    between.model <- .twolevel.replace.one(between.model, "S", Sb0, Sb,
                                           startvalues2, lbound, ubound)
    between.model <- .twolevel.replace.one(between.model, "M", Mb0, Mb,
                                           startvalues2, lbound, ubound)
  }

  ## Add cluster-constant definition variables used at the between level,
  ## including references introduced by replace.constraints substitution.
  between.vars <- all.vars(parse(text=c(Ab, Sb, Mb)))
  between.defvars <- between.vars[startsWith(between.vars, "data.")]
  def.cols <- setdiff(sub("^data\\.", "", between.defvars), cluster)
  missing.cols <- setdiff(def.cols, colnames(data))
  if (length(missing.cols) > 0) {
    stop("Between-level definition variable column(s) missing from 'data': ",
         paste(missing.cols, collapse=", "), call.=FALSE)
  }
  for (v in def.cols) {
    by.cluster <- split(data[[v]], data[[cluster]], drop=TRUE)
    if (any(vapply(by.cluster, function(x) length(unique(x)) > 1L, logical(1))))
      stop("Between-level definition variable '", v,
           "' must be constant within each cluster.")
    between.data[[v]] <- data[[v]][!duplicated(data[[cluster]])]
  }
  if (length(def.cols))
    between.model <- mxModel(between.model,
                             mxData(between.data, type="raw", primaryKey=cluster))

  mx.model <- mxModel(model.name, type="RAM", between.model,
                      manifestVars=manifest.vars, latentVars=latent.vars,
                      mxData(data, type="raw"),
                      within.build$paths, join.paths)

  ## Within-level replacements land on "mx.model" itself (the top model
  ## IS the "within" level here), so this can only happen after it exists.
  if (length(replacements)) {
    mx.model <- .twolevel.replace.one(mx.model, "A", Aw0, Aw, startvalues2, lbound, ubound)
    mx.model <- .twolevel.replace.one(mx.model, "S", Sw0, Sw, startvalues2, lbound, ubound)
    mx.model <- .twolevel.replace.one(mx.model, "M", Mw0, Mw, startvalues2, lbound, ubound)
  }

  ## mxCI over the union of free-parameter labels across both levels. When
  ## replace.constraints applied, some ORIGINAL labels (e.g. a parameter
  ## that got replaced by an algebra) are no longer free parameters at
  ## all, so this is recomputed from the (possibly substituted) matrices
  ## instead -- exactly mirroring how sem()'s single-group path computes
  ## its own "new.para.labels" (all.vars() over the substituted matrices,
  ## excluding "data." definition-variable references).
  if (length(replacements)) {
    all.para.labels <- unique(c(Aw, Sw, Mw, Ab, Sb, Mb))
    all.para.labels <- all.vars(parse(text=all.para.labels))
    all.para.labels <- all.para.labels[!startsWith(all.para.labels, "data.")]
  } else {
    all.para.labels <- unique(c(between.build$para.labels, within.build$para.labels))
  }
  ## mxCI() rejects an empty reference vector outright rather than treating
  ## it as "no CIs requested" (confirmed by testing) -- reachable whenever
  ## every within- and between-level free parameter has been fixed via
  ## replace.constraints=TRUE, a valid, if unusual, model to fit.
  if (length(all.para.labels) > 0) {
    mx.model <- mxModel(mx.model, mxCI(all.para.labels))
  }

  ## mxAlgebra/mxConstraint from lavaan2RAM()'s level==0 rows (fn1 := ...,
  ## a == b, etc.), attached at the container level -- see
  ## .lavaan2RAM.twolevel()'s comment on why these live at RAM$mxalgebras
  ## directly (not nested under "within"/"between"). Constraints already
  ## consumed by the replace.constraints substitution above have already
  ## been removed from "mxalgebras" at this point.
  mxalgebras.ci <- NULL
  if (!is.null(mxalgebras)) {
    for (i in seq_along(mxalgebras)) {
      mx.model <- mxModel(mx.model, mxalgebras[[i]])
      name.mxalgebra <- names(mxalgebras)[i]
      if (!inherits(mxalgebras[[i]], "MxConstraint")) {
        mx.model <- mxModel(mx.model, mxCI(c(name.mxalgebra)))
        mxalgebras.ci <- c(mxalgebras.ci, name.mxalgebra)
      }
    }
  }

  if (!is.null(mxModel.Args)) {
    for (i in seq_along(mxModel.Args)) {
      mx.model <- mxModel(mx.model, mxModel.Args[[i]])
    }
  }

  warning.msg <- error.msg <- NULL

  if (run) {
    fallback <- .run.mxTryHard.fallback(mx.model, intervals.type, silent=silent, ...)
    mx.fit <- fallback$mx.fit
    warning.msg <- fallback$warning.msg
    error.msg <- fallback$error.msg
  } else {
    mx.fit <- .tighten.feasibility(mx.model)
  }

  out <- list(mx.fit=mx.fit, RAM=RAM, data=data, mxalgebras=mxalgebras.ci,
              intervals.type=intervals.type, warning=warning.msg,
              error=error.msg,
              cluster=cluster,
              within.aliases=attr(within.RAM, "shadow.aliases"))
  class(out) <- "mxsem"
  out
}


## Equality-only symbols are not automatically free OpenMx parameters.
.check.replacement.parameters <- function(RAM, replace.constraints, extras=NULL) {
  if (isTRUE(replace.constraints) || is.null(RAM)) return(invisible(NULL))
  extra.labels <- function(x, path=character(0)) {
    qualify <- function(name) if (!length(path)) name else
      paste(paste(path,collapse="."),name,sep=".")
    if (inherits(x,"MxMatrix"))
      return(c(x$labels,qualify(x$name),qualify(c(x$labels))))
    if (inherits(x,c("MxAlgebra","MxConstraint"))) return(qualify(x$name))
    if (inherits(x,"MxModel")) {
      nested <- c(path,x$name)
      return(c(qualify(x$name),unlist(lapply(x@matrices,extra.labels,path=nested)),
               unlist(lapply(x@algebras,extra.labels,path=nested)),
               unlist(lapply(x@submodels,extra.labels,path=nested))))
    }
    if (is.list(x)) return(unlist(lapply(x,extra.labels,path=path)))
    NULL
  }
  collect <- function(x) {
    if (!is.list(x)) return(list(matrices=list(), algebras=list()))
    if (all(c("A", "S", "F") %in% names(x)))
      return(list(matrices=Filter(Negate(is.null), x[c("A", "S", "M")]),
                  algebras=x$mxalgebras))
    parts <- lapply(x, collect)
    list(matrices=unlist(lapply(parts, `[[`, "matrices"), recursive=FALSE),
         algebras=c(attr(x, "mxalgebras"), x$mxalgebras,
                    unlist(lapply(parts, `[[`, "algebras"), recursive=FALSE)))
  }
  parts <- collect(RAM)
  if (!length(parts$algebras)) return(invisible(NULL))
  matrices <- lapply(parts$matrices,function(x) as.mxMatrix(x,name="check"))
  labels <- unique(c(unlist(lapply(matrices,function(m) c(m$labels))),extra.labels(extras)))
  introduced <- unique(unlist(lapply(parts$algebras,function(x) {
    e <- x$formula
    if (inherits(x,"MxConstraint") && identical(e[[1]],as.name("==")) &&
        deparse1(e[[2]]) %in% labels) all.vars(e[[3]]) else NULL
  })))
  unknown <- setdiff(as.character(introduced),c(labels,names(parts$algebras),
    "Amatrix","Smatrix","Mmatrix","Fmatrix","Id"))
  ## Group model components can be referenced as e.g. B.Amatrix[1,2].
  ## Their final group names are assigned by the dispatchers; let OpenMx
  ## resolve these qualified references instead of treating them as labels.
  components <- c("Amatrix", "Smatrix", "Mmatrix", "Fmatrix", "Id", "Id_A", "expCov", "expMean")
  qualified.component <- vapply(strsplit(unknown, ".", fixed=TRUE), function(parts)
    length(parts) > 1L && parts[length(parts)] %in% components, logical(1))
  unknown <- unknown[!startsWith(unknown,"data.") & !qualified.component]
  if (length(unknown)) stop("Parameters introduced only in equality expressions (",
    paste(unknown,collapse=", "),") require 'replace.constraints=TRUE' or an explicit parameter matrix.",
    call.=FALSE)
  invisible(NULL)
}
