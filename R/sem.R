## NOTE (previously a FIXME here): replace.constraints=TRUE combined with
## intervals.type="LB" was originally believed to error for single-group
## models, with a matching @note on sem() itself. Retested after this
## session's chained-constraint fix (.resolve.chained.constraints()) --
## which is the likeliest original cause, since an unresolved chain used
## to leave an orphaned, disconnected free parameter whose name would not
## have matched what mxCI() was told to expect -- and confirmed working
## correctly now, including with a chained constraint. See sem()'s @note
## for the currently-accurate wording.


#' Fit a structural equation model using OpenMx
#' 
#' It fits a structural equation model by creating a mxModel from a RAM object.
#' 
#' 
#' @aliases sem create.mxModel
#' @param model.name A string for the model name in
#' \code{\link[OpenMx]{mxModel}}.
#' @param RAM A RAM object including a list of matrices of the model returned
#' from \code{\link[metaSEM]{lavaan2RAM}}.
#' @param data A data frame or matrix of data.
#' @param Cov A covariance matrix may also be used if \code{data}==NULL.
#' @param means A named vector of means (optional) if \code{Cov} is used.
#' @param numObs If \code{Cov} is used, a sample size must be provided.
#' @param intervals.type Either \code{z} (default if missing) or \code{LB}. If
#' it is \code{z}, it calculates the 95\% confidence intervals (CIs) based on
#' the estimated standard error. If it is \code{LB}, it calculates the 95\%
#' likelihood-based CIs on the parameter estimates.
#' @param startvalues A list of named starting values of the free parameters,
#' e.g., list(a=1, b=2)
#' @param lbound A list of lower bound of the free parameters. If it is not
#' provided, all free parameters are assumed \code{NA}.
#' @param ubound A list of upper bound of the free parameters. If it is not
#' provided, all free parameters are assumed \code{NA}.
#' @param replace.constraints Logical. If \code{TRUE}, the parameters on the
#' left hand side will be replaced by the constraints on the right hand side.
#' That is, the parameters on the left hand side are no longer parameters in
#' the model.
#' @param mxModel.Args A list of arguments passed to
#' \code{\link[OpenMx]{mxModel}}.
#' @param run Logical. If \code{FALSE}, only return the mx model without
#' running the analysis.
#' @param silent Logical. An argument is passed to either
#' \code{\link[OpenMx]{mxRun}} or \code{\link[OpenMx]{mxTryHard}}
#' @param group A character string naming the column in \code{data}
#' identifying group membership. Required (and only used) when \code{RAM} is
#' a multiple-group RAM object, i.e., created by \code{\link[metaSEM]{lavaan2RAM}}
#' with \code{ngroups > 1}. \code{data} must be a single data frame in this
#' case; it is split internally by the unique values of \code{data[[group]]}
#' in \emph{first-appearance order in the data} (not sorted, and not a
#' factor's own \code{levels()} order, even for a factor column -- this
#' matches how \code{lavaan} itself numbers groups), which must correspond,
#' in order, to the groups used when the model was specified (e.g.,
#' \code{lavaan}'s \code{c(...)} syntax numbers groups the same way). A
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
#' mismatch. \code{Cov}/\code{numObs}-based input is not yet supported for
#' multiple-group models. \code{replace.constraints=TRUE} IS supported
#' (substituting a constraint into whichever group's \code{A}/\code{S}/
#' \code{M} matrices contain it, including a constraint that spans
#' multiple groups, e.g. constraining one group's variance to be a
#' multiple of another's); a definition variable (\code{data.*}) inside
#' such a constraint resolves against whichever group's own data it is
#' used in, since (unlike two-level's between level) every group already
#' carries its own full data slice.
#' @param cluster A character string naming the column in \code{data}
#' identifying level-2 (cluster) membership. Required (and only used) when
#' \code{RAM} is a two-level RAM object, i.e., created by
#' \code{\link[metaSEM]{lavaan2RAM}} from \code{lavaan}'s two-level
#' (\code{level: 1} / \code{level: 2}) syntax. \code{data} must be a single
#' data frame of individual-level (level-1) rows in this case.
#' \code{Cov}/\code{numObs}-based input is not yet supported for two-level
#' models, and neither is combining two-level and multiple-group models in
#' the same call. \code{replace.constraints=TRUE} IS supported for
#' two-level models (substituting a constraint into whichever of the
#' within-/between-level \code{A}/\code{S}/\code{M} matrices contains it,
#' including a constraint that spans both levels); a definition variable
#' (\code{data.*}) may be used in such a constraint for a level-1 (within)
#' parameter, but not yet for a genuine level-2-only (between, cluster-level)
#' covariate, since \code{data} passed to the level-2 submodel currently
#' carries only the cluster ID column.
#' @param \dots Further arguments will be passed to either
#' \code{\link[OpenMx]{mxRun}} or \code{\link[OpenMx]{mxTryHard}}
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

  intervals.type <- match.arg(intervals.type)

  ## Multiple-group RAM (from lavaan2RAM(..., ngroups>1)) and two-level RAM
  ## (from lavaan2RAM() on "level: 1"/"level: 2" syntax): dispatch each to
  ## a dedicated builder rather than threading extra code paths through the
  ## single-group/single-level logic below (kept untouched to avoid
  ## regressions).
  if (.is.RAM.multigroup(RAM)) {
    ## Combined two-level + multi-group is not supported (see
    ## lavaan2RAM()); a "cluster" argument here could not possibly do
    ## anything, so reject it explicitly rather than silently ignore it --
    ## passing both could otherwise mislead a caller into thinking both
    ## dimensions were modeled when only "group" was.
    if (!is.null(cluster)) {
      stop("'cluster' has no effect for a multiple-group RAM object ",
          "(class \"RAM_multigroup\"); combined two-level and ",
          "multiple-group models are not supported. Did you mean 'group'?")
    }
    return(.sem.multigroup(model.name=model.name, RAM=RAM, data=data,
                           group=group, intervals.type=intervals.type,
                           startvalues=startvalues, lbound=lbound, ubound=ubound,
                           replace.constraints=replace.constraints,
                           mxModel.Args=mxModel.Args, run=run, silent=silent, ...))
  }
  if (.is.RAM.twolevel(RAM)) {
    if (!is.null(group)) {
      stop("'group' has no effect for a two-level RAM object (class ",
          "\"RAM_twolevel\"); combined two-level and multiple-group ",
          "models are not supported. Did you mean 'cluster'?")
    }
    return(.sem.twolevel(model.name=model.name, RAM=RAM, data=data,
                         cluster=cluster, intervals.type=intervals.type,
                         startvalues=startvalues, lbound=lbound, ubound=ubound,
                         replace.constraints=replace.constraints,
                         mxModel.Args=mxModel.Args, run=run, silent=silent, ...))
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
  
  ## Any names of the constraints == parameters?
  ## If yes, these parameters are replaced by the constraints
  index <- sapply(mxalgebras, function(x) {
    ## Convert R language to a vector string
    ## form[1]: "=="
    ## form[2]: "m"
    ## form[3]: "p1 * cos(p2 * data.x) + p2 * sin(p1 * data.x)"
    form <- as.character(x$formula)
    if (form[1]=="==" & form[2]%in% para.labels) TRUE else FALSE
  })
    
  ############################################# 
  ## Need to replace parameters with mxalgebras,
  ## if any TRUE and replace.constraints==TRUE
  if (any(index) & replace.constraints) {

    ## Extract constraints that needed to be replaced
    mxalgebras.const <- mxalgebras[index]

    ## Resolve any constraint whose RHS references another constraint's
    ## OWN label to a fixed point BEFORE applying any of them -- otherwise
    ## a chain (e.g. "residW == exp(a0); residB == 2*residW") would apply
    ## each constraint's un-expanded RHS text independently and silently
    ## reintroduce the eliminated label as a brand-new, disconnected free
    ## parameter (see .resolve.chained.constraints()).
    const.labels <- vapply(mxalgebras.const, function(x) as.character(x$formula)[2],
                           character(1))
    const.rhs <- vapply(mxalgebras.const, function(x) as.character(x$formula)[3],
                        character(1))
    const.rhs <- .resolve.chained.constraints(const.labels, const.rhs)

    for (i in seq_along(const.labels)) {
      lab <- const.labels[i]
      rhs <- const.rhs[i]

      ## Replace the A matrix
      if (any(grep(lab, A))) {
        A[which(lab==A)] <- rhs
      }

      ## Replace the S matrix
      if (any(grep(lab, S))) {
        S[which(lab==S)] <- rhs
      }

      ## Replace the M matrix
      if (any(grep(lab, M))) {
        M[which(lab==M)] <- rhs
      }
    }

    ## Remove the constraints so they won't be added again
    mxalgebras[index] <- NULL

    ## Any REMAINING ":="/"=="/"<"/">" entry that still references one of
    ## the just-eliminated labels needs the same substitution applied to
    ## its own formula (see .substitute.remaining.mxalgebras()).
    subs <- setNames(lapply(as.list(const.rhs), function(r) parse(text=r)[[1]]),
                     const.labels)
    mxalgebras <- .substitute.remaining.mxalgebras(mxalgebras, subs)
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
  new.para.labels <- unique(c(A, S, M))
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
      ## Check if the name constains constraint1, constraint2, ...,
      ## If no, they are mxalgebra, not mxconstraints. Include them in mxCI. 
      if (!grepl("^constraint[0-9]", name.mxalgebra)) {
        mx.model <- mxModel(mx.model, mxCI(c(name.mxalgebra)))
        mxalgebras.ci <- c(mxalgebras.ci, name.mxalgebra)
      }
    }
  }
    
  ## Initialise here so they are available regardless of whether run=TRUE
  warning.msg <- error.msg <- NULL

  if (run) {
    fallback <- .run.mxTryHard.fallback(mx.model, intervals.type, ...)
    mx.fit <- fallback$mx.fit
    warning.msg <- fallback$warning.msg
    error.msg <- fallback$error.msg
  } else {
    mx.fit <- mx.model
  }

  out <- list(mx.fit=mx.fit, RAM=RAM, data=data, mxalgebras=mxalgebras.ci,
              intervals.type=intervals.type, warning=warning.msg,
              error=error.msg)
  class(out) <- "mxsem"
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


summary.mxsem <- function(object, robust=FALSE, ...) {
  if (!is.element("mxsem", class(object)))
    stop("\"object\" must be an object of class \"mxsem\".")
  .check.mx.fit(object)

  ## imxRobustSE() cannot compute row-wise gradients through a submodel
  ## unless the top model's fit function is mxFitFunctionMultigroup -- true
  ## for the multi-group path, but not for the two-level (relational,
  ## primaryKey/joinKey) path, where "between" is a submodel of a plain
  ## mxFitFunctionML() model. Caught here with a clearer message instead of
  ## letting OpenMx's own ("...please use an MxFitFunctionMultigroup...")
  ## error, which doesn't mention two-level models at all, surface as-is.
  if (isTRUE(robust) && .is.RAM.twolevel(object$RAM)) {
    stop("robust=TRUE (robust/sandwich standard errors) is not supported ",
        "for two-level 'sem' fits: OpenMx's imxRobustSE() cannot compute ",
        "row-wise gradients through the relational (primaryKey/joinKey) ",
        "'between' submodel.")
  }

  # calculate coefficients    
  my.mx <- summary(object$mx.fit)
  ## Exclude lbound ubound etc
  my.para <- my.mx$parameters[, 1:6, drop=FALSE]   

  # Determine if CIs on parameter estimates are present
  if (object$intervals.type=="z") {

    ## Replace the SEs with robust SEs
    if (robust) {
      my.robust <- suppressMessages(imxRobustSE(object$mx.fit))
      my.para[, "Std.Error"] <- my.robust[my.para$name]
    }
        
    my.para$lbound <- with(my.para, Estimate - qnorm(.975)*Std.Error)
    my.para$ubound <- with(my.para, Estimate + qnorm(.975)*Std.Error)
    coefficients <- my.para[, -c(1:4), drop=FALSE]
    dimnames(coefficients)[[1]] <- my.para$name
      
  } else {
 
    ## Convert a data frame with length of 0 in my.mx$CI and remove the last column "note"
    my.ci <- my.mx$CI
    if (length(my.ci)==0) my.ci <- NULL else my.ci <- my.ci[, 1:3, drop=FALSE]        
    
    ## Select the elements matched my.para (excluded I2)  
    my.ci <- my.ci[row.names(my.ci) %in% my.para$name, ]
      
    my.ci <- data.frame(name=row.names(my.ci), my.ci)
    my.para <- merge(my.para, my.ci, by=c("name"))      
    coefficients <- my.para[, -c(1:4,8)]
    dimnames(coefficients)[[1]] <- my.para$name
    # NA for LBCI
    coefficients$Std.Error <- NA
  }
  
  coefficients$"z value" <- coefficients$Estimate/coefficients$Std.Error
  coefficients$"Pr(>|z|)" <- 2*(1-pnorm(abs(coefficients$"z value")))

  informationCriteria <- my.mx$informationCriteria
  ## Better column names
  colnames(informationCriteria) <- c("df Penalty", "Parameters Penalty",
                                     "Sample-Size Adjusted")

  ## Get the mxalgebras
  if (!is.null(object$mxalgebras)) {
    if (object$intervals.type=="z") {
      estimate <- eval(parse(text=paste0("mxEval(rbind(",
                                         paste(object$mxalgebras, collapse=","),
                                         "), object$mx.fit)")))
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
        robust.cov <- suppressMessages(imxRobustSE(object$mx.fit, details=TRUE))$cov
        SE <- eval(parse(text=paste0("mxSE(rbind(",
                                     paste(object$mxalgebras, collapse=","),
                                     "), model=object$mx.fit, cov=robust.cov, silent=TRUE)")))
      } else {
        SE <- eval(parse(text=paste0("mxSE(rbind(",
                                     paste(object$mxalgebras, collapse=","),
                                     "), model=object$mx.fit, silent=TRUE)")))
      }
      mxalgebras <- cbind(lbound=estimate - 1.96*SE,
                          estimate=estimate,
                          ubound=estimate + 1.96*SE)
      dimnames(mxalgebras) <- list(object$mxalgebras,
                                   c("lbound", "estimate", "ubound"))
    } else {
      my.ci <- my.mx$CI
      index <- NULL
      for (i in seq_along(object$mxalgebras)) {
        ## Get the names of the mxalgebras combined with the model name.
        ## fixed=TRUE: the same "." -as-regex-wildcard mistake found
        ## elsewhere in this function -- without it, a CI row name that
        ## merely resembles "<model.name>.<algebra>" (any single
        ## character in place of the literal ".") would also match.
        index <- c(index, grep(paste(object$mx.fit$name, object$mxalgebras[i], sep="."),
                               rownames(my.mx$CI), fixed=TRUE))
      }
      mxalgebras <- my.mx$CI[index, c("lbound", "estimate", "ubound")]
      dimnames(mxalgebras) <- list(object$mxalgebras,
                                   c("lbound", "estimate", "ubound"))
    }
  } else {
    mxalgebras <- NULL
  }

  out <- list(coefficients=coefficients, mxalgebras=mxalgebras,
              intervals.type=object$intervals.type,
              robust=robust, no.studies=my.mx$numObs,
              obsStat=my.mx$observedStatistics,
              estPara=my.mx$estimatedParameters, df=my.mx$degreesOfFreedom,
              Minus2LL=my.mx$Minus2LogLikelihood,
              Mx.status1=object$mx.fit@output$status[[1]],
              informationCriteria=informationCriteria)
    class(out) <- "summary.mxsem"
    out
}

print.summary.mxsem <- function(x, ...) {
    if (!is.element("summary.mxsem", class(x))) {
      stop("\"x\" must be an object of class \"summary.mxsem\".")
    }
    
    cat("95% confidence intervals: ")
    switch(x$intervals.type,
           z = cat("z statistic approximation (robust=", x$robust, ")", sep=""),
           LB = cat("Likelihood-based statistic") )

    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients, P.values=TRUE, ...)

    if (!is.null(x$mxalgebras)) {
      cat("\nMxalgebras:\n")
      print(x$mxalgebras)
    }
    
    cat("\nInformation Criteria:\n")
    print(x$informationCriteria)
   
    cat("\nNumber of subjects (or studies):", x$no.studies)
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

vcov.mxsem <- function(object, robust=FALSE, ...) {
  if (!is.element("mxsem", class(object)))
    stop("\"object\" must be an object of class \"mxsem\".")
  .check.mx.fit(object)

  ## Same guard as summary.mxsem(robust=TRUE) -- see its comment for why
  ## imxRobustSE() cannot handle the two-level (relational) path. Originally
  ## only added there; missing here was an oversight caught on review.
  if (isTRUE(robust) && .is.RAM.twolevel(object$RAM)) {
    stop("robust=TRUE (robust/sandwich standard errors) is not supported ",
        "for two-level 'sem' fits: OpenMx's imxRobustSE() cannot compute ",
        "row-wise gradients through the relational (primaryKey/joinKey) ",
        "'between' submodel.")
  }

  if (robust) {
    suppressMessages(imxRobustSE(object$mx.fit, details=TRUE)$cov)
  } else {
    vcov(object$mx.fit)
  } 
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

  ## combine=TRUE: plot every submodel (all groups, or within+between) as
  ## panels of one figure via par(mfrow=), the same manual approach
  ## semPlot's own examples use for multi-group lavaan fits (semPaths()'s
  ## panelGroups=TRUE only recognises objects it builds internally, not
  ## semPlotModel objects built here via ramModel()).
  if (isTRUE(combine) && .is.RAM.multigroup(RAM)) {
    sub.names <- names(x$mx.fit@submodels)
    ## Title each panel with the user's own group value (e.g. "site-2"),
    ## not the sanitized submodel name (e.g. "site_2") -- translated back
    ## via group.map, which .sem.multigroup() stores for exactly this.
    panel.titles <- sub.names
    if (!is.null(x$group.map)) {
      orig <- names(x$group.map)[match(sub.names, x$group.map)]
      panel.titles <- ifelse(is.na(orig), sub.names, orig)
    }
    op <- graphics::par(mfrow=c(1, length(sub.names)))
    on.exit(graphics::par(op))
    for (i in seq_along(sub.names)) {
      .plot.mxsem.one(x$mx.fit[[sub.names[i]]], RAM[[i]], x$data, character(0),
                      what=what, nCharNodes=nCharNodes, nCharEdges=nCharEdges,
                      layout=layout, sizeMan=sizeMan, sizeLat=sizeLat,
                      edge.label.cex=edge.label.cex, color=color,
                      weighted=weighted, main=panel.titles[i],
                      main.line=main.line, main.adj=main.adj, ...)
    }
    return(invisible(NULL))
  }
  if (isTRUE(combine) && .is.RAM.twolevel(RAM)) {
    op <- graphics::par(mfrow=c(1, 2))
    on.exit(graphics::par(op))
    shadow.vars <- intersect(colnames(RAM$between$A), rownames(RAM$within$F))
    .plot.mxsem.one(x$mx.fit, RAM$within, x$data, character(0),
                    what=what, nCharNodes=nCharNodes, nCharEdges=nCharEdges,
                    layout=layout, sizeMan=sizeMan, sizeLat=sizeLat,
                    edge.label.cex=edge.label.cex, color=color,
                    weighted=weighted, main="within",
                    main.line=main.line, main.adj=main.adj, ...)
    .plot.mxsem.one(x$mx.fit$between, RAM$between, x$data, shadow.vars,
                    what=what, nCharNodes=nCharNodes, nCharEdges=nCharEdges,
                    layout=layout, sizeMan=sizeMan, sizeLat=sizeLat,
                    edge.label.cex=edge.label.cex, color=color,
                    weighted=weighted, main="between",
                    main.line=main.line, main.adj=main.adj, ...)
    return(invisible(NULL))
  }

  mx.sub <- x$mx.fit
  shadow.vars <- character(0)

  ## Multiple-group and two-level mxsem objects hold several submodels, not
  ## one flat RAM -- when not combining, plot ONE submodel's diagram at a
  ## time, picked via "group"/"level" (mirroring how one would look at one
  ## group, or one level, of the full model).
  if (.is.RAM.multigroup(RAM)) {
    ## names(RAM) are just generic group indices ("1", "2", ...); the real,
    ## user-meaningful group names are the fitted submodel names (the
    ## first-appearance-order unique values of data[[group]] used to split
    ## the data in .sem.multigroup(), sanitized into legal mxModel names).
    ## Match "group" against the *original* (pre-sanitizing) group values
    ## via group.map FIRST, falling back to matching sub.names directly
    ## only when group.map doesn't recognise it. This order matters: a
    ## sanitized name can collide with a different group's own literal
    ## value (e.g. original groups "a.b" and "a_b" sanitize/de-duplicate to
    ## submodels "a_b" and "a_b_1" -- matching sub.names first would let
    ## plot(fit, group="a_b") silently select "a.b"'s submodel instead of
    ## the group literally named "a_b", since "a_b" is *also* a valid
    ## submodel name in its own right). Checking group.map first means a
    ## user's own data values always take priority over a coincidental
    ## sanitized-name collision.
    sub.names <- names(x$mx.fit@submodels)
    idx <- .resolve.plot.group(group, sub.names, x$group.map)
    group <- sub.names[idx]
    mx.sub <- x$mx.fit[[group]]
    RAM <- RAM[[idx]]
  } else if (.is.RAM.twolevel(RAM)) {
    if (is.null(level)) level <- "within"
    level <- match.arg(level, c("within", "between"))
    ## The between level has no true manifest variables in this scope (see
    ## dev/PLAN-two-level-multigroup-sem.md section 8) -- semPlot::ramModel()
    ## requires at least one or errors ("invalid 'type' (list) of argument"
    ## for an all-latent model). Draw "shadow" variables (between-level
    ## components that back a within-level manifest variable) as manifest
    ## for this plot only; genuine level-2-only latent factors stay latent.
    shadow.vars <- character(0)
    if (level=="between") {
      mx.sub <- x$mx.fit$between
      shadow.vars <- intersect(colnames(RAM$between$A), rownames(RAM$within$F))
    }
    RAM <- RAM[[level]]
  }

  .plot.mxsem.one(mx.sub, RAM, x$data, shadow.vars, what=what,
                  nCharNodes=nCharNodes, nCharEdges=nCharEdges, layout=layout,
                  sizeMan=sizeMan, sizeLat=sizeLat, edge.label.cex=edge.label.cex,
                  color=color, weighted=weighted, main=NULL, ...)
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
                            main.line=2.5, main.adj=NULL, ...) {

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

  ## A plain mxMatrix's $values already carries dimnames; an mxAlgebra's
  ## $result never does (confirmed by testing) -- fall back to "RAM"'s own
  ## dimnames for that specific matrix, which are unaffected by
  ## replace.constraints and already match "mx.sub" here (RAM is the
  ## caller's per-group/per-level slice, same as mx.sub).
  get.mat <- function(name, dimnames.fallback) {
    if (name %in% names(mx.sub@matrices)) return(mx.sub@matrices[[name]]$values)
    values <- mx.sub@algebras[[name]]$result
    if (!is.null(values) && is.null(dimnames(values))) dimnames(values) <- dimnames.fallback
    values
  }
  A <- get.mat(mat.names[["A"]], dimnames(RAM$A))
  S <- get.mat(mat.names[["S"]], dimnames(RAM$S))
  F <- get.mat(mat.names[["F"]], dimnames(RAM$F))
  M <- get.mat(mat.names[["M"]], dimnames(RAM$M))
  ## Fixed when M is NULL, i.e., no mean structure
  if (is.null(M)) M <- matrix(0, nrow=1, ncol=ncol(A))

  ## When there are definition variables, data in the first role are used in
  ## the output. Better to replace it with their means. "data" is the whole
  ## fitted mxsem object's data (all groups combined, for a multi-group
  ## plot) as an approximation -- the group actually plotted isn't tracked
  ## separately on the mxsem object, so this is a display-only
  ## simplification and does not affect the fitted model.
  for (i in seq_len(nrow(S)))
    for (j in seq_len(ncol(S))) {
      ## fixed=TRUE: matches the literal "data." substring, consistent
      ## with the strsplit() below it -- without this, "." is a regex
      ## wildcard and a literal-looking label such as "database" (data +
      ## b + ase) would also match, then fail the strsplit(fixed=TRUE)
      ## right after it (no literal "data." to split on, so tmp is NA).
      if (grepl("data.", RAM$S[i, j], fixed=TRUE)) {
        tmp <- strsplit(RAM$S[i, j], "data.", fixed=TRUE)[[1]][2]
        S[i, j] <- eval(parse(text=paste0("mean(data$", tmp, ", na.rm=TRUE)")))
      }
    }

  for (i in seq_len(nrow(A)))
    for (j in seq_len(ncol(A))) {
      if (grepl("data.", RAM$A[i, j], fixed=TRUE)) {
        tmp <- strsplit(RAM$A[i, j], "data.", fixed=TRUE)[[1]][2]
        A[i, j] <- eval(parse(text=paste0("mean(data$", tmp, ", na.rm=TRUE)")))
      }
    }

  for (j in seq_len(ncol(M))) {
    if (grepl("data.", RAM$M[1, j], fixed=TRUE)) {
      tmp <- strsplit(RAM$M[1, j], "data.", fixed=TRUE)[[1]][2]
      M[1, j] <- eval(parse(text=paste0("mean(data$", tmp, ", na.rm=TRUE)")))
    }
  }

  ## index of observed variables
  allNames <- colnames(A)
  if (length(shadow.vars) > 0) {
    ## F itself (not just manNames/latNames below) determines how many
    ## "manifest" nodes semPlot builds -- the fitted between-level F has
    ## zero rows, so it must be rebuilt to actually select shadow.vars, not
    ## just relabelled.
    F <- create.Fmatrix(allNames %in% shadow.vars, as.mxMatrix=FALSE)
    dimnames(F) <- list(shadow.vars, allNames)
    manNames <- shadow.vars
    latNames <- setdiff(allNames, shadow.vars)
  } else {
    index_obs <- (apply(F, 2, sum)==1)
    manNames <- allNames[index_obs]
    latNames <- allNames[!index_obs]
  }

  sem.plot <- semPlot::ramModel(A=A, S=S, F=F, M=M, manNames=manNames, latNames=latNames)

  out <- semPlot::semPaths(sem.plot, what=what, nCharNodes=nCharNodes,
                           nCharEdges=nCharEdges, layout=layout,
                           sizeMan=sizeMan, sizeLat=sizeLat,
                           edge.label.cex=edge.label.cex, color=color,
                           weighted=weighted, ...)
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

## Run "mx.model" via mxRun(); if that errors, retry via mxTryHard() and a
## final mxRun(). Shared by sem() (single-group), .sem.multigroup(), and
## .sem.twolevel(), which used to each inline their own copy of this logic.
##
## Fixes a bug present in all three former copies: they wrapped the
## mxTryHard()/mxRun() retry in tryCatch(..., warning=function(w) {...}),
## and that handler's return value (from cat(), used only for its
## printing side effect) was implicitly NULL -- not an object classed
## "error". Two consequences followed. (1) tryCatch()'s warning handler
## unwinds and returns on the FIRST warning raised anywhere during the
## wrapped call, so a single warning from deep inside mxTryHard()
## immediately discarded the whole in-progress retry and replaced it with
## NULL -- defeating the very purpose of mxTryHard(), which exists to
## push through exactly this kind of warning across many attempts. (2)
## When the model genuinely could not be fit at all (e.g. a missing data
## column), the final mx.fit ended up as plain NULL, which is not classed
## "error"; the subsequent `inherits(mx.fit, "error")` check was therefore
## always FALSE, so the intended `warning("Error in running mxModel.\n")`
## never fired and the caller had no indication anything had failed until
## summary()/plot() crashed later with unrelated-looking internal errors
## (e.g. "$ operator is invalid for atomic vectors").
##
## Fixed here by using withCallingHandlers() for warnings -- which records
## the message but resumes execution via muffleWarning instead of
## unwinding, so mxTryHard() can actually finish its job -- and by having
## the error handler return the error condition itself (keeping it classed
## "error") instead of the implicit NULL from cat().
.run.mxTryHard.fallback <- function(mx.model, intervals.type, ...) {
  warning.msg <- error.msg <- NULL

  mx.fit <- tryCatch(mxRun(mx.model, intervals=(intervals.type=="LB"),
                           suppressWarnings=TRUE, silent=TRUE, ...),
                     error=function(e) e)

  if (inherits(mx.fit, "error")) {
    mx.fit <- withCallingHandlers(
      tryCatch(mxTryHard(mx.model, extraTries=50, intervals=FALSE, silent=TRUE),
               error=function(e) {
                 error.msg <<- conditionMessage(e)
                 cat(error.msg)
                 e
               }),
      warning=function(w) {
        warning.msg <<- conditionMessage(w)
        cat(warning.msg)
        invokeRestart("muffleWarning")
      })

    if (!inherits(mx.fit, "error")) {
      mx.fit <- withCallingHandlers(
        tryCatch(mxRun(mx.fit, intervals=(intervals.type=="LB"),
                       suppressWarnings=TRUE, silent=TRUE, ...),
                 error=function(e) {
                   error.msg <<- conditionMessage(e)
                   cat(error.msg)
                   e
                 }),
        warning=function(w) {
          warning.msg <<- conditionMessage(w)
          cat(warning.msg)
          invokeRestart("muffleWarning")
        })
    }

    if (inherits(mx.fit, "error")) {
      warning("Error in running mxModel.\n")
      mx.fit <- NULL
    }
  }

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
             else if (!is.null(object$warning)) object$warning
             else "unknown reason"
    stop("The model was not fitted successfully, so there is no fitted ",
        "model to use here. Reported reason: ", reason, call.=FALSE)
  }
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

  checkRAM(Amatrix, Smatrix, cor.analysis=FALSE)

  var.names <- colnames(Fmatrix$values)

  mx.data <- mxData(observed=data, type="raw")
  expFun <- mxExpectationRAM(A="Amatrix", S="Smatrix", F="Fmatrix",
                             M="Mmatrix", dimnames=var.names)
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
  A0 <- A; S0 <- S; M0 <- M
  if (!is.null(replacements)) {
    for (lab in names(replacements)) {
      rhs <- replacements[[lab]]
      if (any(lab==A)) A[which(lab==A)] <- rhs
      if (any(lab==S)) S[which(lab==S)] <- rhs
      if (any(lab==M)) M[which(lab==M)] <- rhs
    }
  }

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
                            intervals.type=c("z","LB"), startvalues=NULL,
                            lbound=NULL, ubound=NULL, replace.constraints=FALSE,
                            mxModel.Args=NULL, run=TRUE, silent=TRUE, ...) {

  intervals.type <- match.arg(intervals.type)

  if (is.null(data) || !is.data.frame(data)) {
    stop("Multiple-group models require a single data frame passed via ",
        "'data', plus 'group' naming the column that identifies group ",
        "membership. Cov/numObs-based input is not yet supported for ",
        "multiple-group models.")
  }
  if (is.null(group) || !(group %in% colnames(data))) {
    stop("'group' must name a column in 'data' identifying group ",
        "membership for multiple-group models.")
  }
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
  empty <- vapply(data.list, nrow, integer(1)) == 0
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
    index <- vapply(mxalgebras, function(x) {
      form <- as.character(x$formula)
      form[1]=="==" && form[2] %in% para.labels
    }, logical(1))

    if (any(index)) {
      mxalgebras.const <- mxalgebras[index]
      const.labels <- vapply(mxalgebras.const, function(x) as.character(x$formula)[2],
                             character(1))
      const.rhs <- vapply(mxalgebras.const, function(x) as.character(x$formula)[3],
                          character(1))
      const.rhs <- .resolve.chained.constraints(const.labels, const.rhs)
      replacements <- const.rhs
      names(replacements) <- const.labels
      mxalgebras[index] <- NULL

      ## Symbols newly introduced via the constraints' own RHS text (e.g.
      ## "a0"/"a1" in "exp(a0+a1*data.x)"), excluding definition-variable
      ## ("data.*") references, which are never free parameters.
      introduced.labels <- all.vars(parse(text=const.rhs))
      introduced.labels <- introduced.labels[!startsWith(introduced.labels, "data.")]

      ## Any REMAINING ":="/"=="/"<"/">" entry that still references one
      ## of the just-eliminated labels needs the same substitution applied
      ## to its own formula (see .substitute.remaining.mxalgebras()).
      subs <- setNames(lapply(as.list(const.rhs), function(r) parse(text=r)[[1]]),
                       const.labels)
      mxalgebras <- .substitute.remaining.mxalgebras(mxalgebras, subs)

      ## Starting-value collision guard, same fix as two-level: a
      ## constraint's RHS can reference an *untouched* ordinary free
      ## parameter that lives in a different group (or the same one).
      ## Collect every group's own existing "value*label" starting values
      ## up front and thread them through as a fallback, so OpenMx doesn't
      ## see two conflicting starting values for the same label (as.
      ## mxAlgebra()'s own Xvars matrix otherwise defaults every symbol it
      ## discovers to a flat starting value of 0, colliding with that
      ## parameter's real one -- confirmed by testing: without this,
      ## "resid2 == 2*resid1" errors outright with "the free parameter
      ## 'resid1' has been assigned multiple starting values").
      ram.starts <- function(x) {
        lab <- gsub("^[^*]*\\*", "", x)
        val <- suppressWarnings(as.numeric(gsub("\\*.*$", "", x)))
        ok <- !is.na(val) & lab!=x & !startsWith(lab, "data.")
        setNames(as.list(val[ok]), lab[ok])
      }
      ## unname() on the outer lapply() result is required: RAM (a
      ## RAM_multigroup list) is itself named "1","2",... by lavaan2RAM(),
      ## and do.call(c, <named list of named lists>) triggers R's
      ## automatic outer.inner name-concatenation -- silently producing
      ## names like "1.sigma1" instead of "sigma1", which never matches
      ## anything in as.mxAlgebra()'s own Xvars$labels lookup, leaving
      ## this guard completely inert. Confirmed by testing.
      existing.starts <- do.call(c, unname(lapply(RAM, function(r) {
        c(ram.starts(r$A), ram.starts(r$S), ram.starts(r$M))
      })))
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
               "Id", "Id_A", "expCov", "expMean", para.labels, introduced.labels)
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
      if (!grepl("^constraint[0-9]", name.mxalgebra)) {
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
    fallback <- .run.mxTryHard.fallback(mx.model, intervals.type, ...)
    mx.fit <- fallback$mx.fit
    warning.msg <- fallback$warning.msg
    error.msg <- fallback$error.msg
  } else {
    mx.fit <- mx.model
  }

  out <- list(mx.fit=mx.fit, RAM=RAM, data=data, mxalgebras=mxalgebras.ci,
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

## Translate a RAM-style character A/S/M matrix trio into a list of
## mxPath() objects (the "local" structure of one level -- no cross-level
## joinKey paths, which are added separately by the caller).
##
## This exists because OpenMx's type="RAM" model construction only reliably
## wires up an mxExpectationRAM() (in particular its mean vector, and any
## cross-model joinKey linkage) when the ENTIRE local structure is expressed
## via mxPath() -- see dev/PLAN-two-level-multigroup-sem.md section 4.1.
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
## "=="/"<"/">" row by the same "^constraint[0-9]+$" naming convention
## already used elsewhere in this file for exactly this purpose.
.substitute.remaining.mxalgebras <- function(mxalgebras, subs) {
  if (is.null(mxalgebras) || length(subs) == 0) return(mxalgebras)

  for (nm in names(mxalgebras)) {
    formula <- mxalgebras[[nm]]$formula
    if (length(intersect(all.vars(formula), names(subs))) == 0) next

    new.text <- deparse1(.substitute.symbols(formula, subs))
    if (grepl("^constraint[0-9]+$", nm)) {
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
## variable). See dev/PLAN-two-level-multigroup-sem.md section 4 and
## dev/debug-two-level-multigroup/03-04 for the validated design.
.sem.twolevel <- function(model.name="sem", RAM, data=NULL, cluster=NULL,
                          intervals.type=c("z","LB"), startvalues=NULL,
                          lbound=NULL, ubound=NULL, replace.constraints=FALSE,
                          mxModel.Args=NULL, run=TRUE, silent=TRUE, ...) {

  intervals.type <- match.arg(intervals.type)

  if (is.null(data) || !is.data.frame(data)) {
    stop("Two-level models require a single data frame of individual-level ",
        "(level-1) rows passed via 'data', plus 'cluster' naming the column ",
        "that identifies level-2 (cluster) membership. Cov/numObs-based ",
        "input is not yet supported for two-level models.")
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

  ## ---- Level 2 ("between"): one row per cluster ----
  between.build <- .RAM2paths(RAM$between$A, RAM$between$S, RAM$between$M,
                              startvalues=startvalues, lbound=lbound, ubound=ubound)
  between.data <- data.frame(unique(data[[cluster]]))
  colnames(between.data) <- cluster
  between.model <- mxModel("between", type="RAM",
                           latentVars=colnames(RAM$between$A),
                           mxData(between.data, type="raw", primaryKey=cluster),
                           between.build$paths)

  ## ---- Level 1 ("within"): individual rows, nests "between" ----
  within.build <- .RAM2paths(RAM$within$A, RAM$within$S, RAM$within$M,
                             startvalues=startvalues, lbound=lbound, ubound=ubound)
  manifest.vars <- rownames(RAM$within$F)
  latent.vars <- setdiff(colnames(RAM$within$A), manifest.vars)
  shadow.vars <- intersect(manifest.vars, colnames(RAM$between$A))

  join.paths <- NULL
  if (length(shadow.vars) > 0) {
    join.paths <- mxPath(from=paste0("between.", shadow.vars), to=shadow.vars,
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
  Aw <- Aw0 <- as.symMatrix(RAM$within$A)
  Sw <- Sw0 <- as.symMatrix(RAM$within$S)
  Mw <- Mw0 <- as.symMatrix(RAM$within$M)
  Ab <- Ab0 <- as.symMatrix(RAM$between$A)
  Sb <- Sb0 <- as.symMatrix(RAM$between$S)
  Mb <- Mb0 <- as.symMatrix(RAM$between$M)

  if (isTRUE(replace.constraints) && !is.null(mxalgebras)) {
    para.labels <- unique(c(within.build$para.labels, between.build$para.labels))
    index <- vapply(mxalgebras, function(x) {
      form <- as.character(x$formula)
      form[1]=="==" && form[2] %in% para.labels
    }, logical(1))

    if (any(index)) {
      mxalgebras.const <- mxalgebras[index]

      ## Resolve any constraint whose RHS references another constraint's
      ## OWN label to a fixed point BEFORE applying any of them -- see
      ## .resolve.chained.constraints() for why (confirmed by testing:
      ## without this, "residW == exp(a0); residB == 2*residW" silently
      ## reintroduced "residW" as a disconnected free parameter instead of
      ## leaving the model expressed purely in terms of "a0").
      const.labels <- vapply(mxalgebras.const, function(x) as.character(x$formula)[2],
                             character(1))
      const.rhs <- vapply(mxalgebras.const, function(x) as.character(x$formula)[3],
                          character(1))
      const.rhs <- .resolve.chained.constraints(const.labels, const.rhs)

      for (i in seq_along(const.labels)) {
        lab <- const.labels[i]
        rhs <- const.rhs[i]
        if (any(lab==Aw)) Aw[which(lab==Aw)] <- rhs
        if (any(lab==Sw)) Sw[which(lab==Sw)] <- rhs
        if (any(lab==Mw)) Mw[which(lab==Mw)] <- rhs
        if (any(lab==Ab)) Ab[which(lab==Ab)] <- rhs
        if (any(lab==Sb)) Sb[which(lab==Sb)] <- rhs
        if (any(lab==Mb)) Mb[which(lab==Mb)] <- rhs
      }
      mxalgebras[index] <- NULL

      ## Any REMAINING ":="/"=="/"<"/">" entry that still references one
      ## of the just-eliminated labels needs the same substitution applied
      ## to its own formula (see .substitute.remaining.mxalgebras()).
      subs <- setNames(lapply(as.list(const.rhs), function(r) parse(text=r)[[1]]),
                       const.labels)
      mxalgebras <- .substitute.remaining.mxalgebras(mxalgebras, subs)

      ## A constraint's RHS can reference a label that is itself an
      ## ORDINARY, untouched free parameter elsewhere in the model (most
      ## often across levels, e.g. "residB == 2*residW" leaves "residW"
      ## alone) -- as.mxAlgebra() (inside .twolevel.replace.one() below)
      ## otherwise defaults every symbol it discovers to a flat starting
      ## value of 0, colliding with that parameter's real starting value
      ## (e.g. "0.5*residW" in the untouched matrix) and raising OpenMx's
      ## "has been assigned multiple starting values" error on the very
      ## first attempt -- self-corrected by the mxTryHard() retry that
      ## follows, but needlessly so (confirmed by testing). Collecting
      ## every already-labelled parameter's own starting value up front and
      ## passing it through to every .twolevel.replace.one() call below
      ## avoids the conflict outright. Any label the caller supplied
      ## through "startvalues" still wins (listed second: as.mxAlgebra()
      ## applies its "startvalues" list in order, so a later entry for the
      ## same label overwrites an earlier one).
      ram.starts <- function(x) {
        lab <- gsub("^[^*]*\\*", "", x)
        val <- suppressWarnings(as.numeric(gsub("\\*.*$", "", x)))
        ok <- !is.na(val) & lab!=x & !startsWith(lab, "data.")
        setNames(as.list(val[ok]), lab[ok])
      }
      existing.starts <- c(ram.starts(RAM$within$A), ram.starts(RAM$within$S),
                           ram.starts(RAM$within$M), ram.starts(RAM$between$A),
                           ram.starts(RAM$between$S), ram.starts(RAM$between$M))
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
  }

  ## Genuine level-2-only (between, cluster-level) definition variables are
  ## not yet supported: "between.data" (built just above) only ever carries
  ## the cluster ID column, never a real cluster-level covariate, so a
  ## "data.*" reference among the (possibly replace.constraints-substituted)
  ## between-level matrices can never resolve -- whether written directly
  ## in the model (e.g. "y ~~ data.z*y" at "level: 2") or introduced via a
  ## constraint substitution (e.g. "residB == exp(c0+c1*data.z)"). Left
  ## unchecked, this reaches mxRun() as a hard R error that sem()'s own
  ## mxRun-then-mxTryHard() fallback mistakes for "maybe just needs a
  ## better start", burning through 50 retries before finally re-raising
  ## OpenMx's own message -- which does not even mention that between-level
  ## definition variables are the actual limitation. Caught here instead,
  ## before any of that, with a message that says so directly. (A
  ## within-level, i.e. "level: 1", definition variable is fully supported
  ## and is not affected by this check.)
  between.vars <- all.vars(parse(text=c(Ab, Sb, Mb)))
  between.defvars <- between.vars[startsWith(between.vars, "data.")]
  missing.cols <- setdiff(sub("^data\\.", "", between.defvars), colnames(between.data))
  if (length(missing.cols) > 0) {
    stop("Genuine level-2-only (between, cluster-level) definition ",
        "variable(s) are not yet supported: ",
        paste0("'data.", missing.cols, "'", collapse=", "),
        " reference column(s) that the level-2 (between) data does not ",
        "carry -- it currently contains only the cluster ID. A within-",
        "level (\"level: 1\") definition variable (e.g. \"data.vi\") is ",
        "supported; a genuine cluster-level covariate referenced at ",
        "\"level: 2\" (directly, or via a replace.constraints substitution) ",
        "is not.", call.=FALSE)
  }

  mx.model <- mxModel(model.name, type="RAM", between.model,
                      manifestVars=manifest.vars, latentVars=latent.vars,
                      mxData(data, type="raw"),
                      within.build$paths, join.paths)

  ## Within-level replacements land on "mx.model" itself (the top model
  ## IS the "within" level here), so this can only happen after it exists.
  if (isTRUE(replace.constraints) && !is.null(RAM$mxalgebras)) {
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
  if (isTRUE(replace.constraints) && !is.null(RAM$mxalgebras)) {
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
      if (!grepl("^constraint[0-9]", name.mxalgebra)) {
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
    fallback <- .run.mxTryHard.fallback(mx.model, intervals.type, ...)
    mx.fit <- fallback$mx.fit
    warning.msg <- fallback$warning.msg
    error.msg <- fallback$error.msg
  } else {
    mx.fit <- mx.model
  }

  out <- list(mx.fit=mx.fit, RAM=RAM, data=data, mxalgebras=mxalgebras.ci,
              intervals.type=intervals.type, warning=warning.msg,
              error=error.msg)
  class(out) <- "mxsem"
  out
}
