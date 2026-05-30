#' Estimate Variance Components in Three-Level Univariate Meta-Analysis with
#' Restricted (Residual) Maximum Likelihood Estimation
#'
#' It estimates the variance components of random-effects in three-level
#' univariate meta-analysis with restricted (residual) maximum likelihood
#' (REML) estimation method.
#'
#' Restricted (residual) maximum likelihood obtains the parameter estimates on
#' the transformed data that do not include the fixed-effects parameters. A
#' transformation matrix \eqn{M=I-X(X'X)^{-1}X}{M=I-X(X'X)^{-1}X'} is created
#' based on the design matrix \eqn{X}{X} which is just a column vector when
#' there is no predictor in \code{x}. The last \eqn{N}{N} redundant rows of
#' \eqn{M}{M} is removed where \eqn{N}{N} is the rank of \eqn{X}{X}. After
#' pre-multiplying by \eqn{M} on \code{y}, the parameters of fixed-effects are
#' removed from the model. Thus, only the parameters of random-effects are
#' estimated.
#'
#' An alternative but equivalent approach is to minimize the
#' -2*log-likelihood function: \deqn{ }{
#' log(det|V+T^2|)+log(det|X'(V+T^2)^{-1}X|)+(y-X\alpha)'(V+T^2)^{-1}(y-X*\alpha)}\deqn{
#' \log(\det|V+T^2|)+\log(\det|X'(V+T^2)^{-1}X|)+(y-X\hat{\alpha})'(V+T^2)^{-1}(y-X\hat{\alpha})}{
#' log(det|V+T^2|)+log(det|X'(V+T^2)^{-1}X|)+(y-X\alpha)'(V+T^2)^{-1}(y-X*\alpha)}
#' where \eqn{V}{V} is the known conditional sampling covariance matrix of
#' \eqn{y}{y}, \eqn{T^2}{T^2} is the variance component combining level-2 and
#' level-3 random effects, and \eqn{\hat{\alpha}=(X'(V+T^2)^{-1}X)^{-1}
#' X'(V+T^2)^{-1}y}{\hat{\alpha}=(t(X)(V+T^2)^{-1}X)^{-1}t(X)(V+T^2)^{-1}y}.
#' \code{reml()} minimizes the above likelihood function to obtain the
#' parameter estimates.
#'
#' @aliases reml3 reml3L
#' @param y A vector of \eqn{k}{k} studies of effect size.
#' @param v A vector of \eqn{k}{k} studies of sampling variance.
#' @param cluster A vector of \eqn{k}{k} characters or numbers indicating the
#' clusters.
#' @param x A predictor or a \eqn{k}{k} x \eqn{m}{m} matrix of level-2 and
#' level-3 predictors where \eqn{m}{m} is the number of predictors.
#' @param data An optional data frame containing the variables in the model.
#' @param RE2.startvalue Starting value for the level-2 variance.
#' @param RE2.lbound Lower bound for the level-2 variance.
#' @param RE3.startvalue Starting value for the level-3 variance.
#' @param RE3.lbound Lower bound for the level-3 variance.
#' @param RE.equal Logical. Whether the variance components at level-2 and
#' level-3 are constrained equally.
#' @param intervals.type Either \code{z} (default if missing) or \code{LB}. If
#' it is \code{z}, it calculates the 95\% Wald confidence intervals (CIs) based
#' on the z statistic. If it is \code{LB}, it calculates the 95\%
#' likelihood-based CIs on the parameter estimates. Note that the z values and
#' their associated p values are based on the z statistic. They are not related
#' to the likelihood-based CIs.
#' @param model.name A string for the model name in
#' \code{\link[OpenMx]{mxModel}}.
#' @param suppressWarnings Logical. If \code{TRUE}, warnings are suppressed. It
#' is passed to \code{\link[OpenMx]{mxRun}}.
#' @param silent Logical. Argument to be passed to \code{\link[OpenMx]{mxRun}}
#' @param run Logical. If \code{FALSE}, only return the mx model without
#' running the analysis.
#' @param \dots Further arguments to be passed to \code{\link[OpenMx]{mxRun}}
#' @return An object of class \code{reml} with a list of
#' \item{call}{Object returned by \code{\link[base]{match.call}}}
#' \item{data}{A data matrix of y, v, and x}
#' \item{mx.fit}{A fitted object returned from \code{\link[OpenMx]{mxRun}}}
#' @note \code{reml} is more computationally intensive than \code{meta}.
#' Moreover, \code{reml} is more likely to encounter errors during
#' optimization. Since a likelihood function is directly employed to obtain the
#' parameter estimates, there is no number of studies and number of observed
#' statistics returned by \code{\link[OpenMx]{mxRun}}. Ad hoc steps are used to
#' modify \code{mx.fit@runstate$objectives[[1]]@numObs} and
#' \code{mx.fit@runstate$objectives[[1]]@numStats}.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{meta3L}}, \code{\link[metaSEM]{reml}},
#' \code{\link[metaSEM]{Cooper03}}, \code{\link[metaSEM]{Bornmann07}}
#' @references Cheung, M. W.-L. (2013). Implementing restricted maximum
#' likelihood estimation in structural equation models. \emph{Structural
#' Equation Modeling}, \bold{20(1)}, 157-167.
#'
#' Cheung, M. W.-L. (2014). Modeling dependent effect sizes with three-level
#' meta-analyses: A structural equation modeling approach. \emph{Psychological
#' Methods}, \bold{19}, 211-229.
#'
#' Mehta, P. D., & Neale, M. C. (2005). People Are Variables Too: Multilevel
#' Structural Equations Modeling. \emph{Psychological Methods}, \bold{10(3)},
#' 259-284.
#'
#' Searle, S. R., Casella, G., & McCulloch, C. E. (1992). \emph{Variance
#' components}. New York: Wiley.
#' @keywords meta-analysis
reml3L <- function(y, v, cluster, x, data, RE2.startvalue=0.1, RE2.lbound=1e-10,
                  RE3.startvalue=RE2.startvalue, RE3.lbound=RE2.lbound, RE.equal=FALSE,
                  intervals.type=c("z", "LB"), model.name="Variance component with REML",
                  suppressWarnings=TRUE, silent=TRUE, run=TRUE, ...) {
  mf <- match.call()
  if (missing(data)) {
    data <- sys.frame(sys.parent())
  } else {
    if (!is.data.frame(data)) {
      data <- data.frame(data)
    }
  }
  my.y <- mf[[match("y", names(mf))]]
  my.v <- mf[[match("v", names(mf))]]
  my.cluster <- mf[[match("cluster", names(mf))]]  
  cluster <- as.character(eval(my.cluster, data, enclos = sys.frame(sys.parent())))
  ## check if there are missing data in cluster
  if (any(is.na(cluster)))
      stop("Missing values are not allowed in \"cluster\".\n")  
  cluster.order <- order(cluster)
  cluster <- cluster[cluster.order]
  y <- eval(my.y, data, enclos = sys.frame(sys.parent()))
  y <- y[cluster.order]
  v <- eval(my.v, data, enclos = sys.frame(sys.parent()))
  v <- v[cluster.order]
  no.studies <- length(y)

  ## Response matrix
  Y <- as.mxMatrix(matrix(y, ncol=1), name="Y")
  ## Intercept
  X <- matrix(1, ncol=1, nrow=no.studies)
  
  if (missing(x)) {
    no.x <- 0
  } else {
    my.x <- mf[[match("x", names(mf))]]
    x <- eval(my.x, data, enclos = sys.frame(sys.parent()))
    if (is.vector(x)) x <- x[cluster.order] else x <- x[cluster.order, ]
    X <- cbind(X, x)
    no.x <- ncol(X)-1
  }
  
  numStats <- no.studies-ncol(X)
  X <- as.mxMatrix(X)
  
  ## maximum no. of data in level-2 unit
  k <- lapply( split(cluster, cluster), length )

  ## Prepare V: fixed and known
  V <- as.mxMatrix( Diag(v), name="V" )

  ## Prepare Tau2 
  Tau_2 <- Diag(rep(paste(RE2.startvalue,"*Tau2_2", sep=""), no.studies))
  re2.lbound <- Diag(RE2.lbound, nrow=no.studies, ncol=no.studies)
  re2.lbound[re2.lbound==0] <- NA
  Tau_2 <- as.mxMatrix(Tau_2, name="Tau_2", lbound=re2.lbound)
  
  ## Prepare Tau3 
  Tau_3 <- lapply(k, function(x) {matrix(paste(RE3.startvalue,"*Tau2_3", sep=""), ncol=x, nrow=x)})
  Tau_3 <- bdiagMat(Tau_3)
  re3.lbound <- lapply(k, function(x) {matrix(RE3.lbound, ncol=x, nrow=x)})
  re3.lbound <- bdiagMat(re3.lbound)
  re3.lbound[re3.lbound==0] <- NA
  Tau_3 <- as.mxMatrix(Tau_3, name="Tau_3", lbound=re3.lbound)
 
  # Inverse of (V+Tau)
  W <- mxAlgebra(solve(V+Tau_2+Tau_3), name="W")
  alpha <- mxAlgebra( solve(t(X)%*%W%*%X) %*% t(X) %*% W %*% Y, name="alpha")
  # -2LL
  obj <- mxAlgebra( ( log(det(V+Tau_2+Tau_3)) + log(det(t(X)%*%W%*%X)) +
                       t(Y-X%*%alpha)%*%W%*%(Y-X%*%alpha) ), name="obj")

  # Creat model for REML
  ## reml.model <- mxModel(model=model.name, X, Y, V, Tau_2, Tau_3, W, alpha, obj,
  ##                       mxAlgebraObjective("obj"), mxCI(c("Tau2_2","Tau2_3")))

  # Constrain equal variances
  if (RE.equal) {
    reml.model <- mxModel(model=model.name, X, Y, V, Tau_2, Tau_3, W, alpha, obj,
                          mxFitFunctionAlgebra("obj"), mxCI(c("Tau2")))
    reml.model <- omxSetParameters(reml.model, labels=c("Tau2_2", "Tau2_3"), newlabels=c("Tau2", "Tau2"),
                                   values=c(RE2.startvalue, RE2.startvalue), lbound=c(RE2.lbound, RE2.lbound))
  } else {
    reml.model <- mxModel(model=model.name, X, Y, V, Tau_2, Tau_3, W, alpha, obj,
                          mxFitFunctionAlgebra(algebra="obj", numObs=no.studies, numStats=numStats), mxCI(c("Tau2_2","Tau2_3")))
  }

  ## ## Assuming NA first
  ## reml0.fit <- NA  
  ## if (no.x==0) {

  ##   ## Calculate I2
  ##   ## Based on Higgins and Thompson (2002), Eq. 9
  ##   sum.w <- sum(1/v)
  ##   sum.w2 <- sum(1/v^2)
  ##   ## Typical V based on Q statistic
  ##   qV <- matrix((no.studies-1)*sum.w/(sum.w^2-sum.w2), nrow=1, ncol=1)
  ##   qV <- as.mxMatrix(qV)
  ##   ## Typical V based on harmonic mean  
  ##   hmV <- matrix(no.studies/sum.w, nrow=1, ncol=1)
  ##   hmV <- as.mxMatrix(hmV)
  ##   ## Typical V based on arithmatic mean
  ##   amV <- matrix(mean(v))
  ##   amV <- as.mxMatrix(amV)

  ##   I2q_2 <- mxAlgebra( Tau2[1,1]/(Tau2[1,1]+Tau3[1,1]+qV), name="I2q_2")
  ##   I2q_3 <- mxAlgebra( Tau3[1,1]/(Tau2[1,1]+Tau3[1,1]+qV), name="I2q_3")
  ##   I2hm_2 <- mxAlgebra( Tau2[1,1]/(Tau2[1,1]+Tau3[1,1]+hmV), name="I2hm_2")
  ##   I2hm_3 <- mxAlgebra( Tau3[1,1]/(Tau2[1,1]+Tau3[1,1]+hmV), name="I2hm_3")  
  ##   I2am_2 <- mxAlgebra( Tau2[1,1]/(Tau2[1,1]+Tau3[1,1]+amV), name="I2am_2")
  ##   I2am_3 <- mxAlgebra( Tau3[1,1]/(Tau2[1,1]+Tau3[1,1]+amV), name="I2am_3") 
  ##   ICC_2 <- mxAlgebra( Tau2[1,1]/(Tau2[1,1]+Tau3[1,1]), name="ICC_2")
  ##   ICC_3 <- mxAlgebra( Tau3[1,1]/(Tau2[1,1]+Tau3[1,1]), name="ICC_3")

  ##   I2 <- match.arg(I2, c("I2q", "I2hm", "I2am", "ICC"), several.ok=TRUE)
  ##   ci <- c(outer(I2, c("_2","_3"), paste, sep=""))

  ##   reml.model <- mxModel(model=model.name, X, Y, V, Tau2, Tau3, W, alpha, obj,
  ##                         qV, hmV, amV, I2q_2, I2q_3, I2hm_2, I2hm_3, I2am_2, I2am_3, ICC_2, ICC_3,
  ##                         mxAlgebraObjective("obj"), mxCI(c("Tau_2","Tau3", ci)))
    
  ##   ## meta3 <- mxModel(model=model.name, mxData(observed=my.wide[,-1], type="raw"), oneRow, Id, Ones,
  ##   ##                  inter, coeff, mydata, Tau2, Tau3, V, expMean, expCov,
  ##   ##                  qV, hmV, amV, I2q_2, I2q_3, I2hm_2, I2hm_3, I2am_2, I2am_3, ICC_2, ICC_3,
  ##   ##                  mxFIMLObjective("expCov","expMean", dimnames=paste("y_", 1:k, sep="")),
  ##   ##                  mxCI(c("inter","coeff","Tau2","Tau3", ci)))
  ## } else {
  ##   ## no.x > 0
  ##   reml.model <- mxModel(model=model.name, X, Y, V, Tau2, Tau3, W, alpha, obj,
  ##                         mxAlgebraObjective("obj"), mxCI(c("Tau_2","Tau3")))
  ##   ## meta3 <- mxModel(model=model.name, mxData(observed=my.wide[,-1], type="raw"), oneRow, Id, Ones,
  ##   ##                  inter, coeff, mydata, Tau2, Tau3, V, expMean, expCov,
  ##   ##                  mxFIMLObjective("expCov","expMean", dimnames=paste("y_", 1:k, sep="")),
  ##   ##                  mxCI(c("inter","coeff","Tau2","Tau3")))

  ##   ## Calculate R2
  ##   if (R2) reml0.fit <- tryCatch( reml3L(y=y, v=v, cluster=cluster, data=my.long, model.name="No predictor",
  ##                                  suppressWarnings=TRUE, silent=TRUE), error = function(e) e )    
  ## }

  ## Return mx model without running the analysis
  if (run==FALSE) return(reml.model)
  
  intervals.type <- match.arg(intervals.type)
  # Default is z
  switch(intervals.type,
    z = mx.fit <- tryCatch( mxRun(reml.model, intervals=FALSE,
                                    suppressWarnings = suppressWarnings, silent=silent, ...), error = function(e) e ),
    LB = mx.fit <- tryCatch( mxRun(reml.model, intervals=TRUE,
                                     suppressWarnings = suppressWarnings, silent=silent, ...), error = function(e) e ) )
 
  if (inherits(mx.fit, "error")) {
    cat("Error in running the mxModel:\n")
    warning(print(mx.fit))
    return(mx.fit)
  }

  ## ## Ad-hoc: Add no. of studies and no. of observed statistics
  ## mx.fit@runstate$objectives[[1]]@numObs <- no.studies
  ## ## Ad-hoc: no. of observed statistics after removing the fixed-effects (p)
  ## mx.fit@runstate$objectives[[1]]@numStats <- numStats
  
  out <- list(call = mf, data=data, mx.model=reml.model, mx.fit=mx.fit, intervals.type=intervals.type, numObs=no.studies, numStats=numStats)
  class(out) <- c("reml", "reml3L")
  return(out)
}

reml3 <- function(y, v, cluster, x, data, RE2.startvalue=0.1, RE2.lbound=1e-10,
                  RE3.startvalue=RE2.startvalue, RE3.lbound=RE2.lbound, RE.equal=FALSE,
                  intervals.type=c("z", "LB"), model.name="Variance component with REML",
                  suppressWarnings=TRUE, silent=TRUE, run=TRUE, ...) {

    .Deprecated("reml3L")
    reml3L(y=y, v=v, cluster=cluster, x=x, data=data, RE2.startvalue=RE2.startvalue,
           RE2.lbound=RE2.lbound, RE3.startvalue=RE3.startvalue, RE3.lbound=RE3.lbound,
           RE.equal=RE.equal, intervals.type=intervals.type, model.name=model.name,
           suppressWarnings=suppressWarnings, silent=silent, run=run, ...)
}
