#' Estimate Variance Components with Restricted (Residual) Maximum Likelihood
#' Estimation
#'
#' It estimates the variance components of random-effects in univariate and
#' multivariate meta-analysis with restricted (residual) maximum likelihood
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
#' An alternative but equivalent approach is to minimize the -2*log-likelihood
#' function: \deqn{ }{
#' log(det|V+T^2|)+log(det|X'(V+T^2)^{-1}X|)+(y-X\alpha)'(V+T^2)^{-1}(y-X*\alpha)}\deqn{
#' \log(\det|V+T^2|)+\log(\det|X'(V+T^2)^{-1}X|)+(y-X\hat{\alpha})'(V+T^2)^{-1}(y-X\hat{\alpha})}{
#' log(det|V+T^2|)+log(det|X'(V+T^2)^{-1}X|)+(y-X\alpha)'(V+T^2)^{-1}(y-X*\alpha)}
#' where \eqn{V}{V} is the known conditional sampling covariance matrix of
#' \eqn{y}{y}, \eqn{T^2}{T^2} is the variance component of the random effects,
#' and \eqn{\hat{\alpha}=(X'(V+T^2)^{-1}X)^{-1}
#' }{\hat{\alpha}=(t(X)(V+T^2)^{-1}X)^{-1}t(X)(V+T^2)^{-1}y}\eqn{
#' X'(V+T^2)^{-1}y}{\hat{\alpha}=(t(X)(V+T^2)^{-1}X)^{-1}t(X)(V+T^2)^{-1}y}.
#' \code{reml()} minimizes the above likelihood function to obtain the
#' parameter estimates.
#'
#' @param y A vector of effect size for univariate meta-analysis or a
#' \eqn{k}{k} x \eqn{p}{p} matrix of effect sizes for multivariate
#' meta-analysis where \eqn{k}{k} is the number of studies and \eqn{p}{p} is
#' the number of effect sizes.
#' @param v A vector of the sampling variance of the effect size for univariate
#' meta-analysis or a \eqn{k}{k} x \eqn{p*}{p*} matrix of the sampling
#' covariance matrix of the effect sizes for multivariate meta-analysis where
#' \eqn{p* = p(p+1)/2 }{p* = p(p+1)/2}. It is arranged by column major as used
#' by \code{\link[OpenMx]{vech}}.
#' @param x A predictor or a \eqn{k}{k} x \eqn{m}{m} matrix of predictors where
#' \eqn{m}{m} is the number of predictors.
#' @param data An optional data frame containing the variables in the model.
#' @param RE.constraints A \eqn{p}{p} x \eqn{p}{p} matrix specifying the
#' variance components of the random effects. If the input is not a matrix, it
#' is converted into a matrix by \code{as.matrix()}. The default is that all
#' covariance/variance components are free. The format of this matrix follows
#' \code{\link[metaSEM]{as.mxMatrix}}. Elements of the variance components can
#' be constrained equally by using the same labels. If a zero matrix is
#' specified, it becomes a fixed-effects meta-analysis.
#' @param RE.startvalues A vector of \eqn{p}{p} starting values on the
#' diagonals of the variance component of the random effects. If only one
#' scalar is given, it will be repeated across the diagonals. Starting values
#' for the off-diagonals of the variance component are all 0. A \eqn{p}{p} x
#' \eqn{p}{p} symmetric matrix of starting values is also accepted.
#' @param RE.lbound A vector of \eqn{p}{p} lower bounds on the diagonals of the
#' variance component of the random effects. If only one scalar is given, it
#' will be repeated across the diagonals. Lower bounds for the off-diagonals of
#' the variance component are set at \code{NA}. A \eqn{p}{p} x \eqn{p}{p}
#' symmetric matrix of the lower bounds is also accepted.
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
#' @param silent Logical. An argument to be passed to
#' \code{\link[OpenMx]{mxRun}}
#' @param run Logical. If \code{FALSE}, only return the mx model without
#' running the analysis.
#' @param \dots Further arguments to be passed to \code{\link[OpenMx]{mxRun}}
#' @return An object of class \code{reml} with a list of \item{call}{Object
#' returned by \code{\link[base]{match.call}}} \item{data}{A data matrix of y,
#' v and x } \item{no.y}{No. of effect sizes} \item{no.x}{No. of predictors}
#' \item{miss.vec}{A vector indicating missing data. Studies will be removed
#' before the analysis if they are \code{TRUE}} \item{mx.fit}{A fitted object
#' returned from \code{\link[OpenMx]{mxRun}}}
#' @note \code{reml} is more computationally intensive than \code{meta}.
#' Moreover, \code{reml} is more likely to encounter errors during
#' optimization. Since a likelihood function is directly employed to obtain the
#' parameter estimates, there is no number of studies and number of observed
#' statistics returned by \code{\link[OpenMx]{mxRun}}. Ad-hoc steps are used to
#' modify \code{mx.fit@runstate$objectives[[1]]@numObs} and
#' \code{mx.fit@runstate$objectives[[1]]@numStats}.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{meta}}, \code{\link[metaSEM]{reml3}},
#' \code{\link[metaSEM]{Hox02}}, \code{\link[metaSEM]{Berkey98}}
#' to See Also as \code{\link{help}}, ~~~
#' @references Cheung, M. W.-L. (2013). Implementing restricted maximum
#' likelihood estimation in structural equation models. \emph{Structural
#' Equation Modeling}, \bold{20(1)}, 157-167.
#'
#' Mehta, P. D., & Neale, M. C. (2005). People Are Variables Too: Multilevel
#' Structural Equations Modeling. \emph{Psychological Methods}, \bold{10(3)},
#' 259-284.
#'
#' Searle, S. R., Casella, G., & McCulloch, C. E. (1992). \emph{Variance
#' components}. New York: Wiley.
#'
#' Viechtbauer, W. (2005). Bias and efficiency of meta-analytic variance
#' estimators in the random-effects model. \emph{Journal of Educational and
#' Behavioral Statistics}, \bold{30(3)}, 261-293.
#' @keywords meta-analysis
reml <- function(y, v, x, data, RE.constraints=NULL, RE.startvalues=0.1, RE.lbound=1e-10,
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
  y <- eval(my.y, data, enclos = sys.frame(sys.parent()))
  v <- eval(my.v, data, enclos = sys.frame(sys.parent()))
  
  if (is.vector(y)) {
    no.y <- 1
    no.studies <- length(y)
  } else {
    no.y <- ncol(y)
    no.studies <- nrow(y)
  }
  if (is.vector(v)) {
    no.v <- 1
  } else {
    no.v <- ncol(v)
  }
  if (missing(x)) {
    no.x <- 0
  } else {
    my.x <- mf[[match("x", names(mf))]]
    x <- eval(my.x, data, enclos = sys.frame(sys.parent()))
    if (is.vector(x)) {
      no.x <- 1
    } else {
      no.x <- ncol(x)
    }
  }
   
  # Check the dimensions of y and v
  if ( no.v != no.y*(no.y+1)/2 )
    stop(paste("The expected no. of columns in v is ", no.y*(no.y+1)/2,
               " while the observed no. of columns in v is ", no.v, ".", sep=""))
  v.labels <- vech(outer(1:no.y, 1:no.y, function(x, y) paste("v", x,"_", y, sep = "")))
  y.labels <- paste("y", 1:no.y, sep="")
  x.labels <- paste("x", 1:no.x, sep="")

  ## Create a design matrix: each row represents 1 study
  ## It is not used in the actural analysis.
  ## p: no. of predictors + intercept; remember to multiply by no.y
  if (missing(x)) {
    x.design <- matrix(1, ncol=1, nrow=no.studies)
    p <- no.y
    no.x <- 0
    input.df <- as.matrix(cbind(y,v))
    dimnames(input.df) <- list(NULL, c(y.labels, v.labels))
  } else {
    x.design <- cbind(1, x)
    no.x <- ncol(x.design)
    p <- no.x*no.y
    input.df <- as.matrix(cbind(y, v, x))
    dimnames(input.df) <- list(NULL, c(y.labels, v.labels, x.labels)) 
  }
  
  ## Create Y, a row vector of effect sizes
  ## nrow(Y)= no.y*no.studies
  Y <- c(t(y))
  ## A vector indicating missing values
  miss.vec <- is.na(Y)
  # miss.vec <- matrix( t(is.na(Y)), ncol=1 ) 
  # Remove missing values
  Y <- Y[!miss.vec]
  # Convert it into a column vector
  Y <- matrix(Y, ncol=1)
  # No. of total effect sizes after removing missing data
  no.es <- length(Y)
  # Ad-hoc: no. of observed statistics after removing the fixed-effects (p)
  numStats <- no.es-p
  Y <- as.mxMatrix(Y)
  
  # Function to create design matrix for a study, e.g., x_1=1,2,3, no.y=2
  # 1 0 2 0 3 0
  # 0 1 0 2 0 3
  fn1 <- function(x, no.y) {
    temp <- lapply(x, function(x, k){Diag(x=x, nrow=k, ncol=k)}, k=no.y)
    do.call(cbind, temp)
  }
  # temp: a list of design matrix per study
  temp <- lapply(split(x.design, 1:nrow(x.design)), fn1, no.y=no.y)
  # Convert the list into a design matrix
  # X: based on no. of effect sizes
  X <- do.call(rbind, temp)
  # A design matrix based on Y as a column vector
  X <- X[!miss.vec, , drop=FALSE ]
  X <- as.mxMatrix(X)
  
  # McCulloch (2003) Generalized linear mixed models. E-book p. 67
  # Searle, Casella, & McCulloch (1992). Variance Components. E-book p. 250
  # I - X (X'X)^-1 X'
  # crossprod(X)=X'X
  ## M <- diag(no.es) - X %&% solve( crossprod(X) )
  ## ##N <- diag(no.es) - X %*% solve( t(X)%*% X ) %*% t(X)
  ## M <- M[1:(no.es-p), ]  # remove redundant columns

  ## ## transformed effect size: A row vector as required
  ## y_star <- t( M %*% Y )  
  ## selVars <- paste("X", 1:(no.es-p), sep="")
  ## dimnames(y_star) <- list(NULL, selVars)
  ## M <- as.mxMatrix(M, name="M")
  
  ## matrix of lbound
  if (is.matrix(RE.lbound)) {
    
    if (!all(dim(RE.lbound)==c(no.y, no.y)))
      stop("Dimensions of \"RE.lbound\" are incorrect.")    
    } else {

      lbound <- matrix(NA, nrow=no.y, ncol=no.y)
      Diag(lbound) <- RE.lbound
    }  
  ## Convert it into a large matrix
  lbound <- bdiagRep(lbound, no.studies)
  lbound <- lbound[!miss.vec, !miss.vec]

  free <- bdiagRep( matrix(TRUE, nrow=no.y, ncol=no.y), no.studies )
  free <- free[!miss.vec, !miss.vec]
  ## free is a vector of logical values
  free <- as.logical(vech(free))
  
  ## Preparing the S matrix for covariance elements
  ## No predictor
  if (is.null(RE.constraints)) {
    
    if (is.matrix(RE.startvalues)) {
      # FIXME: test symmetry
      if (!all(dim(RE.startvalues)==c(no.y, no.y)))
        # stop instead of warning here
        stop("Dimensions of \"RE.startvalues\" are incorrect.")                    
      #values <- matrix(c(RE.startvalues), nrow=no.y, ncol=no.y)
    } else {      
      values <- Diag(x=RE.startvalues, nrow=no.y, ncol=no.y)
    }    
    # Large matrix
    values <- bdiagRep(values, no.studies)
    values <- values[!miss.vec, !miss.vec]
    
    Tau.labels <- vech(outer(1:no.y, 1:no.y, function(x,y) { paste("Tau2_",x,"_",y,sep="")}))
    # Large matrix
    Tau.labels <- bdiagRep(vec2symMat(Tau.labels), no.studies)
    Tau.labels <- Tau.labels[!miss.vec, !miss.vec]
    # Replace off diagonals with NA
    Tau.labels[Tau.labels=="0"] <- NA
    
    Tau <- mxMatrix("Symm", ncol=no.es, nrow=no.es, free=free, labels=vech(Tau.labels),
                    lbound=vech(lbound), values=vech(values), name="Tau")
  } else {
    ## Convert RE.constraints into a column matrix if it is not a matrix
    if (!is.matrix(RE.constraints))
      RE.constraints <- as.matrix(RE.constraints)
    
    if (!all(dim(RE.constraints)==c(no.y, no.y)))
      stop("Dimensions of \"RE.constraints\" are incorrect.")

    Tau <- bdiagRep(RE.constraints, no.studies)
    Tau <- Tau[!miss.vec, !miss.vec]
    Tau <- as.mxMatrix(Tau, lbound=lbound, name="Tau")
  }

  ## Known sampling variance matrix
  if (no.y==1) {
    V <- Diag(x=v, nrow=no.studies, ncol=no.studies)
  } else {
    V <- matrix2bdiag(v)
  }
  V <- V[!miss.vec, !miss.vec]
  V <- as.mxMatrix(V)

  # Inverse of (V+Tau)
  W <- mxAlgebra(solve(V+Tau), name="W")
  alpha <- mxAlgebra( solve(t(X)%*%W%*%X) %*% t(X) %*% W %*% Y, name="alpha")
  # -2LL
  obj <- mxAlgebra( ( log(det(V+Tau)) + log(det(t(X)%*%W%*%X)) +
                       t(Y-X%*%alpha)%*%W%*%(Y-X%*%alpha) ), name="obj")

  # Creat model for REML
  reml.model <- mxModel(model=model.name, X, Y, V, W, Tau, alpha, obj,
                   #mxData(observed=y_star, type="raw"),
                   mxFitFunctionAlgebra(algebra="obj", numObs=no.studies, numStats=numStats), mxCI("Tau"))

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

  ## Ad-hoc: Add no. of studies and no. of observed statistics
  ## mx.fit@runstate$objectives[[1]]@numObs <- no.studies
  ## mx.fit@runstate$objectives[[1]]@numStats <- numStats
  
  out <- list(call = mf, data=input.df, no.y=no.y, no.x=no.x, miss.vec=miss.vec, mx.model=reml.model,
              mx.fit=mx.fit, intervals.type=intervals.type, numObs=no.studies, numStats=numStats)
  class(out) <- "reml"
  return(out)
}
