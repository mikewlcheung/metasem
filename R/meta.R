#' Univariate and Multivariate Meta-Analysis with Maximum Likelihood Estimation
#' 
#' It conducts univariate and multivariate meta-analysis with maximum
#' likelihood estimation method. Mixed-effects meta-analysis can be conducted
#' by including study characteristics as predictors. Equality constraints on
#' intercepts, regression coefficients, and variance components can be easily
#' imposed by setting the same labels on the parameter estimates.
#' 
#' 
#' @aliases meta metaFIML
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
#' @param av An auxiliary variable or a \eqn{k}{k} x \eqn{m}{m} matrix of
#' auxiliary variables where \eqn{m}{m} is the number of auxiliary variables.
#' @param data An optional data frame containing the variables in the model.
#' @param intercept.constraints A \eqn{1}{1} x \eqn{p}{p} matrix specifying
#' whether the intercepts of the effect sizes are fixed or free. If the input
#' is not a matrix, the input is converted into a \eqn{1}{1} x \eqn{p}{p}
#' matrix with \code{t(as.matrix(intercept.constraints))}. The default is that
#' the intercepts are free. When there is no predictor, these intercepts are
#' the same as the pooled effect sizes. The format of this matrix follows
#' \code{\link[metaSEM]{as.mxMatrix}}. The intercepts can be constrained
#' equally by using the same labels.
#' @param coef.constraints A \eqn{p}{p} x \eqn{m}{m} matrix specifying how the
#' predictors predict the effect sizes. If the input is not a matrix, it is
#' converted into a matrix by \code{as.matrix()}. The default is that all
#' \eqn{m}{m} predictors predict all \eqn{p}{p} effect sizes. The format of
#' this matrix follows \code{\link[metaSEM]{as.mxMatrix}}. The regression
#' coefficients can be constrained equally by using the same labels.
#' @param RE.constraints A \eqn{p}{p} x \eqn{p}{p} matrix specifying the
#' variance components of the random effects. If the input is not a matrix, it
#' is converted into a matrix by \code{as.matrix()}. The default is that all
#' covariance/variance components are free. The format of this matrix follows
#' \code{\link[metaSEM]{as.mxMatrix}}. Elements of the variance components can
#' be constrained equally by using the same labels. If a zero matrix is
#' specified, it becomes a fixed-effects meta-analysis.
#' @param RE.startvalues A vector of \eqn{p}{p} starting values on the
#' diagonals of the variance component of the random effects. If only one
#' scalar is given, it will be duplicated across the diagonals. Starting values
#' for the off-diagonals of the variance component are all 0. A \eqn{p}{p} x
#' \eqn{p}{p} symmetric matrix of starting values is also accepted.
#' @param RE.lbound A vector of \eqn{p}{p} lower bounds on the diagonals of the
#' variance component of the random effects. If only one scalar is given, it
#' will be duplicated across the diagonals. Lower bounds for the off-diagonals
#' of the variance component are set at \code{NA}. A \eqn{p}{p} x \eqn{p}{p}
#' symmetric matrix of the lower bounds is also accepted.
#' @param intervals.type Either \code{z} (default if missing) or \code{LB}. If
#' it is \code{z}, it calculates the 95\% Wald confidence intervals (CIs) based
#' on the z statistic. If it is \code{LB}, it calculates the 95\%
#' likelihood-based CIs on the parameter estimates. Note that the z values and
#' their associated p values are based on the z statistic. They are not related
#' to the likelihood-based CIs.
#' @param I2 Possible options are \code{"I2q"}, \code{"I2hm"} and
#' \code{"I2am"}. They represent the \code{I2} calculated by using a typical
#' within-study sampling variance from the Q statistic, the harmonic mean and
#' the arithmetic mean of the within-study sampling variances (Xiong, Miller, &
#' Morris, 2010). More than one options are possible. If
#' \code{intervals.type="LB"}, 95\% confidence intervals on the heterogeneity
#' indices will be constructed.
#' @param R2 Logical. If \code{TRUE} and there are predictors, R2 is calculated
#' (Raudenbush, 2009).
#' @param model.name A string for the model name in
#' \code{\link[OpenMx]{mxModel}}.
#' @param suppressWarnings Logical. If \code{TRUE}, warnings are suppressed.
#' The argument to be passed to \code{\link[OpenMx]{mxRun}}.
#' @param silent Logical. An argument to be passed to
#' \code{\link[OpenMx]{mxRun}}
#' @param run Logical. If \code{FALSE}, only return the mx model without
#' running the analysis.
#' @param \dots Further arguments to be passed to \code{\link[OpenMx]{mxRun}}
#' @return An object of class \code{meta} with a list of \item{call}{Object
#' returned by \code{\link[base]{match.call}}} \item{data}{A data matrix of y,
#' v and x } \item{no.y}{No. of effect sizes} \item{no.x}{No. of predictors}
#' \item{miss.x}{A vector indicating whether the predictors are missing.
#' Studies will be removed before the analysis if they are \code{TRUE}}
#' \item{I2}{Types of I2 calculated} \item{R2}{Logical} \item{mx.fit}{A fitted
#' object returned from \code{\link[OpenMx]{mxRun}}} \item{mx0.fit}{A fitted
#' object without any predictor returned from \code{\link[OpenMx]{mxRun}}}
#' @note Missing values (NA) in y and their related elements in v will be
#' removed automatically. When there are missing values in v but not in y,
#' missing values will be replaced by 1e5. Effectively, these effect sizes will
#' have little impact on the analysis. \code{metaFIML()} uses FIML to handle
#' missing covariates in X. It is experimental. It may not be stable.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{reml}}, \code{\link[metaSEM]{Hox02}},
#' \code{\link[metaSEM]{Berkey98}}, \code{\link[metaSEM]{wvs94a}}
#' @references Cheung, M. W.-L. (2008). A model for integrating fixed-,
#' random-, and mixed-effects meta-analyses into structural equation modeling.
#' \emph{Psychological Methods}, \bold{13}, 182-202.
#' 
#' Cheung, M. W.-L. (2009). Constructing approximate confidence intervals for
#' parameters with structural equation models. \emph{Structural Equation
#' Modeling}, \bold{16}, 267-294.
#' 
#' Cheung, M. W.-L. (2013). Multivariate meta-analysis as structural equation
#' models. \emph{Structural Equation Modeling}, \bold{20}, 429-454.
#' 
#' Cheung, M. W.-L. (2015). \emph{Meta-analysis: A structural equation modeling
#' approach}. Chichester, West Sussex: John Wiley & Sons, Inc.
#' 
#' Hardy, R. J., & Thompson, S. G. (1996). A likelihood approach to
#' meta-analysis with random effects. \emph{Statistics in Medicine}, \bold{15},
#' 619-629.
#' 
#' Neale, M. C., & Miller, M. B. (1997). The use of likelihood-based confidence
#' intervals in genetic models. \emph{Behavior Genetics}, \bold{27}, 113-120.
#' 
#' Raudenbush, S. W. (2009). Analyzing effect sizes: random effects models. In
#' H. M. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of
#' research synthesis and meta-analysis} (2nd ed., pp. 295-315). New York:
#' Russell Sage Foundation.
#' 
#' Xiong, C., Miller, J. P., & Morris, J. C. (2010). Measuring study-specific
#' heterogeneity in meta-analysis: application to an antecedent biomarker study
#' of Alzheimer's disease. \emph{Statistics in Biopharmaceutical Research},
#' \bold{2(3)}, 300-309. doi:10.1198/sbr.2009.0067
#' @keywords meta-analysis
meta <- function(y, v, x, data, intercept.constraints=NULL, coef.constraints=NULL,
                 RE.constraints=NULL, RE.startvalues=0.1, RE.lbound=1e-10,
                 intervals.type=c("z", "LB"), I2="I2q", R2=TRUE,
                 model.name="Meta analysis with ML",
                 suppressWarnings=TRUE, silent=TRUE, run=TRUE, ...) {
  mf <- match.call()
  if (missing(data)) {
    data <- sys.frame(sys.parent())
  } else {
    if (!is.data.frame(data)) data <- data.frame(data)
  }
  my.y <- mf[[match("y", names(mf))]]
  my.v <- mf[[match("v", names(mf))]]    
  y <- eval(my.y, data, enclos = sys.frame(sys.parent()))
  v <- eval(my.v, data, enclos = sys.frame(sys.parent()))
   
  if (is.vector(y)) no.y <- 1 else no.y <- ncol(y)  
  if (is.vector(v)) no.v <- 1 else no.v <- ncol(v)
  if (missing(x)) no.x <- 0 else {
    my.x <- mf[[match("x", names(mf))]]
    x <- eval(my.x, data, enclos = sys.frame(sys.parent()))
    if (is.vector(x)) no.x <- 1 else no.x <- ncol(x)
  }
  
  if ( no.v != no.y*(no.y+1)/2 )
    stop(paste("The expected no. of columns in v is ", no.y*(no.y+1)/2,
               " while the observed no. of columns in v is ", no.v, ".", sep=""))
  v.labels <- vech(outer(1:no.y, 1:no.y, function(x, y) paste("v", x,"_", y, sep = "")))
  y.labels <- paste("y", 1:no.y, sep="")
  x.labels <- paste("x", 1:no.x, sep="")

  ## If is.na(v), convert y into NA. NA in y will be handled automatically.
  ## Since NA in v (definition variable) is not allowed. Convert v into 1e10.
  ## Select variances only
  ## FIXME: how about NA in covariances?
  if (no.y==1) {
    y[is.na(v)] <- NA
  } else {
    index <- matrix(0, nrow=no.y, ncol=no.y)
    index[lower.tri(index, diag=TRUE)] <- seq(1, no.y*(no.y+1)/2)
    index <- Diag(index)
    y[is.na(v[, index])] <- NA
  }
    v[is.na(v)] <- 1e10
  
  ## FIXME: It is better to modify miss.x that includes regression coefficients
  if (no.x==0) {
    ## x <- NULL
    input.df <- as.matrix(cbind(y, v))
    dimnames(input.df) <- list(NULL, c(y.labels, v.labels))
    # No missing value in x
    miss.x <- rep(FALSE, nrow(input.df))
  } else {    
    input.df <- as.matrix(cbind(y, v, x))
    dimnames(input.df) <- list(NULL, c(y.labels, v.labels, x.labels))
    if (no.x==1) miss.x <- is.na(x) else miss.x <- apply(is.na(x), 1, any)
  }
  ## Remove missing data; my.df is used in the actual data analysis
  ## Missing y is automatically handled by OpenMx
  my.df <- input.df[!miss.x, ]

  ## Fix a bug reported by Noel Card that the number of statistics are incorrect and with negative dfs.
  ## It is due to my.df is a matrix, whereas OpenMx expects a data frame. ???
  my.df <- as.data.frame(my.df)
  
  ## Preparing the Beta1 matrix for the intercept vector
  ## Inter is a 1 by no.y row vector
  if (is.null(intercept.constraints)) {
    Inter <- matrix( paste("0*Intercept", 1:no.y, sep=""), nrow=1, ncol=no.y )
  } else {
    ## Convert intercept.constraints into a row matrix if it is not a matrix
    if (!is.matrix(intercept.constraints))
      intercept.constraints <- t(as.matrix(intercept.constraints))
    
    if (!all(dim(intercept.constraints)==c(1, no.y)))
      stop("Dimensions of \"intercept.constraints\" are incorrect.")
    Inter <- intercept.constraints
  }  
  Inter <- as.mxMatrix(t(Inter), name="Inter")
 
  ## Without predictors
  ## X: a 1 by (1+no.x) row vector
  if (no.x==0) {
    X <- mxMatrix("Unit", nrow=1, ncol=1, name="X")
    ## No predictor
    Beta1 <- mxAlgebra(Inter, name="Beta1")
    ## Not used; just make sure Beta is present in mxModel()
    Beta <- mxMatrix("Zero", nrow=1, ncol=1, name="Beta")
  } else {
    
    if (is.null(coef.constraints)) {
      yVar <- paste("y", seq(1,no.y), sep="", collapse="+")
      xVar <- paste("x", seq(1,no.x), sep="", collapse="+")
      # Use lm() coefficients as starting values
      startValues <- tryCatch( eval(parse(text=paste("t(coefficients(lm(cbind(",
                                                     yVar, ")~", xVar,", data=my.df)))", sep=""))),
                               error = function(e) e )
      # If error, replace it with 0. Added a column of intercepts
      # Fixed a minor bug that no starting value on the last predictor
      # when intercept.constraints=0
      if ( inherits(startValues, "error") | !is.null(intercept.constraints) )
        startValues <- matrix(0, nrow=no.y, ncol=(no.x+1))
  
      A.labels <- outer(1:no.y, 1:no.x, function(y, x) paste("*Slope", y,"_", x, sep = ""))
      Beta <- matrix( paste(startValues[,-1], A.labels, sep=""), nrow=no.y, ncol=no.x )
    } else {
    ## Convert coef.constraints into a column matrix if it is not a matrix
    if (!is.matrix(coef.constraints))
      coef.constraints <- as.matrix(coef.constraints)
      
      coef.dim <- dim(coef.constraints)
      if (!coef.dim[1]==no.y | !(coef.dim[2] %in% c(no.x, no.x+no.y)))
          stop("Dimensions of \"coef.constraints\" are incorrect.")
       Beta <- coef.constraints
    }
    Beta <- as.mxMatrix(Beta)
    Beta1 <- mxAlgebra( cbind(Inter, Beta), name="Beta1")  
    ## X.matrix <- paste("mxMatrix(\"Full\", nrow=1, ncol=(1+no.x), free=FALSE, values=c(1,",
    ##                   paste("data.x",1:no.x,sep="", collapse=","), "), name=\"X\")", sep="")
    ## eval(parse(text = X.matrix))
    X <- mxMatrix("Full", nrow=1, ncol=(1+no.x), free=FALSE, values=c(1, rep(NA, no.x)),
                  labels=c(NA, paste("data.x",1:no.x,sep="")), name="X")
  }
  
  expMean <- mxAlgebra( X %*% t(Beta1), name="expMean")
  
  ## Fixed a bug in 0.5-0 that lbound is not added into Tau
  ## when RE.constraints is used.
  ## lbound in variance component of the random effects
  if (is.matrix(RE.lbound)) {
    if (!all(dim(RE.lbound)==c(no.y, no.y)))
      warning("Dimensions of \"RE.lbound\" are incorrect.")
      # FIXME: need to handle unequal dimensions better
      lbound <- RE.lbound
      ## lbound is a matrix
    } else {
      lbound <- matrix(NA, nrow=no.y, ncol=no.y)
      Diag(lbound) <- RE.lbound
      ## lbound is a matrix      
  }  
 
  ## Preparing the S matrix for covariance elements
  #  No predictor
  if (is.null(RE.constraints)) {
    # Better to use starting values based on diagonal matrix rather than the UMM
    if (is.matrix(RE.startvalues)) {
      if (!all(dim(RE.startvalues)==c(no.y, no.y)))
        warning("Dimensions of \"RE.startvalues\" are incorrect.")
      values <- vech(RE.startvalues)
    } else {
      values <- vech(Diag(x=RE.startvalues, nrow=no.y, ncol=no.y))
    }
    Tau.labels <- vech(outer(1:no.y, 1:no.y, function(x,y) { paste("Tau2_",x,"_",y,sep="")}))
    Tau <- mxMatrix("Symm", ncol=no.y, nrow=no.y, free=TRUE, labels=Tau.labels,
                    lbound=vech(lbound), values=values, name="Tau")      
  } else {
    ## Convert RE.constraints into a column matrix if it is not a matrix
    if (!is.matrix(RE.constraints))
      RE.constraints <- as.matrix(RE.constraints)
    
    if (!all(dim(RE.constraints)==c(no.y, no.y)))
      stop("Dimensions of \"RE.constraints\" are incorrect.")
    ## Fixed a bug that reads lbound improperly
    ## Since as.mxMatrix expects a full matrix, lbound=vech(lbound) is incorrect
    Tau <- as.mxMatrix(RE.constraints, lbound=c(lbound), name="Tau")
  }
  V <- mxMatrix("Symm", ncol=no.y, nrow=no.y, free=FALSE,
                 labels=paste("data.", v.labels, sep=""), name="V")
  expCov <- mxAlgebra(V+Tau, name="expCov")
                
  ## Assuming NA first
  mx0.fit <- NA
  if (no.x==0) {

    I2 <- match.arg(I2, c("I2q", "I2hm", "I2am"), several.ok=TRUE)
    ## Select variances and exclude covariances
    v_het <- input.df[, paste("v", 1:no.y, "_", 1:no.y, sep=""), drop=FALSE]

    ## Calculate I2
    ## Based on Higgins and Thompson (2002), Eq. 9
    sum.w <- apply(v_het, 2, function(x) sum(1/x))
    sum.w2 <- apply(v_het, 2, function(x) sum(1/x^2))
    ## NA in v has been replaced by 1e10
    no.studies <- apply(v_het, 2, function(x) sum(x<1e9))
    ## Typical V based on Q statistic
    qV <- matrix((no.studies-1)*sum.w/(sum.w^2-sum.w2), nrow=1)
    ## Typical V based on harmonic mean  
    hmV <- matrix(no.studies/sum.w, nrow=1)
    ## Typical V based on arithmatic mean
    amV <- apply(v_het, 2, function(x) mean(x[x<1e9]))
    amV <- matrix(amV, nrow=1)
    V_het <- rbind(qV, hmV, amV)

    ## Select the heter.indices
    ## Before selection: V_het is a c("I2q","I2hm","I2am") by c(y1, y2, y3) matrix
    ## After selecting: A column vector of I2q(y1, y2, y3), I2hm(y1, y2, y3), I2am(y1, y2, y3)
    V_het <- matrix( t( V_het[c("I2q","I2hm","I2am")%in%I2, ] ), ncol=1 )
    V_het <- as.mxMatrix(V_het)

    One <- mxMatrix("Unit", nrow=length(I2), ncol=1, name="One")
    Tau_het <- mxAlgebra( One %x% diag2vec(Tau), name="Tau_het")    
    I2_values <- mxAlgebra( Tau_het/(Tau_het+V_het), name="I2_values")

    ## Modified for OpenMx 2.0
    mx.model <- mxModel(model=model.name, mxData(observed=my.df, type="raw"),
                        mxExpectationNormal(covariance="expCov", means="expMean", dimnames=y.labels),
                        mxFitFunctionML(),
                        Inter, Beta, Beta1, expMean, X, expCov, Tau, V, One, V_het, Tau_het, I2_values,
                        mxCI(c("Tau","Inter","I2_values")))
  } else {
    ## no.x > 0
      
    ## Modified for OpenMx 2.0
    mx.model <- mxModel(model=model.name, mxData(observed=my.df, type="raw"),
                        mxExpectationNormal(covariance="expCov", means="expMean", dimnames=y.labels),
                        mxFitFunctionML(),
                        Inter, Beta, Beta1, expMean, X, expCov, Tau, V, mxCI(c("Tau","Inter","Beta")))

    ## Calculate R2
    if (R2) mx0.fit <- tryCatch( meta(y=y, v=v, data=my.df, model.name="No predictor",
                                      suppressWarnings=TRUE, silent=TRUE), error = function(e) e )
  }

  
  ## meta <- mxModel(model=model.name, mxData(observed=my.df, type="raw"),
  ##                 mxFIMLObjective( covariance="S", means="M", dimnames=y.labels),
  ##                 Beta1, M, X, S, Tau, V, mxCI(c("Tau","Beta1")))

  ## Return mx model without running the analysis
  if (run==FALSE) return(mx.model)
  
  intervals.type <- match.arg(intervals.type)
  # Default is z
  switch(intervals.type,
    z = mx.fit <- tryCatch( mxRun(mx.model, intervals=FALSE, suppressWarnings=suppressWarnings,
                                  silent=silent, ...), error = function(e) e ),
    LB = mx.fit <- tryCatch( mxRun(mx.model, intervals=TRUE, suppressWarnings=suppressWarnings,
                                   silent=silent, ...), error = function(e) e ) )
 
  if (inherits(mx.fit, "error")) {
    cat("Error in running mxModel:\n")
    warning(print(mx.fit))
    return(mx.fit)
  } else {
     out <- list(call=mf, data=input.df, no.y=no.y, no.x=no.x, miss.x=miss.x, mx.model=mx.model,
                 I2=I2, R2=R2, mx.fit=mx.fit, mx0.fit=mx0.fit, intervals.type=intervals.type)
     class(out) <- "meta"
  }
  
  return(out)
}
