meta <- function(y, v, x, data, intercept.constraints, coef.constraints,
                 RE.constraints, RE.startvalues=0.1, RE.lbound=1e-10,
                 intervals.type=c("z", "LB"), model.name="Meta analysis with ML",
                 suppressWarnings = TRUE, ...) {
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
  
  ## Replace NA in v with 1e5
  ## Users should manually remove elements in y and v when there are NA in v only.
  v[is.na(v)] <- 1e5
  
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

  if (no.x==0) {
    ## x <- NULL
    input.df <- as.matrix(cbind(y, v))
    dimnames(input.df) <- list(NULL, c(y.labels, v.labels))
    # No missing value in x
    miss.x <- rep(FALSE, nrow(input.df))
  } else {    
    input.df <- as.matrix(cbind(y, v, x))
    dimnames(input.df) <- list(NULL, c(y.labels, v.labels, x.labels))
    if (no.x==1) {
      ## miss.x: any one in x is missing
      miss.x <- is.na(x)
    } else {
      miss.x <- apply(is.na(x), 1, any)
    }
  }
  ## Remove missing data; my.df is used in the actual data analysis
  ## Missing y is automatically handled by OpenMx
  my.df <- input.df[!miss.x, ]

  ## Preparing the Beta matrix for the intercept vector
  ## Beta is a 1 by no.y row vector
  if (missing(intercept.constraints)) {
    Beta <- matrix( paste("0*Intercept", 1:no.y, sep=""), nrow=1, ncol=no.y )
  } else {
    if (!all(dim(intercept.constraints)==c(1, no.y)))
      stop("Dimensions of \"intercept.constraints\" are incorrect.")
    Beta <- intercept.constraints
  }  
  Beta <- t(Beta)
 
  ## Without predictors
  ## X: a 1 by (1+no.x) row vector
  if (no.x==0) {
    X <- mxMatrix("Unit", nrow=1, ncol=1, name="X") 
  } else {
    if (missing(coef.constraints)) {
      yVar <- paste("y", seq(1,no.y), sep="", collapse="+")
      xVar <- paste("x", seq(1,no.x), sep="", collapse="+")
      # Use lm() coefficients as starting values
      startValues <- tryCatch( eval(parse(text=paste("t(coefficients(lm(cbind(",
                               yVar, ")~", xVar,", data=data.frame(my.df))))", sep=""))) )
      # If error, replace it with 0. Added a column of intercepts
      if (inherits(startValues, "error"))
        startValues <- matrix(0, nrow=no.y, ncol=(no.x+1))
      
      A.labels <- outer(1:no.y, 1:no.x, function(y, x) paste("*Slope", y,"_", x, sep = ""))
      A <- matrix( paste(startValues[,-1], A.labels, sep=""), nrow=no.y, ncol=no.x )
    } else {
      coef.dim <- dim(coef.constraints)
      if (!coef.dim[1]==no.y | !(coef.dim[2] %in% c(no.x, no.x+no.y)))
          stop("Dimensions of \"coef.constraints\" are incorrect.")
       A <- coef.constraints
    }
   ## 
    Beta <- cbind(Beta, A)
    ## X.matrix <- paste("mxMatrix(\"Full\", nrow=1, ncol=(1+no.x), free=FALSE, values=c(1,",
    ##                   paste("data.x",1:no.x,sep="", collapse=","), "), name=\"X\")", sep="")
    ## eval(parse(text = X.matrix))
    X <- mxMatrix("Full", nrow=1, ncol=(1+no.x), free=FALSE, values=c(1, rep(NA, no.x)),
                  labels=c(NA, paste("data.x",1:no.x,sep="")), name="X")
  }
  Beta <- as.mxMatrix(Beta)
  M <- mxAlgebra( X %*% t(Beta), name="M")
  
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
      diag(lbound) <- RE.lbound
      ## lbound is a matrix      
  }  
 
  ## Preparing the S matrix for covariance elements
  #  No predictor
  if (missing(RE.constraints)) {
    # Better to use starting values based on diagonal matrix rather than the UMM
    if (is.matrix(RE.startvalues)) {
      if (!all(dim(RE.startvalues)==c(no.y, no.y)))
        warning("Dimensions of \"RE.startvalues\" are incorrect.")
      values <- vech(RE.startvalues)
    } else {
      values <- vech(diag(x=RE.startvalues, nrow=no.y, ncol=no.y))
    }
    Tau.labels <- vech(outer(1:no.y, 1:no.y, function(x,y) { paste("Tau2_",x,"_",y,sep="")}))
    Tau <- mxMatrix("Symm", ncol=no.y, nrow=no.y, free=TRUE, labels=Tau.labels,
                    lbound=vech(lbound), values=values, name="Tau")      
  } else {
    if (!all(dim(RE.constraints)==c(no.y, no.y)))
      stop("Dimensions of \"RE.constraints\" are incorrect.")    
    Tau <- as.mxMatrix(RE.constraints, lbound=vech(lbound), name="Tau")
  }
  V <- mxMatrix("Symm", ncol=no.y, nrow=no.y, free=FALSE,
                 labels=paste("data.", v.labels, sep=""), name="V")
  S <- mxAlgebra(V+Tau, name="S")
  
  meta <- mxModel(model=model.name, mxData(observed=my.df, type="raw"),
                  mxFIMLObjective( covariance="S", means="M", dimnames=y.labels),
                  Beta, M, X, S, Tau, V, mxCI(c("Tau","Beta")))

  intervals.type <- match.arg(intervals.type)
  # Default is z
  switch(intervals.type,
    z = meta.fit <- tryCatch( mxRun(meta, intervals=FALSE,
                                    suppressWarnings = suppressWarnings, ...), error = function(e) e ),
    LB = meta.fit <- tryCatch( mxRun(meta, intervals=TRUE,
                                     suppressWarnings = suppressWarnings, ...), error = function(e) e ) )
 
  if (inherits(meta.fit, "error")) {
    cat("Error in running mxModel:\n")
    warning(print(meta.fit))
  }
  
  out <- list(call = mf, data=input.df, no.y=no.y, no.x=no.x,
              miss.x=miss.x, meta.fit=meta.fit)
  class(out) <- "meta"
  return(out)
}
