meta <- function(y, v, x, intercept.constraints, coeff.constraints,
                 RE.constraints, RE.startvalues=0.1, RE.lbound=1e-10,
                 intervals.type=c("z", "LB"), model.name="Meta analysis with ML",
                 suppressWarnings = TRUE, ...) {
  if (is.vector(y)) no.y <- 1 else no.y <- ncol(y)  
  if (is.vector(v)) no.v <- 1 else no.v <- ncol(v)
  if (missing(x)) no.x <- 0 else {
    if (is.vector(x)) no.x <- 1 else no.x <- ncol(x)
  }
  
  if ( no.v != no.y*(no.y+1)/2 )
    stop(paste("The expected no. of columns in v is ", no.y*(no.y+1)/2,
               " while the observed no. of columns in v is ", no.v, ".", sep=""))
  v.labels <- vech(outer(1:no.y, 1:no.y, function(x, y) paste("v", x,"_", y, sep = "")))
  y.labels <- paste("y", 1:no.y, sep="")
  x.labels <- paste("x", 1:no.x, sep="")

  if (no.x==0) {
    x <- NULL
    input.df <- as.matrix(cbind(y, v))
    dimnames(input.df) <- list(NULL, c(y.labels, v.labels))
    # No missing value in x
    miss.x <- rep(FALSE, nrow(input.df))
  } else {    
    input.df <- as.matrix(cbind(y, v, x))
    dimnames(input.df) <- list(NULL, c(y.labels, v.labels, x.labels))
    if (no.x==1) miss.x <- is.na(x) else miss.x <- apply(is.na(x), 1, any)
  }
  # Remove missing data; my.df is used in the actual data analysis
  my.df <- input.df[!miss.x, ]
    
  if (no.x==0) {
    # No predictor
    A <- mxMatrix("Zero", nrow=no.y, ncol=no.y, name="A")
  } else {
    if (missing(coeff.constraints)) {
      yVar <- paste("y", seq(1,no.y), sep="", collapse="+")
      xVar <- paste("x", seq(1,no.x), sep="", collapse="+")
      # Use lm() coefficients as starting values
      startValues <- tryCatch( eval(parse(text=paste("t(coefficients(lm(cbind(",
                               yVar, ")~", xVar,", data=data.frame(my.df))))", sep=""))) )
      # If error, replace it with 0. Added a column of intercepts
      if (inherits(startValues, "error"))
        startValues <- matrix(0, nrow=no.y, ncol=(no.x+1))
      
      # A1 includes regresson coefficients among ys
      # The default is no prediction among ys
      A.labels <- outer(1:no.y, 1:no.x, function(y, x) paste("*Slope", y,"_", x, sep = ""))
      A1 <- cbind(matrix(0, nrow=no.y, ncol=no.y),
                  matrix( paste(startValues[,-1], A.labels, sep=""), nrow=no.y, ncol=no.x ))
    } else {
      coeff.dim <- dim(coeff.constraints)
      if (!coeff.dim[1]==no.y | !(coeff.dim[2] %in% c(no.x, no.x+no.y)))
          stop("The dimensions of \"coeff.constraints\" are incorrect.")
      # Regression among ys have not been defined yet
      if (no.x==ncol(coeff.constraints)) {
        A1 <- cbind( matrix(0, nrow=no.y, ncol=no.y), coeff.constraints)
      } else {
        A1 <- coeff.constraints
      }
    }
    A <- rbind( A1,
                matrix(0, nrow=no.x, ncol=(no.x+no.y)))
    A <- as.mxMatrix(A, name="A")
  }

  ## Fixed a bug in 0.5-0 that lbound is not added into Tau
  ## when RE.constraints is used.
  ## lbound in variance component of the random effects
  if (is.matrix(RE.lbound)) {
    if (!all(dim(RE.lbound)==c(no.y, no.y)))
      warning("The dimensions of \"RE.lbound\" are incorrect.")
      # FIXME: need to handle unequal dimensions better
      # lbound is a matrix
      ## lbound <- vech(RE.lbound)
    } else {
      lbound <- matrix(NA, nrow=no.y, ncol=no.y)
      diag(lbound) <- RE.lbound
      # lbound is a matrix      
      ## lbound <- vech(lbound)
  }  
  
  ## Preparing the S matrix for covariance elements
  #  No predictor
  if (missing(RE.constraints)) {
    # Better to use starting values based on diagonal matrix rather than the UMM
    if (is.matrix(RE.startvalues)) {
      if (!all(dim(RE.startvalues)==c(no.y, no.y)))
        warning("The dimensions of \"RE.startvalues\" are incorrect.")
      values <- vech(RE.startvalues)
    } else {
      values <- vech(diag(x=RE.startvalues, nrow=no.y, ncol=no.y))
    }

    Tau.labels <- vech(outer(1:no.y, 1:no.y, function(x,y) { paste("Tau2_",x,"_",y,sep="")}))
    Tau <- mxMatrix("Symm", ncol=no.y, nrow=no.y, free=TRUE, labels=Tau.labels,
                    lbound=vech(lbound), values=values, name="Tau")      
  } else {
    if (!all(dim(RE.constraints)==c(no.y, no.y)))
      stop("The dimensions of \"RE.constraints\" are incorrect.")
    
    Tau <- as.mxMatrix(RE.constraints, lbound=lbound, name="Tau")
  }
  V <- mxMatrix("Symm", ncol=no.y, nrow=no.y, free=FALSE,
                 labels=paste("data.", v.labels, sep=""), name="V")
   
  if (no.x == 0) {
    S <- mxAlgebra(V+Tau, name="S")
  } else {
    S1 <- mxAlgebra(V+Tau, name="S1")
    S2 <- mxMatrix("Full", nrow=no.y, ncol=no.x, free=FALSE, name="S2")
    S3 <- mxMatrix("Full", nrow=no.x, ncol=no.y, free=FALSE, name="S3")
    S4 <- mxMatrix("Symm", nrow=no.x, ncol=no.x, free=TRUE,
                   values=vech(var(x, use="complete.obs")), name="S4")
    S <- mxAlgebra(rbind(cbind(S1,S2),
                         cbind(S3,S4)), name="S")
  }

  ## Preparing the M matrix for the mean vector 
  if (missing(intercept.constraints)) {
    M <- matrix( paste("0*Intercept", 1:no.y, sep=""), nrow=1, ncol=no.y )
    #M <- matrix("0*", nrow=1, ncol=no.y)
  } else {
    if (!all(dim(intercept.constraints)==c(1, no.y)))
      stop("The dimensions of \"intercept.constraints\" are incorrect.")
    M <- intercept.constraints
  }
  
  if (no.x==0) {
    M <- as.mxMatrix(M, name="M")
  } else {
    M <- cbind(M, matrix(paste(round(colMeans(matrix(x, ncol=no.x), na.rm=TRUE),2),"*", sep="")))
    M <- as.mxMatrix(M, name="M")
  }
 
  if (no.x==0) {
    F <- mxMatrix("Iden", nrow=(no.y), ncol=(no.y), free=FALSE, name="F",
                dimnames=list( paste("y", 1:no.y, sep=""), paste("y", 1:no.y, sep="") ))
  } else {
    F <- mxMatrix("Iden", nrow=(no.y+no.x), ncol=(no.y+no.x), free=FALSE, name="F",
                dimnames=list( c(y.labels, x.labels), c(y.labels, x.labels)) )
  }
  
  if (no.x==0) {
    meta <- mxModel(model=model.name, mxData(observed=my.df, type="raw"),
                    mxRAMObjective("A", "S", "F", "M"), A, S, F, M, Tau, V, 
                    mxCI(c("Tau","M")))
  } else {
    meta <- mxModel(model=model.name, mxData(observed=my.df, type="raw"),
                    mxRAMObjective("A", "S", "F", "M"), A, S, F, M, Tau, V, 
                    S1, S2, S3, S4, mxCI(c("Tau","M", "A")))
  }

  intervals.type <- match.arg(intervals.type)
  # Default is z
  switch(intervals.type,
    z = meta.fit <- tryCatch( mxRun(meta, intervals=FALSE,
                                    suppressWarnings = suppressWarnings, ...), error = function(e) e ),
    LB = meta.fit <- tryCatch( mxRun(meta, intervals=TRUE,
                                     suppressWarnings = suppressWarnings, ...), error = function(e) e ) )
 
  if (inherits(meta.fit, "error")) {
    cat("Error in running the mxModel:\n")
    stop(print(meta.fit))
  }
  
  out <- list(call = match.call(), data=input.df, no.y=no.y, no.x=no.x,
              miss.x=miss.x, meta.fit=meta.fit)
  class(out) <- "meta"
  return(out)
}
