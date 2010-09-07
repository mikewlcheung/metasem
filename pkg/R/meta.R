meta <- function(y, v, x, intercept.constraints, coeff.constraints,
                 RE.constraints, RE.lbound=1e-10, intervals=FALSE, ...) {
  if (is.vector(y)) no.y <- 1 else no.y <- ncol(y)  
  if (is.vector(v)) no.v <- 1 else no.v <- ncol(v)
  if (missing(x)) no.x <- 0 else {
    if (is.vector(x)) no.x <- 1 else no.x <- ncol(x)
  }
  
  if ( no.v != no.y*(no.y+1)/2 )
    stop(paste("The expected no. of columns in v is ", no.y*(no.y+1)/2,
               " while the observed no. of columns in v is ", no.v, ".", sep=""))
  v.labels <- vech(outer(1:no.y, 1:no.y, function(x, y) paste("v", x, y, sep = "")))
  y.labels <- paste("y", 1:no.y, sep="")
  x.labels <- paste("x", 1:no.x, sep="")

  if (no.x==0) {
    x <- NULL
    my.df <- as.matrix(cbind(y, v))
    dimnames(my.df) <- list(NULL, c(y.labels, v.labels))
    # No missing value in x
    miss.x <- rep(FALSE, nrow(my.df))
  } else {
    if (no.x==1) miss.x <- is.na(x) else miss.x <- apply(is.na(x), 1, any)
    my.df <- as.matrix(cbind(y, v, x))
    # Remove missing data
    my.df <- my.df[!miss.x, ] 
    dimnames(my.df) <- list(NULL, c(y.labels, v.labels, x.labels))
  }
  
  if (no.x==0) {
    # No predictor
    A <- mxMatrix("Zero", nrow=no.y, ncol=no.y, name="A")
  } else {
    if (missing(coeff.constraints)) {
      yVar <- paste("y", seq(1,no.y), sep="", collapse="+")
      xVar <- paste("x", seq(1,no.x), sep="", collapse="+")
      # Use lm() coefficients as starting values
      startValues <- eval(parse(text=paste("t(coefficients(lm(cbind(",
                                  yVar, ")~", xVar,", data=data.frame(my.df))))", sep="")))
      # A1 includes regresson coefficients among ys
      # The default is no prediction among ys
      A1 <- cbind(matrix(0, nrow=no.y, ncol=no.y),
                  matrix(paste(round(startValues[,-1],2), "*", sep=""), nrow=no.y, ncol=no.x))
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
    
  
  ## Preparing the S matrix for covariance elements
  #  No predictor
  if (missing(RE.constraints)) {
    # Simplify some calculations with only 1 y
    if (no.y==1) {
      # Unweighted method of moments as starting value
      values <- max(RE.lbound, (var(y, na.rm=TRUE)-mean(v, na.rm=TRUE)))
      Tau <- mxMatrix("Symm", ncol=no.y, nrow=no.y, free=TRUE, lbound=RE.lbound,
                      values=values, name="Tau")
    } else {
      lbound <- vech(diag(x=RE.lbound, nrow=no.y, ncol=no.y))
      values <- matrix(0, nrow=no.y, ncol=no.y)
      values[lower.tri(values, diag=TRUE)] <- apply(v[!miss.x, ], 2, mean, na.rm=TRUE)
      values <- var(y[!miss.x, ], na.rm=TRUE) - values+t(values)-diag(diag(values))
      if (is.pd(values)) {
        values <- vech(values)
      } else {
        values <- vech(diag(x=0.0001, nrow=no.y, ncol=no.y))
      }
       Tau <- mxMatrix("Symm", ncol=no.y, nrow=no.y, free=TRUE, lbound=lbound, values=values, name="Tau")      
    }

  } else {
    if (!all(dim(RE.constraints)==c(no.y, no.y)))
      stop("The dimensions of \"RE.constraints\" are incorrect.")
    
    Tau <- as.mxMatrix(RE.constraints, name="Tau")
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
    M <- matrix("0*", nrow=1, ncol=no.y)
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
    meta <- mxModel("Meta analysis", mxData(observed=my.df, type="raw"),
                    mxRAMObjective("A", "S", "F", "M"), A, S, F, M, Tau, V, 
                    mxCI(c("Tau","M")))
  } else {
    meta <- mxModel("Meta analysis", mxData(observed=my.df, type="raw"),
                    mxRAMObjective("A", "S", "F", "M"), A, S, F, M, Tau, V, 
                    S1, S2, S3, S4, mxCI(c("Tau","M", "A")))
  }
  meta.fit <- mxRun(meta, intervals=intervals, ...)
  ## out <- list(call = match.call(), y=y, v=v, x=x, miss.x=miss.x, meta.fit=meta.fit)
  ## class(out) <- "meta"
  ## return(out)
}
