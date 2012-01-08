meta3 <- function(y, v, cluster, x, data, intercept.constraints, coef.constraints, 
                  RE2.constraints, RE2.lbound=1e-10,
                  RE3.constraints, RE3.lbound=1e-10,
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
  my.cluster <- mf[[match("cluster", names(mf))]]  
  y <- eval(my.y, data, enclos = sys.frame(sys.parent()))
  v <- eval(my.v, data, enclos = sys.frame(sys.parent()))
  cluster <- eval(my.cluster, data, enclos = sys.frame(sys.parent()))

  ## maximum no. of data in level-2 unit
  k <- max(sapply( split(cluster, cluster), length))
  ## data in long format
  my.long <- data.frame(y, v, cluster)

  ## Add level-2 and level-3 predictors into the data set
  if (missing(x)) {
    no.x <- 0
    miss.x <- rep(FALSE, nrow(my.long))
    } else {
    my.x <- mf[[match("x", names(mf))]]
    x <- eval(my.x, data, enclos = sys.frame(sys.parent()))
    if (is.vector(x)) no.x <- 1 else no.x <- ncol(x)
    old.labels <- names(my.long)
    my.long <- data.frame(my.long, x)
    names(my.long) <- c(old.labels, paste("x_", 1:no.x, sep=""))
    if (no.x==1) miss.x <- is.na(x) else miss.x <- apply(is.na(x), 1, any)
  }
  ## Remove missing data. Missing y is automatically handled by OpenMx.
  my.long <- my.long[!miss.x, ]
  ## index to order cluster
  index <- order(my.long$cluster)  
  my.long <- my.long[index, ]
  
  ## c() is required to convert matrix to vector when the data are balanced.
  ## c() is not required when the data are unbalanced.
  my.long$time <- c(unlist(sapply(split(my.long$y, my.long$cluster), function(x) 1:length(x))))
  my.wide <- reshape(my.long, timevar="time", idvar=c("cluster"), direction="wide")

  ## Replace "." with "_" since OpenMx does not allow "." as variable names
  names(my.wide) <- sub("\\.", "_", names(my.wide))
  ## Replace NA with 0 in v as NA is not allowed in definition variables
  temp <- my.wide[, paste("v_", 1:k, sep="")]
  temp[is.na(temp)] <- 0
  my.wide[, paste("v_", 1:k, sep="")] <- temp
  ## Missing indicator on y
  miss.y <- is.na(my.wide[, paste("y_", 1:k, sep="")])
  
  ## Prepare matrices
  if (missing(intercept.constraints)) {
    inter <- mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=0, labels="Intercept", name="inter")
  } else {
    if (!all(dim(intercept.constraints)==c(1, 1)))
      stop("Dimensions of \"intercept.constraints\" are incorrect.\n")
    inter <- as.mxMatrix(intercept.constraints, name="inter")
  }     
  oneRow <- mxMatrix("Full", nrow=1, ncol=k, free=FALSE, values=1, name="oneRow")
  
  if ( no.x==0 ) {
    ## matrices with all zeroes; not actually used in the calculations
    coeff <- mxMatrix("Full", nrow=1, ncol=1, free=FALSE, values=0, name="coeff")
    mydata <- mxMatrix("Full", nrow=1, ncol=k, free=FALSE, values=0, name="mydata")
  } else {
    ## NA is not available for definition variable.
    ## Replace NA with 0. Since y is missing, 0 does not affect the results.
    for (i in 1:no.x) {
      temp <- my.wide[, paste("x", i, 1:k, sep="_")]
      temp[miss.y] <- 0
      my.wide[, paste("x", i, 1:k, sep="_")] <- temp
    }    
    mydata <- mxMatrix("Full", nrow=no.x, ncol=k, free=FALSE, name="mydata",
                      labels=c(outer(1:no.x, 1:k, function(x, y) paste("data.x_", x,"_", y, sep = ""))))

    if (missing(coef.constraints)) {
      coeff <- mxMatrix("Full", nrow=1, ncol=no.x, free=TRUE, values=0,
                        labels=paste("Slope_", 1:no.x, sep=""), name="coeff")
    } else {
      coeff <- as.mxMatrix(coef.constraints, name="coeff")
    }    
  }
  
  ## if ( no.x3==0 ) {
  ##   ## matrices with all zeroes; not actually used in the calculations
  ##   coeff3 <- mxMatrix("Full", nrow=1, ncol=1, free=FALSE, values=0, name="coeff3")
  ##   data3 <- mxMatrix("Full", nrow=1, ncol=k, free=FALSE, values=0, name="data3")
  ## } else {
  ##   data3 <- mxMatrix("Full", nrow=no.x3, ncol=k, free=FALSE, name="data3",
  ##                      labels=rep( paste("data.x3_", 1:no.x3, sep=""), k ))
    
  ##   if (missing(coeff3.constraints)) {
  ##     coeff3 <- mxMatrix("Full", nrow=1, ncol=no.x3, free=TRUE, values=0,
  ##                        labels=paste("Slope3_", 1:no.x3, sep=""), name="coeff3")
  ##   } else {
  ##     coeff3 <- coeff3.constraints
  ##   }    
  ## }

  if ( length(RE2.lbound) != 1 ) {
    warning("\"RE2.lbound\" should be a scalar.")
    RE2.lbound <- 1e-10
  }
  if ( length(RE3.lbound) != 1 ) {
    warning("\"RE3.lbound\" should be a scalar.")
    RE3.lbound <- 1e-10
  }

  if ( missing(RE2.constraints) ) {
    Tau2 <- mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=0.01, labels="Tau2_2",                   
                     lbound=RE2.lbound, name="Tau2")
  } else {
    Tau2 <- as.mxMatrix(RE2.constraints, name="Tau2", lbound=RE2.lbound)
  }

  if ( missing(RE3.constraints) ) {
    Tau3 <- mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=0.01, labels="Tau2_3",                   
                     lbound=RE3.lbound, name="Tau3")
  } else {
    Tau3 <- as.mxMatrix(RE3.constraints, name="Tau3", lbound=RE3.lbound)
  }  
  
  Id <- mxMatrix("Iden", nrow=k, ncol=k, name="Id")
  Ones <- mxMatrix("Full", nrow=k, ncol=k, free=FALSE, values=1, name="Ones") 
  ## conditional sampling variances
  V <- mxMatrix("Diag", nrow=k, ncol=k, free=FALSE, labels=paste("data.v_", 1:k, sep=""), name="V")
  ## expMean <- mxAlgebra( oneRow %x% inter + coeff2 %*% data2 + coeff3 %*% data3, name="expMean")
  expMean <- mxAlgebra( oneRow %x% inter + coeff %*% mydata, name="expMean")     
  expCov <- mxAlgebra( Ones %x% Tau3 + Id %x% Tau2 + V, name="expCov")

  meta3 <- mxModel(model=model.name, mxData(observed=my.wide[,-1], type="raw"), oneRow, Id, Ones,
                   inter, coeff, mydata, Tau2, Tau3, V, expMean, expCov,
                   mxFIMLObjective("expCov","expMean", dimnames=paste("y_", 1:k, sep="")),           
                   mxCI(c("inter","coeff","Tau2","Tau3")))
  
  intervals.type <- match.arg(intervals.type)
  # Default is z
  switch(intervals.type,
    z = meta.fit <- tryCatch( mxRun(meta3, intervals=FALSE,
                                    suppressWarnings = suppressWarnings, ...), error = function(e) e ),
    LB = meta.fit <- tryCatch( mxRun(meta3, intervals=TRUE,
                                     suppressWarnings = suppressWarnings, ...), error = function(e) e ) )
 
  if (inherits(meta.fit, "error")) {
    cat("Error in running the mxModel:\n")
    stop(print(meta.fit))
  }
  
  out <- list(call = mf, data.wide=my.wide, data=my.long, no.y=1, no.x=no.x, miss.x=rep(FALSE, nrow(my.long)), 
              meta.fit=meta.fit)
  class(out) <- "meta"
  return(out)
}
