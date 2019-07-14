metaFIML <- function(y, v, x, av, data, intercept.constraints=NULL,
                     coef.constraints=NULL, RE.constraints=NULL,
                     RE.startvalues=0.1, RE.lbound=1e-10,
                     intervals.type=c("z", "LB"), R2=TRUE,
                     model.name="Meta analysis with FIML",
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

    ## FIXME: Allow missing x with av
    if (missing(x)) {
        stop("\"x\" cannot be missing.\n")
    } else {
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
    l.labels <- paste("l", 1:no.x, sep="")
    all.labels <- c(y.labels,x.labels,l.labels)

    if (missing(av)) {
        no.av <- 0    
    } else {
        my.av <- mf[[match("av", names(mf))]]
        av <- eval(my.av, data, enclos = sys.frame(sys.parent()))
        if (is.vector(av)) no.av <- 1 else no.av <- ncol(av)
        av.labels <- paste("av", 1:no.av, sep="")
        all.labels <- c(all.labels, av.labels)
    }     

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

    if (no.av==0) {
        ## x <- NULL
        my.df <- as.matrix(cbind(y, v, x))
        dimnames(my.df) <- list(NULL, c(y.labels, v.labels, x.labels))
    } else {    
        my.df <- as.matrix(cbind(y, v, x, av))
        dimnames(my.df) <- list(NULL, c(y.labels, v.labels, x.labels, av.labels))
    }

    ## M: intercept of y
    if (is.null(intercept.constraints)) {
        My <- matrix( paste("0*Intercept", 1:no.y, sep=""), nrow=1, ncol=no.y )
    } else {
        if (!all(dim(intercept.constraints)==c(1, no.y)))
            stop("Dimensions of \"intercept.constraints\" are incorrect.")
        My <- as.mxMatrix(t(intercept.constraints), name="My")
    }

    ## Intercept of x and av
    Mx <- matrix(0, nrow=1, ncol=no.x)
    Ml <- matrix( paste(round(apply(my.df[, -(1:(no.y+no.v)), drop=FALSE], 2, mean, na.rm=TRUE), 2),
                  paste("MeanX", 1:(no.x+no.av), sep=""), sep="*" ), nrow=1 )
    Mmatrix <- as.mxMatrix(cbind(My,Mx,Ml), name="Mmatrix")

    ## y1, y2, m1, m2, m3, l1, l2, l3, av1, av2              
    if (is.null(coef.constraints)) {
        yVar <- paste("y", seq(1,no.y), sep="", collapse="+")
        xVar <- paste("x", seq(1,no.x), sep="", collapse="+")
        ## Use lm() coefficients as starting values
        startValues <- tryCatch( eval(parse(text=paste("t(coefficients(lm(cbind(",
                                      yVar, ")~", xVar,", data=data.frame(my.df))))", sep=""))) )
        ## If error, replace it with 0. Added a column of intercepts
        if (inherits(startValues, "error") | !is.null(intercept.constraints))
            startValues <- matrix(0, nrow=no.y, ncol=(no.x+1))
      
        A.labels <- outer(1:no.y, 1:no.x, function(y, x) paste("*Slope", y,"_", x, sep = ""))
        A <- matrix( paste(startValues[,-1], A.labels, sep=""), nrow=no.y, ncol=no.x )
    } else {
        coef.dim <- dim(coef.constraints)
        if (!coef.dim[1]==no.y | !(coef.dim[2] %in% c(no.x, no.x+no.y)))
            stop("Dimensions of \"coef.constraints\" are incorrect.")
        A <- coef.constraints
    }
    
    ## cbind(y&m on y&m, rbind(A, l on m, l on l)) 
    Amatrix <- cbind( matrix(0, ncol=(no.y+no.x), nrow=(no.y+2*no.x)),
                      rbind(A, diag(no.x), matrix(0, ncol=no.x, nrow=no.x)) )
    ## There are av: rbind(cbind(A, l on all), all on av)
    if (no.av!=0) {
        Amatrix <- rbind( cbind(Amatrix, matrix(0, ncol=no.av, nrow=(no.y+2*no.x))),
                          matrix(0, nrow=no.av, ncol=(no.y+2*no.x+no.av)) )
    }
    Amatrix <- as.mxMatrix(Amatrix)

    ## starting values for both x and av
    CovX.st <- var(my.df[, -(1:(no.y+no.v))], na.rm=TRUE) 
    CovX.labels <- outer(1:(no.x+no.av), 1:(no.x+no.av), function(x, y) paste("CovX", x,"_X", y, sep = ""))
    CovX.labels[upper.tri(CovX.labels)] <- t(CovX.labels)[upper.tri(CovX.labels)]
    CovX <- matrix(paste(round(CovX.st,2), CovX.labels, sep="*"), ncol=(no.x+no.av), nrow=(no.x+no.av))
    CovX <- bdiagMat( list(matrix(0, ncol=no.x, nrow=no.x), CovX) )
    CovX <- as.mxMatrix(CovX)
  
    ## lbound in variance component of the random effects
    if (is.matrix(RE.lbound)) {
        if (!all(dim(RE.lbound)==c(no.y, no.y)))
            warning("Dimensions of \"RE.lbound\" are incorrect.")
        ## FIXME: need to handle unequal dimensions better
        lbound <- RE.lbound
        ## lbound is a matrix
    } else {
        lbound <- matrix(NA, nrow=no.y, ncol=no.y)
        Diag(lbound) <- RE.lbound
        ## lbound is a matrix      
    }  
 
    ## Preparing the S matrix for covariance elements
    ##  No predictor
    if (is.null(RE.constraints)) {
        ## Better to use starting values based on diagonal matrix rather than the UMM
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

        Tau <- as.mxMatrix(RE.constraints, lbound=c(lbound), name="Tau")
    }

    ## V known
    V <- mxMatrix("Symm", ncol=no.y, nrow=no.y, free=FALSE,
                  labels=paste("data.", v.labels, sep=""), name="V")
    Vy <- mxAlgebra(V+Tau, name="Vy")

    ## CovYX: X includes both x and av for convenience
    CovYX <- matrix(0, nrow=(2*no.x), ncol=no.y)
    if (no.av==0) {
        CovYX <- as.mxMatrix(CovYX)
    } else {
        ## av is not correlated with y
        Covyav <- outer(1:no.av, 1:no.y, function(x, y) paste("0*CovX", x+no.x,"_Y", y, sep = ""))
        ##Covyav <- matrix(0, nrow=no.av, ncol=no.y)
        CovYX <- as.mxMatrix(rbind(CovYX, Covyav), name="CovYX")
    }

    Smatrix <- mxAlgebra( cbind( rbind(Vy,CovYX), rbind(t(CovYX),CovX) ), name="Smatrix" )
    Fmatrix <- create.Fmatrix(c(rep(1, no.y+no.x), rep(0, no.x), rep(1, no.av)))
  
    mx.model <- mxModel(model=model.name, mxData(observed=my.df, type="raw"),
                        mxExpectationRAM(A="Amatrix", S="Smatrix", F="Fmatrix", M="Mmatrix", dimnames=all.labels),
                        mxFitFunctionML(),
                        Tau, CovX, CovYX, Amatrix, Smatrix, Fmatrix, Mmatrix, V, Vy,
                        mxCI(c("Amatrix","Mmatrix","Tau", "CovYX")))

    ## Calculate R2
    mx0.fit <- NA
    if (R2) mx0.fit <- tryCatch(metaFIML(y=y, v=v, x=x, av=av, data=my.df, model.name="No predictor",
                                         coef.constraints=matrix(0, nrow=no.y, ncol=no.x),
                                         suppressWarnings=TRUE, R2=FALSE, silent=TRUE),
                                error = function(e) e)    

    ## Return mx model without running the analysis
    if (run==FALSE) return(mx.model)
  
    intervals.type <- match.arg(intervals.type)

    ## Default is z
    switch(intervals.type,
           z = mx.fit <- tryCatch(mxRun(mx.model, intervals=FALSE,
                                  suppressWarnings = suppressWarnings, silent=silent, ...),
                                  error = function(e) e),
           LB = mx.fit <- tryCatch(mxRun(mx.model, intervals=TRUE,
                                  suppressWarnings = suppressWarnings, silent=silent, ...),
                                  error = function(e) e ))

    if (inherits(mx.fit, "error")) {
        cat("Error in running mxModel:\n")
        warning(print(mx.fit))
        return(mx.fit)
    } else {
        out <- list(call=mf, data=my.df, miss.x=rep(0, nrow(my.df)), no.y=no.y, no.x=no.x,
                    no.av=no.av, R2=R2, mx0.fit=mx0.fit, mx.fit=mx.fit, intervals.type=intervals.type)
        class(out) <- "meta"
        return(out)
    }

}

