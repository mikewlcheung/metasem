## library(metaSEM)

## x <- c(-2,-1.64,-1.33,-0.7,0,0.45,1.2,1.64,2.32,2.9)
## y <- c(0.699369,0.700462,0.695354,1.03905,1.97389,2.41143,1.91091,0.919576,-0.730975,-1.42001)*10
## v <- rnorm(10, mean=5, sd=0.1)
## data <- data.frame(x=x, y=y, v=0)

## model <- "y ~ m*1
##           y ~~ s*y
##           m == p1*cos(p2*data.x) + p2*sin(p1*data.x)
##           s == exp(VarErr)
##           "
## RAM <- lavaan2RAM(model, obs.variables = "y")
## startvalues <- list(p1=1.88, p2=0.7, VarErr=-5)

## fit <- create.mxModel(RAM=RAM, data=data, startvalues=startvalues)

## model <- "y ~ m*1
##           y ~~ s*y
##           m == p1*cos(p2*data.x) + p2*sin(p1*data.x)
##           fn := p1*p2
##           "
## RAM <- lavaan2RAM(model, obs.variables = "y")
## startvalues <- list(p1=1.88, p2=0.7, VarErr=-5)

## fit <- create.mxModel(RAM=RAM, data=data, startvalues=startvalues)
## summary(fit)

##     m1 <- "fy =~ 1*r
##            r ~~ data.r_v*r
##            fx =~ 1*JP_alpha
##            JP_alpha ~~ 0*JP_alpha
##            fy ~ Slope1_1*fx
##            fy ~~ Tau2_1_1*fy
##            fx ~~ CovX1_X1*fx
##            fx ~ MeanX1*1
##            fy ~ Intercept1*1"

## RAM <- lavaan2RAM(m1, obs.variables = c("r", "JP_alpha"), std.lv=FALSE)
## data <- Jaramillo05
##     fit1b <- create.mxModel(RAM=RAM, data=Jaramillo05)

create.mxModel <- function(model.name="mxModel", RAM=NULL, data=NULL,
                           Cov=NULL, means=NULL, numObs,
                           intervals.type = c("z", "LB"), startvalues=NULL, 
                           mxModel.Args=NULL, run=TRUE, mxTryHard=FALSE,
                           silent=TRUE, ...) {

    intervals.type <- match.arg(intervals.type)
    switch(intervals.type,
           z = intervals <- FALSE,
           LB = intervals <- TRUE)
    
    Amatrix <- as.mxMatrix(RAM$A, name="Amatrix")
    Smatrix <- as.mxMatrix(RAM$S, name="Smatrix")
    Fmatrix <- as.mxMatrix(RAM$F, name="Fmatrix")
    Mmatrix <- as.mxMatrix(RAM$M, name="Mmatrix")
    
    ## Some basic checking in RAM
    checkRAM(Amatrix, Smatrix, cor.analysis=FALSE)

    ## Extract the dimnames from Fmatrix$values
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
    
    ## Collate the starting values from RAM and add them to startvalues
    para.labels <- c(Amatrix$labels[Amatrix$free], Smatrix$labels[Smatrix$free],
                     Mmatrix$labels[Mmatrix$free])
    para.values <- c(Amatrix$values[Amatrix$free], Smatrix$values[Smatrix$free],
                     Mmatrix$values[Mmatrix$free])
    names(para.values) <- para.labels
    para.values <- as.list(para.values)
    ## Remove starting values from para.values if they are overlapped with startvalues    
    para.values[names(para.values) %in% names(startvalues)] <- NULL
    startvalues <- c(startvalues, para.values)

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
    ## Need to replace parameters with mxalgebras, if any TRUE
    if (any(index)) {

        ## Extract constraints that needed to be replaced
        mxalgebras.const <- mxalgebras[index]

        for (i in seq_along(mxalgebras.const)) {
            form <- as.character(mxalgebras.const[[i]]$formula)

            ## Replace the A matrix
            if (any(grep(form[2], A))) {
                A[which(form[2]==A)] <- form[3]
            }

            ## Replace the S matrix
            if (any(grep(form[2], S))) {
                S[which(form[2]==S)] <- form[3]
            }

            ## Replace the M matrix
            if (any(grep(form[2], M))) {
                M[which(form[2]==M)] <- form[3]
            }
        }

        ## Remove the constraints so they won't be added again
        mxalgebras[index] <- NULL
    }

    ## Check whether there are replacements
    ## Remove the starting values before comparisons
    if (all(A==as.symMatrix(RAM$A))) {
        mx.model <- mxModel(mx.model, Amatrix)
    } else {
        A <- as.mxAlgebra(A, startvalues=startvalues, name="Amatrix")
        mx.model <- mxModel(mx.model, A$mxalgebra, A$parameters, A$list)
    }

    if (all(S==as.symMatrix(RAM$S))) {
        mx.model <- mxModel(mx.model, Smatrix)
    } else {
        S <- as.mxAlgebra(S, startvalues=startvalues, name="Smatrix")
        mx.model <- mxModel(mx.model, S$mxalgebra, S$parameters, S$list)
    }
         
    ## Create an identity matrix from the no. of columens of Fmatrix,
    ## including all latent and observed variables
    Id <- as.mxMatrix(diag(ncol(Fmatrix$values)), name="Id")

    ## Expected covariance matrix and means of the observed and latent variables
    Id_A <- mxAlgebra(solve(Id - Amatrix), name="Id_A")
    expCov <- mxAlgebra(Id_A %&% Smatrix, name="expCov")   
    
    ## Add the mean structure only if there are means
    if (!is.null(data) | !is.null(means)) {

        if (all(M==as.symMatrix(RAM$M))) {
            mx.model <- mxModel(mx.model, Mmatrix)
        } else {
            M <- as.mxAlgebra(M, startvalues=startvalues, name="Mmatrix")
            mx.model <- mxModel(mx.model, M$mxalgebra, M$parameters, M$list)
        }
        
        expMean <- mxAlgebra(Mmatrix %*% t(Id_A), name="expMean")
        mx.model <- mxModel(mx.model, Fmatrix, Id, Id_A, expCov, expMean,
                            mxCI(c("Amatrix", "Smatrix", "Mmatrix")))
    } else {
        ## No mean structure
        mx.model <- mxModel(mx.model, Fmatrix, Id, Id_A, expCov, 
                            mxCI(c("Amatrix", "Smatrix")))
    }
    
    ## Add additional arguments to mxModel
    if (!is.null(mxModel.Args)) {
        for (i in seq_along(mxModel.Args)) {
            mx.model <- mxModel(mx.model, mxModel.Args[[i]])
        }
    }
    
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
                }
            }
    }
    
    if (run==TRUE) {
        if (mxTryHard==TRUE) {
            out <- mxTryHard(mx.model, intervals=intervals, silent=silent, ...)
        } else {
            out <- mxRun(mx.model, intervals=intervals, silent=silent, ...)
        }
    } else {
        out <- mx.model
    }
    out
}
