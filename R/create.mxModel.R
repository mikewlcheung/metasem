## Easy creation of mx models
## Note. Dimension names of variables are assumed stored in RAM$F;
## otherwise, they have to be specified in var.names. 
create.mxModel <- function(model.name="mxModel", RAM=NULL, Amatrix=NULL,
                           Smatrix=NULL, Fmatrix=NULL, Mmatrix=NULL,
                           Vmatrix=NULL, data, intervals.type = c("z", "LB"),
                           mx.algebras=NULL, mxModel.Args=NULL,
                           mxRun.Args=NULL, var.names=NULL,
                           suppressWarnings=TRUE,
                           silent=TRUE, run=TRUE, ...) {

    intervals.type <- match.arg(intervals.type)
    switch(intervals.type,
           z = intervals <- FALSE,
           LB = intervals <- TRUE)

    ## Read RAM first. If it is not specified, read individual matrices
    if (!is.null(RAM)) {
        Amatrix <- as.mxMatrix(RAM$A, name="Amatrix")
        Smatrix <- as.mxMatrix(RAM$S, name="Smatrix")
        Fmatrix <- as.mxMatrix(RAM$F, name="Fmatrix")
        Mmatrix <- as.mxMatrix(RAM$M, name="Mmatrix")
    } else {
        if (is.matrix(Amatrix)) {
            Amatrix <- as.mxMatrix(Amatrix, name="Amatrix")
        } else {    
            Amatrix@name <- "Amatrix"
        }
        if (is.matrix(Smatrix)) {
            Smatrix <- as.mxMatrix(Smatrix, name="Smatrix")
        } else {    
            Smatrix@name <- "Smatrix"
        }
        if (is.matrix(Fmatrix)) {
            Fmatrix <- as.mxMatrix(Fmatrix, name="Fmatrix")
        } else {    
            Fmatrix@name <- "Fmatrix"
        }
        if (is.matrix(Mmatrix)) {
            Mmatrix <- as.mxMatrix(Mmatrix, name="Mmatrix")
        } else {    
            Mmatrix@name <- "Mmatrix"
        }        
    }

    ## Some basic checking in RAM
    checkRAM(Amatrix, Smatrix, cor.analysis=FALSE)

    ## Extract the dimnames from Fmatrix$values
    if (is.null(var.names)) {
        var.names <- colnames(Fmatrix$values)
    }

    ## Add Vmatrix if provided
    if (!is.null(Vmatrix)) {
        if (is.matrix(Vmatrix)) {
            Vmatrix <- as.mxMatrix(Vmatrix)
            Vmatrix@name <- "Vmatrix"
        }
        
        Sfull <- mxAlgebra(Smatrix+Vmatrix, name="Sfull")
    } else {
        Sfull <- mxAlgebra(Smatrix, name="Sfull")
    }

    ## If matrix or data.frame is provided, setup mxData
    ## Otherwise, users have to setup mxData()
    if (is.data.frame(data) | is.matrix(data)) {
        mx.data <- mxData(observed=data, type="raw")
    } else {
        mx.data <- data
    }    

    ## Create an identity matrix from the dimensions of Amatrix
    Id <- as.mxMatrix(diag(nrow(Amatrix$values)), name="Id")

    ## Expected covariance matrix and means of the observed and latent variables
    Id_A <- mxAlgebra(solve(Id - Amatrix), name="Id_A")
    expCov <- mxAlgebra(Id_A %&% Smatrix, name="expCov")
    expMean <- mxAlgebra(Mmatrix %*% t(Id_A), name="expMean")
    
    mx.model <- mxModel(model.name, Amatrix, Smatrix, Fmatrix, Mmatrix,
                        Vmatrix, Sfull, Id, Id_A, expCov, expMean,
                        mx.data, mxFitFunctionML(),
                        mxCI(c("Amatrix", "Smatrix", "Mmatrix")),
                        mxExpectationRAM(A="Amatrix", S="Sfull", F="Fmatrix",
                                         M="Mmatrix",
                                         dimnames=var.names))
    
    ## Add additional arguments to mxModel
    if (!is.null(mxModel.Args)) {
        for (i in seq_along(mxModel.Args)) {
            mx.model <- mxModel(mx.model, mxModel.Args[[i]])
        }
    }

    ## Add additional mxAlgebras
    if (!is.null(mx.algebras)) {
        for (i in seq_along(mx.algebras)) {
            mx.model <- mxModel(mx.model, mx.algebras[[i]])
        }
        mx.model <- mxModel(mx.model,
                            mxCI(c("Amatrix", "Smatrix", "Mmatrix",
                                   names(mx.algebras))))
    }
    
    if (run==TRUE) {
        out <- do.call(mxRun,
                       c(list(mx.model, intervals=intervals,
                              suppressWarnings=suppressWarnings,
                              silent=silent), mxRun.Args))
    } else {
        out <- mx.model
    }
    out
}
